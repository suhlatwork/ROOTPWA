///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2015 Sebastian Uhl (TUM)
//
//    This file is part of ROOTPWA
//
//    ROOTPWA is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    ROOTPWA is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with ROOTPWA.  If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//-------------------------------------------------------------------------
//
// Description:
//      implementation of the wrapper around BAT
//
//-------------------------------------------------------------------------


#include "minimizerBat.h"

#include <boost/tokenizer.hpp>

#include <BAT/BCParameter.h>

#include <fileUtils.hpp>
#include <reportingUtils.hpp>

#include "components.h"
#include "fsmd.h"
#include "function.h"
#include "model.h"


rpwa::resonanceFit::minimizerBat::functionAdaptor::functionAdaptor(const rpwa::resonanceFit::functionConstPtr& fitFunction)
	: BCModel(),
	  _fitFunction(fitFunction)
{
}


double
rpwa::resonanceFit::minimizerBat::functionAdaptor::LogAPrioriProbability(const std::vector<double>& /*parameters*/)
{
	return 0.;
}


double
rpwa::resonanceFit::minimizerBat::functionAdaptor::LogLikelihood(const std::vector<double>& parameters)
{
	return _fitFunction->logLikelihood(parameters);
}


rpwa::resonanceFit::minimizerBat::minimizerBat(const rpwa::resonanceFit::modelConstPtr& fitModel,
                                               const rpwa::resonanceFit::functionConstPtr& fitFunction,
                                               const std::string& outFileName)
	: _fitModel(fitModel),
	  _functionAdaptor(fitFunction),
	  _outFileName(outFileName)
{
}


std::map<std::string, double>
rpwa::resonanceFit::minimizerBat::minimize(std::vector<std::string>& freeParameters,
                                           rpwa::resonanceFit::parameters& fitParameters,
                                           rpwa::resonanceFit::parameters& /*fitParametersError*/,
                                           TMatrixT<double>& /*covarianceMatrix*/,
                                           rpwa::resonanceFit::cache& cache)
{
	// in case the freeParameters vector is empty, free all parameters
	if(freeParameters.size() == 0) {
		printWarn << "using default release order of parameters." << std::endl;
		freeParameters.push_back("*");
	} else {
		// FIXME
	}

	if(not initParameters(fitParameters,
	                      freeParameters[0])) {
		printErr << "error while setting start parameters." << std::endl;
		return std::map<std::string, double>();
	}

	_functionAdaptor.MarginalizeAll();

	_functionAdaptor.FindMode(_functionAdaptor.GetBestFitParameters());

	std::string plotFileName(_outFileName);
	if(rpwa::extensionFromPath(plotFileName) == "root") {
		plotFileName = rpwa::changeFileExtension(plotFileName, "plot.pdf");
	} else {
		plotFileName += ".plot.pdf";
	}
	_functionAdaptor.PrintAllMarginalized(plotFileName.c_str());

	std::string updateFileName(_outFileName);
	if(rpwa::extensionFromPath(updateFileName) == "root") {
		updateFileName = rpwa::changeFileExtension(updateFileName, "update.pdf");
	} else {
		updateFileName += ".update.pdf";
	}
	_functionAdaptor.PrintKnowledgeUpdatePlots(updateFileName.c_str());

	_fitModel->importParameters(_functionAdaptor.GetBestFitParameters().data(), fitParameters, cache);

	// store and print information about fit quality
	std::map<std::string, double> fitQuality;

	return fitQuality;
}


bool
rpwa::resonanceFit::minimizerBat::initParameters(const rpwa::resonanceFit::parameters& fitParameters,
                                                 const std::string& freeParameters)
{
	// tokenize freeParameters string (default separators also include '*' and ',')
	boost::char_separator<char> separators(" \t\n");
	boost::tokenizer<boost::char_separator<char> > tokenizeFreeParameters(freeParameters, separators);

	// changes status of variables (fixed/released)
	// * couplings have to be freed explicitely by adding 'coupling' to freeParameters
	// * branchings also have to be freed explicitely with the keyword 'branching'
	// * additional parameters can be freed with freeParameters
	// * fixed values from config remain fixed

	size_t parcount=0;
	// first add all couplings
	for(size_t idxComponent = 0; idxComponent < _fitModel->getNrComponents(); ++idxComponent) {
		const rpwa::resonanceFit::componentConstPtr& comp = _fitModel->getComponent(idxComponent);
		for(size_t idxCoupling=0; idxCoupling<comp->getNrCouplings(); ++idxCoupling) {
			const rpwa::resonanceFit::component::channel& channel = comp->getChannelFromCouplingIdx(idxCoupling);
			const std::vector<size_t>& bins = channel.getBins();
			for(size_t i = 0; i < bins.size(); ++i) {
				const size_t idxBin = bins[i];
				std::ostringstream prefixBin;
				prefixBin << "coupling__bin"
				          << idxBin;
				std::ostringstream prefixName;
				prefixName << prefixBin.str()
				           << "__"
				           << comp->getName()
				           << "__";
				if(comp->getNrBranchings() > 1) {
					const std::string waveQN = channel.getWaveName().substr(0, channel.getWaveName().find("="));
					prefixName << waveQN;
				} else {
					prefixName << channel.getWaveName();
				}

				bool free = false;
				if(find(tokenizeFreeParameters.begin(), tokenizeFreeParameters.end(), "*") != tokenizeFreeParameters.end()
				   or find(tokenizeFreeParameters.begin(), tokenizeFreeParameters.end(), "coupling") != tokenizeFreeParameters.end()
				   or find(tokenizeFreeParameters.begin(), tokenizeFreeParameters.end(), prefixBin.str()) != tokenizeFreeParameters.end()
				   or find(tokenizeFreeParameters.begin(), tokenizeFreeParameters.end(), prefixName.str()) != tokenizeFreeParameters.end()) {
					free = true;
				}
				bool fix = not free;

				const std::complex<double> parameter = fitParameters.getCoupling(idxComponent, idxCoupling, idxBin);

				if (fix) {
					printInfo << "parameter " << parcount << " ('" << (prefixName.str() + "__real") << "') fixed to " << parameter.real() << std::endl;
					BCParameter batParameterReal(prefixName.str() + "__real",
					                             parameter.real(),
					                             parameter.real());
					batParameterReal.Fix(parameter.real());
					_functionAdaptor.AddParameter(batParameterReal);
					++parcount;

					if(not channel.isAnchor(idxBin)) {
						printInfo << "parameter " << parcount << " ('" << (prefixName.str() + "__imag") << "') fixed to " << parameter.imag() << std::endl;
						BCParameter batParameterImag(prefixName.str() + "__imag",
						                             parameter.imag(),
						                             parameter.imag());
						batParameterImag.Fix(parameter.imag());
						_functionAdaptor.AddParameter(batParameterImag);
						++parcount;
					}
				} else {
					printInfo << "parameter " << parcount << " ('" << (prefixName.str() + "__real") << "') is unlimited" << std::endl;
					_functionAdaptor.AddParameter(prefixName.str() + "__real",
					                              channel.isAnchor(idxBin) ? 0. : -std::numeric_limits<double>::max(),
					                              std::numeric_limits<double>::max());
					++parcount;

					if(not channel.isAnchor(idxBin)) {
						printInfo << "parameter " << parcount << " ('" << (prefixName.str() + "__imag") << "') is unlimited" << std::endl;
						_functionAdaptor.AddParameter(prefixName.str() + "__imag",
						                              -std::numeric_limits<double>::max(),
						                              std::numeric_limits<double>::max());
						++parcount;
					}
				}
			}
		} // end loop over channels
	} // end loop over components

	// second eventually add all branchings
	for(size_t idxComponent = 0; idxComponent < _fitModel->getNrComponents(); ++idxComponent) {
		const rpwa::resonanceFit::componentConstPtr& comp = _fitModel->getComponent(idxComponent);
		for(size_t idxBranching = 0; idxBranching < comp->getNrBranchings(); ++idxBranching) {
			// skip branchings that are always real and fixed to 1
			if(comp->isBranchingFixed(idxBranching)) {
				continue;
			}

			const rpwa::resonanceFit::component::channel& channel = comp->getChannelFromBranchingIdx(idxBranching);
			const std::string waveQN = channel.getWaveName().substr(0, channel.getWaveName().find("="));
			const std::string waveDecay = channel.getWaveName().substr(channel.getWaveName().find("=")+1);
			std::ostringstream prefixName;
			prefixName << "branching__"
			           << comp->getName()
			           << "__"
			           << waveQN
			           << "__"
			           << waveDecay;

			bool free = false;
			if(find(tokenizeFreeParameters.begin(), tokenizeFreeParameters.end(), "*") != tokenizeFreeParameters.end()
			   or find(tokenizeFreeParameters.begin(), tokenizeFreeParameters.end(), "branching") != tokenizeFreeParameters.end()
			   or find(tokenizeFreeParameters.begin(), tokenizeFreeParameters.end(), prefixName.str()) != tokenizeFreeParameters.end()) {
				free = true;
			}
			bool fix = not free;

			const std::complex<double> parameter = fitParameters.getBranching(idxComponent, idxBranching);

			if(fix) {
				printInfo << "parameter " << parcount << " ('" << (prefixName.str() + "__real") << "') fixed to " << parameter.real() << std::endl;
				BCParameter batParameterReal(prefixName.str() + "__real",
				                             parameter.real(),
				                             parameter.real());
				batParameterReal.Fix(parameter.real());
				_functionAdaptor.AddParameter(batParameterReal);
				++parcount;

				printInfo << "parameter " << parcount << " ('" << (prefixName.str() + "__imag") << "') fixed to " << parameter.imag() << std::endl;
				BCParameter batParameterImag(prefixName.str() + "__imag",
				                             parameter.imag(),
				                             parameter.imag());
				batParameterImag.Fix(parameter.imag());
				_functionAdaptor.AddParameter(batParameterImag);
				++parcount;
			} else {
				printInfo << "parameter " << parcount << " ('" << (prefixName.str() + "__real") << "') is unlimited" << std::endl;
				_functionAdaptor.AddParameter(prefixName.str() + "__real",
				                              -std::numeric_limits<double>::max(),
				                              std::numeric_limits<double>::max());
				++parcount;

				printInfo << "parameter " << parcount << " ('" << (prefixName.str() + "__imag") << "') is unlimited" << std::endl;
				_functionAdaptor.AddParameter(prefixName.str() + "__imag",
				                              -std::numeric_limits<double>::max(),
				                              std::numeric_limits<double>::max());
				++parcount;
			}
		} // end loop over channels
	} // end loop over components

	// third add parameters of the components, i.e. mass and width
	for(size_t idxComponent = 0; idxComponent < _fitModel->getNrComponents(); ++idxComponent) {
		const rpwa::resonanceFit::componentConstPtr& comp = _fitModel->getComponent(idxComponent);
		for(size_t idxParameter=0; idxParameter<comp->getNrParameters(); ++idxParameter) {
			const rpwa::resonanceFit::parameter& parameter = comp->getParameter(idxParameter);
			const std::string name = comp->getName() + "__" + parameter.name();

			bool free = false;
			if(find(tokenizeFreeParameters.begin(), tokenizeFreeParameters.end(), "*") != tokenizeFreeParameters.end()
			   or find(tokenizeFreeParameters.begin(), tokenizeFreeParameters.end(), comp->getName()) != tokenizeFreeParameters.end()
			   or find(tokenizeFreeParameters.begin(), tokenizeFreeParameters.end(), parameter.name()) != tokenizeFreeParameters.end()
			   or find(tokenizeFreeParameters.begin(), tokenizeFreeParameters.end(), name) != tokenizeFreeParameters.end()) {
				free = true;
			}
			bool fix = not free;
			if(parameter.fixed()) {
				fix = true;
			}

			const double startValue = fitParameters.getParameter(idxComponent, idxParameter);

			if(fix) {
				printInfo << "parameter " << parcount << " ('" << name << "') fixed to " << startValue << std::endl;
				BCParameter batParameter(name,
				                         startValue,
				                         startValue);
				batParameter.Fix(startValue);
				_functionAdaptor.AddParameter(batParameter);
			} else if(parameter.limitedLower() and parameter.limitedUpper()) {
				printInfo << "parameter " << parcount << " ('" << name << "') "
				          << "limited between " << parameter.limitLower()
				          << " and " << parameter.limitUpper() << std::endl;
				_functionAdaptor.AddParameter(name,
				                              parameter.limitLower(),
				                              parameter.limitUpper());
			} else if(parameter.limitedLower()) {
				printInfo << "parameter " << parcount << " ('" << name << "') "
				          << "limited larger than " << parameter.limitLower() << std::endl;
				_functionAdaptor.AddParameter(name,
				                              parameter.limitLower(),
				                              std::numeric_limits<double>::max());
			} else if(parameter.limitedUpper()) {
				printInfo << "parameter " << parcount << " ('" << name << "') "
				          << "limited smaller than " << parameter.limitUpper() << std::endl;
				_functionAdaptor.AddParameter(name,
				                              -std::numeric_limits<double>::max(),
				                              parameter.limitUpper());
			} else {
				printInfo << "parameter " << parcount << " ('" << name << "') is unlimited" << std::endl;
				_functionAdaptor.AddParameter(name,
				                              -std::numeric_limits<double>::max(),
				                              std::numeric_limits<double>::max());
			}
			++parcount;
		}
	} // end loop over components

	// set parameters for final-state mass-dependence
	if(_fitModel->getFsmd()) {
		const rpwa::resonanceFit::fsmdConstPtr& fsmd = _fitModel->getFsmd();
		const size_t maxNrBins = fsmd->isSameFunctionForAllBins() ? 1 : fsmd->getNrBins();
		for(size_t idxBin = 0; idxBin < maxNrBins; ++idxBin) {
			for(size_t idxParameter = 0; idxParameter < fsmd->getNrParameters(idxBin); ++idxParameter) {
				const rpwa::resonanceFit::parameter& parameter = fsmd->getParameter(idxBin, idxParameter);
				std::ostringstream name;
				name << "fsmd__bin"
				     << idxBin
				     << "__"
				     << parameter.name();

				const bool fix = parameter.fixed();

				const double startValue = fitParameters.getParameter(_fitModel->getNrComponents(), fsmd->getParameterIndex(idxBin)+idxParameter);

				if(fix) {
					printInfo << "parameter " << parcount << " ('" << name.str() << "') fixed to " << startValue << std::endl;
					BCParameter batParameter(name.str(),
					                         startValue,
					                         startValue);
					batParameter.Fix(startValue);
					_functionAdaptor.AddParameter(batParameter);
				} else if(parameter.limitedLower() and parameter.limitedUpper()) {
					printInfo << "parameter " << parcount << " ('" << name.str() << "') "
					          << "limited between " << parameter.limitLower()
					          << " and " << parameter.limitUpper() << std::endl;
					_functionAdaptor.AddParameter(name.str(),
					                              parameter.limitLower(),
					                              parameter.limitUpper());
				} else if(parameter.limitedLower()) {
					printInfo << "parameter " << parcount << " ('" << name.str() << "') "
					          << "limited larger than " << parameter.limitLower() << std::endl;
					_functionAdaptor.AddParameter(name.str(),
					                              parameter.limitLower(),
					                              std::numeric_limits<double>::max());
				} else if(parameter.limitedUpper()) {
					printInfo << "parameter " << parcount << " ('" << name.str() << "') "
					          << "limited smaller than " << parameter.limitUpper() << std::endl;
					_functionAdaptor.AddParameter(name.str(),
					                              -std::numeric_limits<double>::max(),
					                              parameter.limitUpper());
				} else {
					printInfo << "parameter " << parcount << " ('" << name.str() << "') is unlimited" << std::endl;
					_functionAdaptor.AddParameter(name.str(),
					                              -std::numeric_limits<double>::max(),
					                              std::numeric_limits<double>::max());
				}
				++parcount;
			}
		}
	}

	return true;
}
