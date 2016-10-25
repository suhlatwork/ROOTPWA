///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2010-2012 Sebastian Neubert (TUM)
//    Copyright 2014,2015 Sebastian Uhl (TUM)
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
//      fit model for resonance fit
//      * calculate production amplitudes for single waves
//      * calculate spin-density matrix elements
//
//-------------------------------------------------------------------------


#ifndef RESONANCEFIT_MODEL_HH
#define RESONANCEFIT_MODEL_HH

#include <complex>
#include <iostream>
#include <limits>
#include <vector>

#include <boost/multi_array.hpp>

#include "forward.h"

namespace rpwa {

	namespace resonanceFit {

		class cache;
		class parameters;

		class model {

		public:

			model();

			void add(const rpwa::resonanceFit::componentPtr& comp);

			bool init(const rpwa::resonanceFit::informationConstPtr& fitInformation,
			          const std::string& anchorWaveName,
			          const std::string& anchorComponentName);

			bool isMappingEqualInAllBins() const { return _mappingEqualInAllBins; }

			size_t getNrParameters() const { return _nrParameters; }
			void importParameters(const double* par,
			                      rpwa::resonanceFit::parameters& parameters,
			                      rpwa::resonanceFit::cache& cache) const;

			size_t getNrComponents() const { return _components.size(); }
			rpwa::resonanceFit::componentConstPtr getComponent(size_t idxComponent) const { return _components[idxComponent]; }

			rpwa::resonanceFit::fsmdConstPtr getFsmd() const { return _fsmd; }
			void setFsmd(const rpwa::resonanceFit::fsmdPtr& fsmd);

			size_t getMaxChannelsInComponent() const { return _maxChannelsInComponent; }
			size_t getMaxParametersInComponent() const { return _maxParametersInComponent; }

			size_t getAnchorWave() const { return _idxAnchorWave; }

			const std::vector<std::pair<size_t, size_t> >& getComponentChannel(const size_t idxBin, const size_t idxWave) const { return _waveComponentChannel[idxBin][idxWave]; }

			std::complex<double> productionAmplitude(const rpwa::resonanceFit::parameters& fitParameters,
			                                         rpwa::resonanceFit::cache& cache,
			                                         const size_t idxWave,
			                                         const size_t idxBin,
			                                         const double mass,
			                                         const size_t idxMass = std::numeric_limits<size_t>::max()) const;
			double intensity(const rpwa::resonanceFit::parameters& fitParameters,
			                 rpwa::resonanceFit::cache& cache,
			                 const size_t idxWave,
			                 const size_t idxBin,
			                 const double mass,
			                 const size_t idxMass = std::numeric_limits<size_t>::max()) const;
			double phaseAbsolute(const rpwa::resonanceFit::parameters& fitParameters,
			                     rpwa::resonanceFit::cache& cache,
			                     const size_t idxWave,
			                     const size_t idxBin,
			                     const double mass,
			                     const size_t idxMass = std::numeric_limits<size_t>::max()) const;
			std::complex<double> spinDensityMatrix(const rpwa::resonanceFit::parameters& fitParameters,
			                                       rpwa::resonanceFit::cache& cache,
			                                       const size_t idxWave,
			                                       const size_t jdxWave,
			                                       const size_t idxBin,
			                                       const double mass,
			                                       const size_t idxMass = std::numeric_limits<size_t>::max()) const;
			double phase(const rpwa::resonanceFit::parameters& fitParameters,
			             rpwa::resonanceFit::cache& cache,
			             const size_t idxWave,
			             const size_t jdxWave,
			             const size_t idxBin,
			             const double mass,
			             const size_t idxMass = std::numeric_limits<size_t>::max()) const;

			std::ostream& print(std::ostream& out = std::cout) const;

		private:

			bool initMapping(const rpwa::resonanceFit::informationConstPtr& fitInformation,
			                 const std::string& anchorWaveName,
			                 const std::string& anchorComponentName);

			bool _mappingEqualInAllBins;

			std::vector<rpwa::resonanceFit::componentPtr> _components;

			rpwa::resonanceFit::fsmdPtr _fsmd;

			size_t _nrParameters;
			size_t _maxChannelsInComponent;
			size_t _maxParametersInComponent;

			size_t _idxAnchorWave;
			size_t _idxAnchorComponent;
			std::vector<size_t> _idxAnchorChannel;

			boost::multi_array<std::vector<std::pair<size_t, size_t> >, 2> _waveComponentChannel;

		};

		std::ostream& operator<< (std::ostream& out, const rpwa::resonanceFit::model& fitModel);

	} // end namespace resonanceFit

} // end namespace rpwa


inline
std::ostream&
rpwa::resonanceFit::operator<< (std::ostream& out, const rpwa::resonanceFit::model& fitModel)
{
	return fitModel.print(out);
}


#endif // RESONANCEFIT_MODEL_HH
