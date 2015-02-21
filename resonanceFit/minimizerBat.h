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
//      wrapper around BAT
//
//-------------------------------------------------------------------------


#ifndef RESONANCEFIT_MINIMIZERBAT_HH
#define RESONANCEFIT_MINIMIZERBAT_HH

#include <string>
#include <vector>

#include <BAT/BCModel.h>

#include "forward.h"
#include "minimizer.h"

namespace rpwa {

	namespace resonanceFit {

		class minimizerBat : public rpwa::resonanceFit::minimizer {

		private:

			class functionAdaptor : public BCModel {

			public:

				functionAdaptor(const rpwa::resonanceFit::functionConstPtr& fitFunction);
				virtual ~functionAdaptor() {}

				double LogAPrioriProbability(const std::vector<double>& parameters);
				double LogLikelihood(const std::vector<double>& parameters);

			private:

				const rpwa::resonanceFit::functionConstPtr _fitFunction;

			};

		public:

			minimizerBat(const rpwa::resonanceFit::modelConstPtr& fitModel,
			             const rpwa::resonanceFit::functionConstPtr& fitFunction,
			             const std::string& outFileName);
			virtual ~minimizerBat() {}

			std::map<std::string, double> minimize(std::vector<std::string>& freeParameters,
			                                       rpwa::resonanceFit::parameters& fitParameters,
			                                       rpwa::resonanceFit::parameters& fitParametersError,
			                                       TMatrixT<double>& covarianceMatrix,
			                                       rpwa::resonanceFit::cache& cache);

		private:

			bool initParameters(const rpwa::resonanceFit::parameters& fitParameters,
			                    const std::string& freeParameters);

			const rpwa::resonanceFit::modelConstPtr _fitModel;

			rpwa::resonanceFit::minimizerBat::functionAdaptor _functionAdaptor;

			const std::string _outFileName;

		};

	} // end namespace resonanceFit

} // end namespace rpwa


#endif // RESONANCEFIT_MINIMIZERBAT_HH
