#include "massDepFitLikeli.h"

#include "massDepFitModel.h"


rpwa::massDepFit::likelihood*
rpwa::massDepFit::likelihood::Clone() const {
	return new likelihood(*this);
}


unsigned int
rpwa::massDepFit::likelihood::NDim() const {
	return _compset->getNrParameters();
}


unsigned int
rpwa::massDepFit::likelihood::NDataPoints() const {
	// calculate data points:
	// * diagonal elements are real numbers
	// * non-diagonal elements are complex numbers
	// * remember (Re,Im) => factor 2
	// * diagonal elements are only checked once, of diagonal elements with
	//   the two different combinations (i,j) and (j,i)
	unsigned int nrPts(0);

	for(size_t idxWave=0; idxWave<_nrWaves; ++idxWave) {
		for(size_t jdxWave=0; jdxWave<_nrWaves; ++jdxWave) {
			nrPts += _wavePairMassBinLimits[idxWave][jdxWave].second - _wavePairMassBinLimits[idxWave][jdxWave].first + 1;
		}
	}

	nrPts *= _nrBins;

	return nrPts;
}


bool
rpwa::massDepFit::likelihood::init(rpwa::massDepFit::model* compset,
                                   const std::vector<double>& massBinCenters,
                                   const boost::multi_array<std::complex<double>, 3>& productionAmplitudes,
                                   const boost::multi_array<double, 6>& productionAmplitudesCovariance,
                                   const boost::multi_array<std::complex<double>, 4>& spinDensityMatrices,
                                   const boost::multi_array<double, 6>& spinDensityCovarianceMatrices,
                                   const boost::multi_array<std::pair<size_t, size_t>, 2>& wavePairMassBinLimits,
                                   bool useCovariance)
{
	_compset = compset;

	_massBinCenters = massBinCenters;

	_productionAmplitudes.resize(std::vector<size_t>(productionAmplitudes.shape(), productionAmplitudes.shape()+productionAmplitudes.num_dimensions()));
	_productionAmplitudes = productionAmplitudes;
	_productionAmplitudesCovariance.resize(std::vector<size_t>(productionAmplitudesCovariance.shape(), productionAmplitudesCovariance.shape()+productionAmplitudesCovariance.num_dimensions()));
	_productionAmplitudesCovariance = productionAmplitudesCovariance;

	_spinDensityMatrices.resize(std::vector<size_t>(spinDensityMatrices.shape(), spinDensityMatrices.shape()+spinDensityMatrices.num_dimensions()));
	_spinDensityMatrices = spinDensityMatrices;
	_spinDensityCovarianceMatrices.resize(std::vector<size_t>(spinDensityCovarianceMatrices.shape(), spinDensityCovarianceMatrices.shape()+spinDensityCovarianceMatrices.num_dimensions()));
	_spinDensityCovarianceMatrices = spinDensityCovarianceMatrices;

	_wavePairMassBinLimits.resize(std::vector<size_t>(wavePairMassBinLimits.shape(), wavePairMassBinLimits.shape()+wavePairMassBinLimits.num_dimensions()));
	_wavePairMassBinLimits = wavePairMassBinLimits;

	_useCovariance = useCovariance;

	_nrBins = _spinDensityMatrices.size();
	_nrMassBins = _massBinCenters.size();
	_nrWaves = _wavePairMassBinLimits.size();

	return true;
}


double
rpwa::massDepFit::likelihood::DoEval(const double* par) const {
	// set parameters for resonances, background and phase space
	_compset->setParameters(par);

	return DoEvalSpinDensityMatrix();
}


double
rpwa::massDepFit::likelihood::DoEvalSpinDensityMatrix() const {
	double chi2=0;

	// loop over bins
	for(unsigned idxBin=0; idxBin<_nrBins; ++idxBin) {
		// loop over mass-bins
		for(unsigned idxMass=0; idxMass<_nrMassBins; ++idxMass) {
			const double mass = _massBinCenters[idxMass];

			// sum over the contributions to chi2 -> rho_ij
			for(size_t idxWave=0; idxWave<_nrWaves; ++idxWave) {
				for(size_t jdxWave=idxWave; jdxWave<_nrWaves; ++jdxWave) {
					// check that this mass bin should be taken into account for this
					// combination of waves
					if(idxMass < _wavePairMassBinLimits[idxWave][jdxWave].first || idxMass > _wavePairMassBinLimits[idxWave][jdxWave].second) {
						continue;
					}

					// calculate target spin density matrix element
					const std::complex<double> rhoFit = _compset->spinDensityMatrix(idxWave, jdxWave, idxBin, mass, idxMass);

					const std::complex<double> rhoDiff = rhoFit - _spinDensityMatrices[idxBin][idxMass][idxWave][jdxWave];

					double dchi;
					if(idxWave==jdxWave) {
						dchi = norm(rhoDiff) / _spinDensityCovarianceMatrices[idxBin][idxMass][idxWave][jdxWave][0][0];
					} else {
						if(_useCovariance) {
							dchi  = rhoDiff.real()*rhoDiff.real() * _spinDensityCovarianceMatrices[idxBin][idxMass][idxWave][jdxWave][1][1];
							dchi -= rhoDiff.real()*rhoDiff.imag() * _spinDensityCovarianceMatrices[idxBin][idxMass][idxWave][jdxWave][0][1];
							dchi -= rhoDiff.real()*rhoDiff.imag() * _spinDensityCovarianceMatrices[idxBin][idxMass][idxWave][jdxWave][1][0];
							dchi += rhoDiff.imag()*rhoDiff.imag() * _spinDensityCovarianceMatrices[idxBin][idxMass][idxWave][jdxWave][0][0];

							dchi /= _spinDensityCovarianceMatrices[idxBin][idxMass][idxWave][jdxWave][0][0]*_spinDensityCovarianceMatrices[idxBin][idxMass][idxWave][jdxWave][1][1]
							        - _spinDensityCovarianceMatrices[idxBin][idxMass][idxWave][jdxWave][0][1]*_spinDensityCovarianceMatrices[idxBin][idxMass][idxWave][jdxWave][1][0];
						} else {
							dchi  = rhoDiff.real()*rhoDiff.real() / _spinDensityCovarianceMatrices[idxBin][idxMass][idxWave][jdxWave][0][0];
							dchi += rhoDiff.imag()*rhoDiff.imag() / _spinDensityCovarianceMatrices[idxBin][idxMass][idxWave][jdxWave][1][1];
						}
					}
					chi2 += dchi;
				} // end loop over jdxWave
			} // end loop over idxWave
		} // end loop over mass-bins
	} // end loop over bins

	return chi2;
}
