#include "fitResult_py.h"

#include <boost/python.hpp>

#include <TTree.h>

#include "boostContainers_py.hpp"
#include "fitResult.h"
#include "rootConverters_py.h"
#include "stlContainers_py.h"

namespace bp = boost::python;


namespace {

	void fitResult_fill_1(rpwa::fitResult&           self,
	                      const unsigned int         nmbEvents,
	                      const unsigned int         normNmbEvents,
	                      const double               massBinCenter,
	                      const double               logLikelihood,
	                      const int                  rank,
	                      const bp::tuple&           pyWaveInfo,
	                      const bp::tuple&           pyProdAmpInfo,
	                      PyObject*                  pyFitParCovMatrix,
	                      const rpwa::complexMatrix& normIntegral,
	                      const rpwa::complexMatrix& acceptedNormIntegral,
	                      const bp::object&          pyPhaseSpaceIntegral,
	                      const bool                 converged,
	                      const bool                 hasHessian)
	{
		boost::tuples::tuple<bp::list, bp::list, bp::list> btWaveInfo;
		if(not rpwa::py::convertBPTupleToTuple<bp::list, bp::list>(pyWaveInfo, btWaveInfo)) {
			PyErr_SetString(PyExc_TypeError, "Got invalid input for prodAmpInfo when executing rpwa::fitResult::fill()");
			bp::throw_error_already_set();
		}
		rpwa::fitResult::waveInfoType waveInfo;
		const bp::list& pyWaveNames = boost::tuples::get<0>(btWaveInfo);
		std::vector<std::string>& waveNames = boost::tuples::get<0>(waveInfo);
		if(not rpwa::py::convertBPObjectToVector<std::string>(pyWaveNames, waveNames)) {
			PyErr_SetString(PyExc_TypeError, "Got invalid input for waveNames when executing rpwa::fitResult::fill()");
			bp::throw_error_already_set();
		}
		const bp::list& pyWaveRefls = boost::tuples::get<1>(btWaveInfo);
		std::vector<int>& waveRefls = boost::tuples::get<1>(waveInfo);
		if(not rpwa::py::convertBPObjectToVector<int>(pyWaveRefls, waveRefls)) {
			PyErr_SetString(PyExc_TypeError, "Got invalid input for waveRefls when executing rpwa::fitResult::fill()");
			bp::throw_error_already_set();
		}
		const bp::list& pyListWaveProdAmpIndices = boost::tuples::get<2>(btWaveInfo);
		std::vector<std::vector<unsigned int> >& waveProdAmpIndices = boost::tuples::get<2>(waveInfo);
		waveProdAmpIndices.resize(bp::len(pyListWaveProdAmpIndices));
		for(int i = 0; i < bp::len(pyListWaveProdAmpIndices); ++i) {
			if(not rpwa::py::convertBPObjectToVector<unsigned int>(pyListWaveProdAmpIndices[i], waveProdAmpIndices[i]))
			{
				std::stringstream strStr;
				strStr<<"Could not convert element "<<i<<" when executing rpwa::fitResult::fill()";
				PyErr_SetString(PyExc_TypeError, strStr.str().c_str());
				bp::throw_error_already_set();
			}
		}
		boost::tuples::tuple<bp::list, bp::list, bp::list, bp::list> btProdAmpInfo;
		if(not rpwa::py::convertBPTupleToTuple<bp::list, bp::list, bp::list, bp::list>(pyProdAmpInfo, btProdAmpInfo)) {
			PyErr_SetString(PyExc_TypeError, "Got invalid input for prodAmpInfo when executing rpwa::fitResult::fill()");
			bp::throw_error_already_set();
		}
		rpwa::fitResult::prodAmpInfoType prodAmpInfo;
		const bp::list& pyProdAmpWaveIndices = boost::tuples::get<0>(btProdAmpInfo);
		std::vector<unsigned int>& prodAmpWaveIndices = boost::tuples::get<0>(prodAmpInfo);
		if(not rpwa::py::convertBPObjectToVector<unsigned int>(pyProdAmpWaveIndices, prodAmpWaveIndices)) {
			PyErr_SetString(PyExc_TypeError, "Got invalid input for prodAmpWaveIndices when executing rpwa::fitResult::fill()");
			bp::throw_error_already_set();
		}
		const bp::list& pyProdAmpRanks = boost::tuples::get<1>(btProdAmpInfo);
		std::vector<unsigned int>& prodAmpRanks = boost::tuples::get<1>(prodAmpInfo);
		if(not rpwa::py::convertBPObjectToVector<unsigned int>(pyProdAmpRanks, prodAmpRanks)) {
			PyErr_SetString(PyExc_TypeError, "Got invalid input for prodAmpRanks when executing rpwa::fitResult::fill()");
			bp::throw_error_already_set();
		}
		const bp::list& pyProdAmps = boost::tuples::get<2>(btProdAmpInfo);
		std::vector<std::complex<double> >& prodAmps = boost::tuples::get<2>(prodAmpInfo);
		if(not rpwa::py::convertBPObjectToVector<std::complex<double> >(pyProdAmps, prodAmps)) {
			PyErr_SetString(PyExc_TypeError, "Got invalid input for prodAmps when executing rpwa::fitResult::fill()");
			bp::throw_error_already_set();
		}
		const bp::list& pyListFitParCovMatrixIndices = boost::tuples::get<3>(btProdAmpInfo);
		std::vector<std::pair<int, int> >& fitParCovMatrixIndices = boost::tuples::get<3>(prodAmpInfo);
		fitParCovMatrixIndices.resize(bp::len(pyListFitParCovMatrixIndices));
		for(int i = 0; i < bp::len(pyListFitParCovMatrixIndices); ++i) {
			if(not rpwa::py::convertBPObjectToPair<int, int>(pyListFitParCovMatrixIndices[i], fitParCovMatrixIndices[i]))
			{
				std::stringstream strStr;
				strStr<<"Could not convert element "<<i<<" when executing rpwa::fitResult::fill()";
				PyErr_SetString(PyExc_TypeError, strStr.str().c_str());
				bp::throw_error_already_set();
			}
		}
		TMatrixT<double>* fitParCovMatrix = rpwa::py::convertFromPy<TMatrixT<double>* >(pyFitParCovMatrix);
		if(not fitParCovMatrix) {
			PyErr_SetString(PyExc_TypeError, "Got invalid input for fitParCovMatrix when executing rpwa::fitResult::fill()");
			bp::throw_error_already_set();
		}
		std::vector<double> phaseSpaceIntegral;
		if(not rpwa::py::convertBPObjectToVector<double>(pyPhaseSpaceIntegral, phaseSpaceIntegral)) {
			PyErr_SetString(PyExc_TypeError, "Got invalid input for phaseSpaceIntegral when executing rpwa::fitResult::fill()");
			bp::throw_error_already_set();
		}
		self.fill(nmbEvents, normNmbEvents, massBinCenter, logLikelihood, rank, waveInfo, prodAmpInfo, *fitParCovMatrix,
		          normIntegral, acceptedNormIntegral, phaseSpaceIntegral, converged, hasHessian);
	}

	void fitResult_fill_2(rpwa::fitResult& self, const rpwa::fitResult& result) {
		self.fill(result);
	}

	bp::list fitResult_evidenceComponents(const rpwa::fitResult& self) {
		return bp::list(self.evidenceComponents());
	}

	PyObject* fitResult_prodAmpCov_1(const rpwa::fitResult& self, const unsigned int prodAmpIndex)
	{
		return rpwa::py::convertToPy<TMatrixT<double> >(self.prodAmpCov(prodAmpIndex));
	}

	PyObject* fitResult_prodAmpCov_2(const rpwa::fitResult& self, const std::vector<unsigned int>& prodAmpIndices)
	{
		return rpwa::py::convertToPy<TMatrixT<double> >(self.prodAmpCov(prodAmpIndices));
	}

	double fitResult_phaseSpaceIntegral_1(const rpwa::fitResult& self, const unsigned int waveIndex)
	{
		return self.phaseSpaceIntegral(waveIndex);
	}

	double fitResult_phaseSpaceIntegral_2(const rpwa::fitResult& self, const std::string& waveName)
	{
		return self.phaseSpaceIntegral(waveName);
	}

	std::complex<double> fitResult_spinDensityMatrixElem_1(const rpwa::fitResult& self,
	                                                       const unsigned int waveIndexA,
	                                                       const unsigned int waveIndexB)
	{
		return self.spinDensityMatrixElem(waveIndexA, waveIndexB);
	}

	PyObject* fitResult_spinDensityMatrixElemCov_1(const rpwa::fitResult& self,
	                                               const unsigned int waveIndexA,
	                                               const unsigned int waveIndexB)
	{
		return rpwa::py::convertToPy<TMatrixT<double> >(self.spinDensityMatrixElemCov(waveIndexA, waveIndexB));
	}

	double fitResult_phase_1(const rpwa::fitResult& self,
	                         const unsigned int waveIndexA,
	                         const unsigned int waveIndexB)
	{
		return self.phase(waveIndexA, waveIndexB);
	}

	double fitResult_phaseErr_1(const rpwa::fitResult& self,
	                            const unsigned int waveIndexA,
	                            const unsigned int waveIndexB)
	{
		return self.phaseErr(waveIndexA, waveIndexB);
	}

	double fitResult_coherence_1(const rpwa::fitResult& self,
	                             const unsigned int waveIndexA,
	                             const unsigned int waveIndexB)
	{
		return self.coherence(waveIndexA, waveIndexB);
	}

	double fitResult_coherenceErr_1(const rpwa::fitResult& self,
	                                const unsigned int waveIndexA,
	                                const unsigned int waveIndexB)
	{
		return self.coherenceErr(waveIndexA, waveIndexB);
	}

	double fitResult_overlap_1(const rpwa::fitResult& self,
	                           const unsigned int waveIndexA,
	                           const unsigned int waveIndexB)
	{
		return self.overlap(waveIndexA, waveIndexB);
	}

	double fitResult_overlapErr_1(const rpwa::fitResult& self,
	                              const unsigned int waveIndexA,
	                              const unsigned int waveIndexB)
	{
		return self.overlapErr(waveIndexA, waveIndexB);
	}

	std::complex<double> fitResult_spinDensityMatrixElem_2(const rpwa::fitResult& self,
	                                                       const std::string& waveNameA,
	                                                       const std::string& waveNameB)
	{
		return self.spinDensityMatrixElem(waveNameA, waveNameB);
	}

	PyObject* fitResult_spinDensityMatrixElemCov_2(const rpwa::fitResult& self,
	                                               const std::string& waveNameA,
	                                               const std::string& waveNameB)
	{
		return rpwa::py::convertToPy<TMatrixT<double> >(self.spinDensityMatrixElemCov(waveNameA, waveNameB));
	}

	double fitResult_phase_2(const rpwa::fitResult& self,
	                         const std::string waveNameA,
	                         const std::string waveNameB)
	{
		return self.phase(waveNameA, waveNameB);
	}

	double fitResult_phaseErr_2(const rpwa::fitResult& self,
	                            const std::string waveNameA,
	                            const std::string waveNameB)
	{
		return self.phaseErr(waveNameA, waveNameB);
	}

	double fitResult_coherence_2(const rpwa::fitResult& self,
	                             const std::string waveNameA,
	                             const std::string waveNameB)
	{
		return self.coherence(waveNameA, waveNameB);
	}

	double fitResult_coherenceErr_2(const rpwa::fitResult& self,
	                                const std::string waveNameA,
	                                const std::string waveNameB)
	{
		return self.coherenceErr(waveNameA, waveNameB);
	}

	double fitResult_overlap_2(const rpwa::fitResult& self,
	                           const std::string waveNameA,
	                           const std::string waveNameB)
	{
		return self.overlap(waveNameA, waveNameB);
	}

	double fitResult_overlapErr_2(const rpwa::fitResult& self,
	                              const std::string waveNameA,
	                              const std::string waveNameB)
	{
		return self.overlapErr(waveNameA, waveNameB);
	}

	double fitResult_intensity_1(const rpwa::fitResult& self, const unsigned int waveIndex) {
		return self.intensity(waveIndex);
	}

	double fitResult_intensityErr_1(const rpwa::fitResult& self, const unsigned int waveIndex) {
		return self.intensityErr(waveIndex);
	}

	double fitResult_intensity_2(const rpwa::fitResult& self, const char* waveNamePattern) {
		return self.intensity(waveNamePattern);
	}

	double fitResult_intensityErr_2(const rpwa::fitResult& self, const char* waveNamePattern) {
		return self.intensityErr(waveNamePattern);
	}

	double fitResult_intensity_3(const rpwa::fitResult& self) {
		return self.intensity();
	}

	double fitResult_intensityErr_3(const rpwa::fitResult& self) {
		return self.intensityErr();
	}

	bp::list fitResult_prodAmps(const rpwa::fitResult& self)
	{
		bp::list retval;
		const std::vector<TComplex>& prodAmps = self.prodAmps();
		for(unsigned int i = 0; i < prodAmps.size(); ++i) {
			retval.append(std::complex<double>(prodAmps[i].Re(), prodAmps[i].Im()));
		}
		return retval;
	}

	bp::list fitResult_prodAmpRanks(const rpwa::fitResult& self)
	{
		return bp::list(self.prodAmpRanks());
	}

	bp::list fitResult_prodAmpWaveIndices(const rpwa::fitResult& self)
	{
		return bp::list(self.prodAmpWaveIndices());
	}

	bp::list fitResult_waveNames(const rpwa::fitResult& self)
	{
		return bp::list(self.waveNames());
	}

	bp::list fitResult_waveRefls(const rpwa::fitResult& self)
	{
		return bp::list(self.waveRefls());
	}

	bp::list fitResult_waveProdAmpIndices(const rpwa::fitResult& self)
	{
		bp::list retval;
		const std::vector<std::vector<unsigned int> >& waveProdAmpIndices = self.waveProdAmpIndices();
		for(unsigned int i = 0; i < waveProdAmpIndices.size(); ++i) {
			retval.append(bp::list(waveProdAmpIndices[i]));
		}
		return retval;
	}

	PyObject* fitResult_fitParCovMatrix(const rpwa::fitResult self)
	{
		return rpwa::py::convertToPy<TMatrixT<double> >(self.fitParCovMatrix());
	}

	bp::list fitResult_fitParCovIndices(const rpwa::fitResult self)
	{
		bp::list retval;
		const std::vector<std::pair<int, int> >& fitParCovIndices = self.fitParCovIndices();
		for(unsigned int i = 0; i < fitParCovIndices.size(); ++i) {
			const std::pair<int, int>& item = fitParCovIndices[i];
			retval.append(bp::make_tuple(item.first, item.second));
		}
		return retval;
	}

	bp::list fitResult_phaseSpaceIntegralVector(const rpwa::fitResult& self)
	{
		return bp::list(self.phaseSpaceIntegralVector());
	}

	bp::tuple fitResult_waveInfo(const rpwa::fitResult& self)
	{
		return bp::make_tuple(fitResult_waveNames(self),
		                      fitResult_waveRefls(self),
		                      fitResult_waveProdAmpIndices(self));
	}

	bp::tuple fitResult_prodAmpInfo(const rpwa::fitResult& self)
	{
		return bp::make_tuple(fitResult_prodAmpWaveIndices(self),
		                      fitResult_prodAmpRanks(self),
		                      fitResult_prodAmps(self),
		                      fitResult_fitParCovIndices(self));
	}

	std::string fitResult_printProdAmps(const rpwa::fitResult self)
	{
		std::stringstream sstr;
		self.printProdAmps(sstr);
		return sstr.str();
	}

	std::string fitResult_printWaves(const rpwa::fitResult self)
	{
		std::stringstream sstr;
		self.printWaves(sstr);
		return sstr.str();
	}

}

void rpwa::py::exportFitResult() {

	bp::def("escapeRegExpSpecialChar", &rpwa::escapeRegExpSpecialChar);
	bp::def("unescapeRegExpSpecialChar", &rpwa::unescapeRegExpSpecialChar);

	bp::class_<rpwa::fitResult>("fitResult")
		.def(bp::init<const rpwa::fitResult&>())
		.def(bp::self_ns::str(bp::self))
		.def("reset", &rpwa::fitResult::reset)
		.def("fill", &fitResult_fill_1)
		.def("fill", &fitResult_fill_2)
		.def("nmbEvents", &rpwa::fitResult::nmbEvents)
		.def("normNmbEvents", &rpwa::fitResult::normNmbEvents)
		.def("massBinCenter", &rpwa::fitResult::massBinCenter)
		.def("logLikelihood", &rpwa::fitResult::logLikelihood)
		.def("evidence", &rpwa::fitResult::evidence)
		.def("evidenceComponents", &fitResult_evidenceComponents)
		.def("rank", &rpwa::fitResult::rank)
		.def("covMatrixValid", &rpwa::fitResult::covMatrixValid)
		.def("converged", &rpwa::fitResult::converged)
		.def("hasHessian", &rpwa::fitResult::hasHessian)
		.def("nmbWaves", &rpwa::fitResult::nmbWaves)
		.def("nmbProdAmps", &rpwa::fitResult::nmbProdAmps)
		.def(
			"waveName"
			, &rpwa::fitResult::waveName
			, bp::return_value_policy<bp::copy_const_reference>()
		)
		.def("waveNameEsc", &rpwa::fitResult::waveNameEsc)
		.def("waveRefl", &rpwa::fitResult::waveRefl)
		.def("waveIndex", &rpwa::fitResult::waveIndex)
		.def(
			"waveNameForProdAmp"
			, &rpwa::fitResult::waveNameForProdAmp
			, bp::return_value_policy<bp::copy_const_reference>()
		)
		.def("waveIndexForProdAmp", &rpwa::fitResult::waveIndexForProdAmp)
		.def("rankOfProdAmp", &rpwa::fitResult::rankOfProdAmp)
		.def("reflOfProdAmp", &rpwa::fitResult::reflOfProdAmp)
		.def("fitParameter", &rpwa::fitResult::fitParameter)
		.def("prodAmp", &fitResult::prodAmp)
		.def("prodAmpCov", &fitResult_prodAmpCov_2)
		.def("prodAmpCov", &fitResult_prodAmpCov_1)
		.def("normIntegral", &rpwa::fitResult::normIntegral)
		.def("acceptedNormIntegral", &rpwa::fitResult::acceptedNormIntegral)
		.def("phaseSpaceIntegral", &fitResult_phaseSpaceIntegral_1)
		.def("phaseSpaceIntegral", &fitResult_phaseSpaceIntegral_2)

		.def("spinDensityMatrixElem", &fitResult_spinDensityMatrixElem_1)
		.def("spinDensityMatrixElemCov", &fitResult_spinDensityMatrixElemCov_1)
		.def("phase", &fitResult_phase_1)
		.def("phaseErr", &fitResult_phaseErr_1)
		.def("coherence", &fitResult_coherence_1)
		.def("coherenceErr", &fitResult_coherenceErr_1)
		.def("overlap", &fitResult_overlap_1)
		.def("overlapErr", &fitResult_overlapErr_1)

		.def("spinDensityMatrixElem", &fitResult_spinDensityMatrixElem_2)
		.def("spinDensityMatrixElemCov", &fitResult_spinDensityMatrixElemCov_2)
		.def("phase", &fitResult_phase_2)
		.def("phaseErr", &fitResult_phaseErr_2)
		.def("coherence", &fitResult_coherence_2)
		.def("coherenceErr", &fitResult_coherenceErr_2)
		.def("overlap", &fitResult_overlap_2)
		.def("overlapErr", &fitResult_overlapErr_2)

		.def("intensity", &fitResult_intensity_1)
		.def("intensityErr", &fitResult_intensityErr_1)
		.def("intensity", &fitResult_intensity_2)
		.def("intensityErr", &fitResult_intensityErr_2)
		.def("intensity", &fitResult_intensity_3)
		.def("intensityErr", &fitResult_intensityErr_3)

		.def("waveInfo", &fitResult_waveInfo)
		.def("prodAmpInfo", &fitResult_prodAmpInfo)
		.def("prodAmps", &fitResult_prodAmps)
		.def("prodAmpRanks", &fitResult_prodAmpRanks)
		.def("prodAmpWaveIndices", &fitResult_prodAmpWaveIndices)
		.def("waveNames", &fitResult_waveNames)
		.def("waveRefls", &fitResult_waveRefls)
		.def("waveProdAmpIndices", &fitResult_waveProdAmpIndices)
		.def("fitParCovMatrix", &fitResult_fitParCovMatrix)
		.def("fitParCovIndices", &fitResult_fitParCovIndices)
		.def(
			"normIntegralMatrix"
			, &rpwa::fitResult::normIntegralMatrix
			, bp::return_value_policy<bp::copy_const_reference>()
		)
		.def(
			"acceptedNormIntegralMatrix"
			, &rpwa::fitResult::acceptedNormIntegralMatrix
			, bp::return_value_policy<bp::copy_const_reference>()
		)
		.def("phaseSpaceIntegralVector", &fitResult_phaseSpaceIntegralVector)
		.def("printProdAmps", &fitResult_printProdAmps)
		.def("printWaves", &fitResult_printWaves)

		.def("setBranchAddress", &rpwa::py::setBranchAddress<rpwa::fitResult*>)
		.def(
			"branch"
			, &rpwa::py::branch<rpwa::fitResult*>
			, (bp::arg("fitResult"),
			   bp::arg("tree"),
			   bp::arg("name"),
			   bp::arg("bufsize")=32000,
			   bp::arg("splitlevel")=99)
		);

	bp::register_ptr_to_python<rpwa::fitResultPtr>();

}
