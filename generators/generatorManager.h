#ifndef GENERATORMANAGER_HH_
#define GENERATORMANAGER_HH_


#include<string>

#include "generatorParameters.hpp"


class TVector3;

namespace rpwa {

	class generator;
	class massAndTPrimePicker;
	class beamAndVertexGenerator;

	class generatorManager {

	  public:

		generatorManager();
		~generatorManager();

		unsigned int event();

		const rpwa::generator& getGenerator() const { return *_generator; }

		bool readReactionFile(const std::string& fileName);

		bool initializeGenerator();

		void overrideMassRange(double lowerLimit, double upperLimit);
		void overrideBeamFile(std::string beamFileName) { _beamFileName = beamFileName; }
		void readBeamfileSequentially(bool readBeamfileSequentially = true);
		void randomizeBeamfileStartingPosition();

		std::ostream& print(std::ostream& out);

		static bool debug() { return _debug; }
		static void setDebug(bool debug = true) { _debug = debug; }

	  private:

		rpwa::Beam _beam;
		rpwa::Target _target;
		rpwa::FinalState _finalState;

		rpwa::beamAndVertexGenerator* _beamAndVertexGenerator;

		rpwa::massAndTPrimePicker* _pickerFunction;

		std::string _beamFileName;

		bool _reactionFileRead;

		rpwa::generator* _generator;

		static bool _debug;

	};

}

#endif

