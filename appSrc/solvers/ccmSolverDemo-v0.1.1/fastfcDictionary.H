// Read in all related dictionaries

IOdictionary APOLLOControlDict
(
    IOobject
    (
		"fastfcControlDict",
		runTime.system(),
		runTime,
		IOobject::MUST_READ,
		IOobject::NO_WRITE
	)
);

const scalar& maxSweep = readScalar(APOLLOControlDict.lookup("maxSweep"));
const scalar& currentTolerance = readScalar(APOLLOControlDict.lookup("currentTolerance"));
const scalar& THV = readScalar(APOLLOControlDict.lookup("THV")) ;
const scalar& slopeStep = readScalar(APOLLOControlDict.lookup("slopeStep"));
