// Read in all related dictionaries

IOdictionary fastControlDict
(
    IOobject
    (
		"fastControlDict",
		runTime.system(),
		runTime,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	)
);

// General Settings
const scalar& maxSweep = readScalar(fastControlDict.lookup("maxSweep"));
const scalar& currentTolerance = readScalar(fastControlDict.lookup("currentTolerance"));
const scalar& slopeStep = readScalar(fastControlDict.lookup("slopeStep"));
const scalar& relaxElectrons = readScalar(fastControlDict.lookup("electronRelaxFactor"));
const scalar& relaxProtons = readScalar(fastControlDict.lookup("protonRelaxFactor"));
const scalar& relaxOxygen = readScalar(fastControlDict.lookup("oxygenRelaxFactor"));
const scalar& relaxHydrogen = readScalar(fastControlDict.lookup("hydrogenRelaxFactor"));
const scalar& relaxWaterVapour = readScalar(fastControlDict.lookup("waterVapourRelaxFactor"));
const scalar& relaxMass = readScalar(fastControlDict.lookup("massRelaxFactor"));
const scalar& relaxAnodeKinetics = readScalar(fastControlDict.lookup("anodeKineticsRelaxFactor"));
const scalar& relaxCathodeKinetics = readScalar(fastControlDict.lookup("cathodeKineticsRelaxFactor"));
// Equation Enabling
const scalar& solveOxygen = readScalar(fastControlDict.lookup("solveOxygenTransport"));
const scalar& solveHydrogen = readScalar(fastControlDict.lookup("solveHydrogenTransport"));
const scalar& solveWaterVapour = readScalar(fastControlDict.lookup("solveWaterVapourTransport"));
const scalar& solveMass = readScalar(fastControlDict.lookup("solveMassTransport"));
// Anode Settings
const scalar& anodeKineticsTHV = readScalar(fastControlDict.lookup("anodeKineticsTHV"));
const scalar& anodeInputECSA = readScalar(fastControlDict.lookup("anodeInputECSA"));
const scalar& anodeSinglePhase = readScalar(fastControlDict.lookup("anodeSinglePhase"));
// Cathode Settings
const scalar& cathodeKineticsAgglomerate = readScalar(fastControlDict.lookup("cathodeKineticsAgg"));
const scalar& cathodeInputECSA = readScalar(fastControlDict.lookup("cathodeInputECSA"));
const scalar& cathodeSinglePhase = readScalar(fastControlDict.lookup("cathodeSinglePhase"));
