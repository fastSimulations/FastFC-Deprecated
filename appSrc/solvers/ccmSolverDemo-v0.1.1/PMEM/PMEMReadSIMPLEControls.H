const dictionary& simple = PMEMRegions[zoneID].solutionDict().subDict("SIMPLE");

const int nNonOrthCorr =
    simple.lookupOrDefault<int>("nNonOrthogonalCorrectors", 0);
