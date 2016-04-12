const wordList PMEMRegionNames(rp["PMEMRegionNames"]);

PtrList<fvMesh> PMEMRegions(PMEMRegionNames.size());

forAll(PMEMRegionNames, zoneID)
{
    Info<< nl << PMEMRegionNames[zoneID] << endl;

    PMEMRegions.set
    (
        zoneID,
        new fvMesh
        (
            IOobject
            (
                PMEMRegionNames[zoneID],
                runTime.timeName(),
                runTime,
                IOobject::MUST_READ
            )
        )
    );
}
