const wordList APTLRegionNames(rp["APTLRegionNames"]);

PtrList<fvMesh> APTLRegions(APTLRegionNames.size());

forAll(APTLRegionNames, zoneID)
{
	Info<< nl << APTLRegionNames[zoneID] << endl;

    APTLRegions.set
    (
        zoneID,
        new fvMesh
        (
            IOobject
            (
                APTLRegionNames[zoneID],
                runTime.timeName(),
                runTime,
                IOobject::MUST_READ
            )
        )
    );
}
