const wordList CCLRegionNames(rp["CCLRegionNames"]);

PtrList<fvMesh> CCLRegions(CCLRegionNames.size());

forAll(CCLRegionNames, zoneID)
{
    Info<< nl << CCLRegionNames[zoneID] << endl;

    CCLRegions.set
    (
        zoneID,
        new fvMesh
        (
            IOobject
            (
                CCLRegionNames[zoneID],
                runTime.timeName(),
                runTime,
                IOobject::MUST_READ
            )
        )
    );
}
