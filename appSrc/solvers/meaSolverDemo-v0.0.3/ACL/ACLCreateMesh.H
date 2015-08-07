const wordList ACLRegionNames(rp["ACLRegionNames"]);

PtrList<fvMesh> ACLRegions(ACLRegionNames.size());

forAll(ACLRegionNames, zoneID)
{
    Info<< nl << ACLRegionNames[zoneID] << endl;

    ACLRegions.set
    (
        zoneID,
        new fvMesh
        (
            IOobject
            (
                ACLRegionNames[zoneID],
                runTime.timeName(),
                runTime,
                IOobject::MUST_READ
            )
        )
    );
}
