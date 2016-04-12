const wordList CPTLRegionNames(rp["CPTLRegionNames"]);

PtrList<fvMesh> CPTLRegions(CPTLRegionNames.size());

forAll(CPTLRegionNames, zoneID)
{
	Info<< nl << CPTLRegionNames[zoneID] << endl;

    CPTLRegions.set
    (
        zoneID,
        new fvMesh
        (
            IOobject
            (
                CPTLRegionNames[zoneID],
                runTime.timeName(),
                runTime,
                IOobject::MUST_READ
            )
        )
    );
}
