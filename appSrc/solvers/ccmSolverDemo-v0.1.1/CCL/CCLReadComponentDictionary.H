// Read in all related dictionaries
IOdictionary transportProperties
(
    IOobject
    (
		"transportProperties",
		runTime.constant(),
		CCLRegions[zoneID],
		IOobject::MUST_READ,
		IOobject::NO_WRITE
	)
);
