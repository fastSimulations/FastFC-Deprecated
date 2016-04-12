// Read in all related dictionaries
IOdictionary transportProperties
(
    IOobject
    (
		"transportProperties",
		runTime.constant(),
		CPTLRegions[zoneID],
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	)
);
