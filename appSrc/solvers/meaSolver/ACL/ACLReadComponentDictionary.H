// Read in all related dictionaries
IOdictionary transportProperties
(
    IOobject
    (
		"transportProperties",
		runTime.constant(),
		ACLRegions[zoneID],
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	)
);
