// Read in all related dictionaries
IOdictionary transportProperties
(
    IOobject
    (
		"transportProperties",
		runTime.constant(),
		APTLRegions[zoneID],
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	)
);
