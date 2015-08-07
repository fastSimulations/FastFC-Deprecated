// Compute the voltage at the boundary of the cathode GDL/Channel (allows for dynamics calc rather than user input)
{
	scalar cathodeVoltage = 0;
	scalar anodeVoltage = 0;
	scalar cellVoltage = 0;

	scalar cellCurrentDensity = 0;

	forAll(CCLRegions, zoneID)
	{
		#include "CCLSetFields.H"
		
		label patchID = mesh.boundaryMesh().findPatchID("ccmCCL");

		const fvPatchField<scalar>& phiEPatchField = potElectron.boundaryField()[patchID];

		cathodeVoltage = gAverage(phiEPatchField);
	}

	forAll(ACLRegions, zoneID)
	{
		#include "ACLSetFields.H"
		
		label patchID = mesh.boundaryMesh().findPatchID("ccmACL");

		const fvPatchField<scalar>& phiEPatchField = potElectron.boundaryField()[patchID];

		anodeVoltage = gAverage(phiEPatchField);
	}

	cellVoltage = cathodeVoltage - anodeVoltage;
	cellCurrentDensity = (anodeCurrentDensity - cathodeCurrentDensity)/2.;


	// Output the polarization curve
	{
		std::ofstream polDat("polarizationData.txt", std::ios_base::app);
		polDat << nl << runTime.timeName() << ", " << cellVoltage << ", " << cellCurrentDensity;
	}
	
	// Output the current density
	//{
	//	fileName outputFile("currentDensity");
	//	OFstream os(runTime.path()/runTime.timeName()/outputFile, ios::app);
	//	os << "2" << nl << "(" << nl;
	//	os << runTime.timeName() << " " << cellCurrentDensity << nl << ");" << endl;
	//}
	// Output the cell voltage
	//{
	//	fileName outputFile("cellVoltage");
	//	OFstream os(runTime.path()/runTime.timeName()/outputFile);
	//	os << "2" << nl << "(" << nl;
	//	os << runTime.timeName() << cellVoltage << nl << ");" << endl;
	//}

}
