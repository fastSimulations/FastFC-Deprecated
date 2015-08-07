// Determine the geometric cross-sectional Area
scalar activeArea = 0.;
scalar electrodeVol = 0.;

// Compute the geometric cross-sectional area (allows for dynamics calc rather than user input)
forAll(ACLRegions, zoneID)
{
	#include "ACLSetFields.H"
	
	label patchID = mesh.boundaryMesh().findPatchID("ACL_to_PMEM");

	Info << "patchID = " << patchID << endl;

	const polyPatch& cPatch = mesh.boundaryMesh()[patchID];

	Info << "cPATCH = " << nl << cPatch << endl;

	const surfaceScalarField& magSf = mesh.magSf();

	scalarField volume = mesh.V();

	activeArea = gSum(magSf.boundaryField()[patchID]);
	electrodeVol = gSum(volume);
}
	
Info << "Active Area = " <<  activeArea << "  [m^2]" << endl;
Info << "Anode Volume = " << electrodeVol << "  [m^3]" << endl;
