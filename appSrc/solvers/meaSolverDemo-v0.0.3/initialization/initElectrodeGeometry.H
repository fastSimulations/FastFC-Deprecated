// Create dimensioned variables for inline determination of thickness and volume for each layer

dimensionedScalar geomActiveArea("geomActiveArea", dimensionSet( 0, 2, 0, 0, 0, 0, 0), 0.);
dimensionedScalar thicknessAPTL("thicknessAPTL", dimensionSet( 0, 1, 0, 0, 0, 0, 0), 0.);
dimensionedScalar thicknessACL("thicknessACL", dimensionSet( 0, 1, 0, 0, 0, 0, 0), 0.);
dimensionedScalar thicknessPMEM("thicknessPMEM", dimensionSet( 0, 1, 0, 0, 0, 0, 0), 0.);
dimensionedScalar thicknessCCL("thicknessCCL", dimensionSet( 0, 1, 0, 0, 0, 0, 0), 0.);
dimensionedScalar thicknessCPTL("thicknessCPTL", dimensionSet( 0, 1, 0, 0, 0, 0, 0), 0.);
dimensionedScalar volumeAPTL("volumeAPTL", dimensionSet(0, 3, 0, 0, 0, 0, 0), 0.);
dimensionedScalar volumeACL("volumeACL", dimensionSet(0, 3, 0, 0, 0, 0, 0), 0.);
dimensionedScalar volumePMEM("volumePMEM", dimensionSet(0, 3, 0, 0, 0, 0, 0), 0.);
dimensionedScalar volumeCCL("volumeCCL", dimensionSet(0, 3, 0, 0, 0, 0, 0), 0.);
dimensionedScalar volumeCPTL("volumeCPTL", dimensionSet(0, 3, 0, 0, 0, 0, 0), 0.);

// Sub routine to determine the inline thickness, volume, and geometric active area
{
	scalar activeArea_ = 0.;
	scalar thicknessAPTL_ = 0.;
	scalar volumeAPTL_ = 0.;
	scalar thicknessACL_ = 0;
	scalar  volumeACL_ = 0.;
	scalar thicknessPMEM_ = 0.;
	scalar volumePMEM_ = 0.;
	scalar thicknessCCL_ = 0.;
	scalar volumeCCL_ = 0.;
	scalar thicknessCPTL_ = 0.;
	scalar volumeCPTL_ = 0.;

	// Compute the geometric cross-sectional area (allows for dynamics calc rather than user input)
	forAll(APTLRegions, zoneID)
	{
		fvMesh& mesh = APTLRegions[zoneID];
		label patchID = mesh.boundaryMesh().findPatchID("anodeChannel");
		const polyPatch& cPatch = mesh.boundaryMesh()[patchID];
		const surfaceScalarField& magSf = mesh.magSf();
		scalarField volume = mesh.V();
		activeArea_ = gSum(magSf.boundaryField()[patchID]);
		volumeAPTL_ = gSum(volume);
		//	Info << "cPATCH = " << nl << cPatch << endl;
		//	Info << "patchID = " << patchID << endl;
	}
	forAll(ACLRegions, zoneID)
	{
		fvMesh& mesh = ACLRegions[zoneID];
		label patchID = mesh.boundaryMesh().findPatchID("ACL_to_PMEM");
		const polyPatch& cPatch = mesh.boundaryMesh()[patchID];
		const surfaceScalarField& magSf = mesh.magSf();
		scalarField volume = mesh.V();
		volumeACL_ = gSum(volume);
	}
	forAll(PMEMRegions, zoneID)
	{
		fvMesh& mesh = PMEMRegions[zoneID];
		label patchID = mesh.boundaryMesh().findPatchID("PMEM_to_CCL");
		const polyPatch& cPatch = mesh.boundaryMesh()[patchID];
		const surfaceScalarField& magSf = mesh.magSf();
		scalarField volume = mesh.V();
		volumePMEM_ = gSum(volume);
	}
	forAll(CCLRegions, zoneID)
	{
		fvMesh& mesh = CCLRegions[zoneID];
		label patchID = mesh.boundaryMesh().findPatchID("CCL_to_CPTL");
		const polyPatch& cPatch = mesh.boundaryMesh()[patchID];
		const surfaceScalarField& magSf = mesh.magSf();
		scalarField volume = mesh.V();
		volumeCCL_ = gSum(volume);
	}
	forAll(CPTLRegions, zoneID)
	{
		fvMesh& mesh = CPTLRegions[zoneID];
		label patchID = mesh.boundaryMesh().findPatchID("cathodeChannel");
		const polyPatch& cPatch = mesh.boundaryMesh()[patchID];
		const surfaceScalarField& magSf = mesh.magSf();
		scalarField volume = mesh.V();
		volumeCPTL_ = gSum(volume);
	}

	geomActiveArea.value() = activeArea_;
	thicknessAPTL.value() = (volumeAPTL_/activeArea_);
	thicknessACL.value() = (volumeACL_/activeArea_);
	thicknessPMEM.value() = (volumePMEM_/activeArea_);
	thicknessCCL.value() = (volumeCCL_/activeArea_);
	thicknessCPTL.value() = (volumeCPTL_/activeArea_);
	volumeAPTL.value() = volumeAPTL_;
	volumeACL.value() = volumeACL_;
	volumePMEM.value() = volumePMEM_;
	volumeCCL.value() = volumeCCL_;
	volumeCPTL.value() = volumeCPTL_;
}

Info << nl << endl;
Info << "     Geometric Active Area = " << geomActiveArea.value() << " [m^2]" << endl;
Info << "       Anode PTL thickness = " << (thicknessAPTL.value()*1e6) << " [micrometers]" << endl;
Info << "       Anode  CL thickness = " << (thicknessACL.value()*1e6) << " [micrometers]" << endl;
Info << "     Polymer MEM thickness = " << (thicknessPMEM.value()*1e6) << " [micrometers]"<< endl;
Info << "     Cathode  CL thickness = " << (thicknessCCL.value()*1e6) << " [micrometers]" << endl;
Info << "     Cathode PTL thickness = " << (thicknessCPTL.value()*1e6) << " [micrometers]" << endl;
