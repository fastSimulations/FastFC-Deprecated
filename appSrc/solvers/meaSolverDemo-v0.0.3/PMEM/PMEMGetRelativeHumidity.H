scalar pmemRH = 0.;

{
	label patchIDPMEMACL =mesh.boundaryMesh().findPatchID("PMEM_to_ACL");
	const fvPatch& primaryPatchPMEMACL = 
		refCast<const fvMesh>
   		(
       		mesh
		).boundary()[patchIDPMEMACL];

	const directMappedPatchBase& mppPMEMACL = 
		refCast<const directMappedPatchBase>
   		(
       		primaryPatchPMEMACL.patch()
   		);

	const polyMesh& neighbourMeshACLPMEM = mppPMEMACL.sampleMesh();

 		const fvPatch& neighbourPatchACLPMEM = 
		refCast<const fvMesh>
   		(
       		neighbourMeshACLPMEM
		).boundary()[mppPMEMACL.samplePolyPatch().index()];

	const fvPatchScalarField& anodeLocalRH = 
		refCast<const fvPatchScalarField>
		(
			neighbourPatchACLPMEM.lookupPatchField<volScalarField, scalar>
			(
				"RH"
			)
		);

	label patchIDPMEMCCL =mesh.boundaryMesh().findPatchID("PMEM_to_CCL");
	const fvPatch& primaryPatchPMEMCCL = 
		refCast<const fvMesh>
   		(
       		mesh
		).boundary()[patchIDPMEMCCL];

	const directMappedPatchBase& mppPMEMCCL = 
		refCast<const directMappedPatchBase>
   		(
       		primaryPatchPMEMCCL.patch()
   		);

	const polyMesh& neighbourMeshCCLPMEM = mppPMEMCCL.sampleMesh();

 		const fvPatch& neighbourPatchCCLPMEM = 
		refCast<const fvMesh>
   		(
       		neighbourMeshCCLPMEM
		).boundary()[mppPMEMCCL.samplePolyPatch().index()];

	const fvPatchScalarField& cathodeLocalRH = 
		refCast<const fvPatchScalarField>
		(
			neighbourPatchCCLPMEM.lookupPatchField<volScalarField, scalar>
			(
				"RH"
			)
		);

	scalar anodeRH_ = gAverage(anodeLocalRH);
	scalar cathodeRH_ = gAverage(cathodeLocalRH);

	pmemRH = (anodeRH_ + cathodeRH_)/2.;
}
