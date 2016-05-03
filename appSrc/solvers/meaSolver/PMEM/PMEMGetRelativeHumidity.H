/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    F ast           	    | FAST-FC: 
	is the		            | The Open Source Analysis and Simulation Toolbox 
    A nalysis and           | for Fuel Cells
	S imulation		        |
	Toolbox for		        | Copyright 2015, David B. Harvey
	F uel                   |         
	C ells                  |                     
-------------------------------------------------------------------------------
License
	FAST-FC is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    FAST-FC is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with FAST-FC.  If not, see <http://www.gnu.org/licenses/>.

	This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Application
    - PMEMGetRelativeHumidity.H

Description
	- Looks up the face value of the RH at the CL/PMEM boundaries and uses 
	that to provide an ave. RH for the polymer membrane
:
\*---------------------------------------------------------------------------*/


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

	scalar anodeRH_ = gAverage(anodeLocalRH);  // Limits are imposed in general fields
	scalar cathodeRH_ = gAverage(cathodeLocalRH); // Limits are imposed in general fields

	pmemRH = (anodeRH_ + cathodeRH_)/2.;
}
