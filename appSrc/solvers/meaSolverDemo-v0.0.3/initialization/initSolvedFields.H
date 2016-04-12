/*---------------------------------------------------------------------------*\
    F ast           	    | FAST-FC: 
	is the		            | The Open Source Analysis and Simulation Toolbox 
    A nalysis and           | for Fuel Cells
	S imulation		        |
	Toolbox for		        | Copyright 2015, David B. Harvey
	F uel                   |         
	C ells                  |                     
-------------------------------------------------------------------------------
License
	FAST-FC and this file are a derivative work of OpenFOAM.

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

Application
   meaSolverDemo 

File
	initSolvedFields.H

Description
	User defined reading of initial guess

Notes
	1/ This file is currently used to list out the initial guess from change dict
	2/ Purpose of the file is to later allow a reading of a predetermined case file 
	time directory for transient

Developers
	David B. Harvey
\*---------------------------------------------------------------------------*/ 
// Start loop for initilization
{
	#include "readTimeControls.H"
//	runTime++;
	// Anode Initizalization Variables
	scalar potElectronAnodeInit_ = 0;
	scalar potProtonAnodeInit_ = 0;
	scalar wH2AnodeInit_ = 0;
	scalar wH2OVapAnodeInit_ = 0;
	scalar TAnodeInit_ = 0;
	scalar pMixAnodeInit_ = 0;
	scalar sLiqAnodeInit_ = 0;

	//Cathode Initilization Variables
	scalar potElectronCathodeInit_ = 0;
	scalar potProtonCathodeInit_ = 0;
	scalar wO2CathodeInit_ = 0;
	scalar wH2OVapCathodeInit_ = 0;
	scalar TCathodeInit_ = 0;
	scalar pMixCathodeInit_ = 0;
	scalar sLiqCathodeInit_ = 0;

	// Anode Initilization Loop
	forAll(APTLRegions, zoneID)
	{
		#include "APTLSetFields.H"	
		label patchID = mesh.boundaryMesh().findPatchID("anodeChannel");
		// Electrons
		const fvPatchField<scalar>& potElectronPatchField = potElectron.boundaryField()[patchID];
		potElectronAnodeInit_ = gAverage(potElectronPatchField);
		// Hydrogen
		const fvPatchField<scalar>& wH2PatchField = wH2.boundaryField()[patchID];
		wH2AnodeInit_ = gAverage(wH2PatchField);
		// Water Vapour
		const fvPatchField<scalar>& wH2OVapPatchField = wH2OVap.boundaryField()[patchID];
		wH2OVapAnodeInit_ = gAverage(wH2OVapPatchField);
		// Temperature
		const fvPatchField<scalar>& TPatchField = T.boundaryField()[patchID];
		TAnodeInit_ = gAverage(TPatchField);
		// Mixture Pressure
		const fvPatchField<scalar>& pMixPatchField = pMix.boundaryField()[patchID];
		pMixAnodeInit_ = gAverage(pMixPatchField);
		// Liquid Water Saturation
		const fvPatchField<scalar>& sLiqPatchField = sLiq.boundaryField()[patchID];
		sLiqAnodeInit_ = gAverage(sLiqPatchField);
	}

	// Cathode Initilization Loop
	forAll(CPTLRegions, zoneID)
	{
		#include "CPTLSetFields.H"	
		label patchID = mesh.boundaryMesh().findPatchID("cathodeChannel");
		// Electrons
		const fvPatchField<scalar>& potElectronPatchField = potElectron.boundaryField()[patchID];
		potElectronCathodeInit_ = gAverage(potElectronPatchField);
		// Hydrogen
		const fvPatchField<scalar>& wO2PatchField = wO2.boundaryField()[patchID];
		wO2CathodeInit_ = gAverage(wO2PatchField);
		// Water Vapour
		const fvPatchField<scalar>& wH2OVapPatchField = wH2OVap.boundaryField()[patchID];
		wH2OVapCathodeInit_ = gAverage(wH2OVapPatchField);
		// Temperature
		const fvPatchField<scalar>& TPatchField = T.boundaryField()[patchID];
		TCathodeInit_ = gAverage(TPatchField);
		// Mixture Pressure
		const fvPatchField<scalar>& pMixPatchField = pMix.boundaryField()[patchID];
		pMixCathodeInit_ = gAverage(pMixPatchField);
		// Liquid Water Saturation
		const fvPatchField<scalar>& sLiqPatchField = sLiq.boundaryField()[patchID];
		sLiqCathodeInit_ = gAverage(sLiqPatchField);
	}

	// Output the Averaged Settings
	{
		Info<< nl
			<<	tab << "Time = " << runTime.timeName() << nl
			<<	tab << "Anode Operating Conditions" << nl
			<<  tab << tab << "wH2 = " << wH2AnodeInit_ << nl
			<<  tab << tab << "wH2OVap = " << wH2OVapAnodeInit_ << nl
			<<  tab << tab << "Temperature = " << TAnodeInit_ << nl
			<<  tab << tab << "Mixture Pressure = " << pMixAnodeInit_ << nl
			<<  tab << tab << "Saturation = " << sLiqAnodeInit_ << nl
			<<	tab << "Cathode Operating Conditions" << nl
			<<  tab << tab << "wO2 = " << wO2CathodeInit_ << nl
			<<  tab << tab << "wH2OVap = " << wH2OVapCathodeInit_ << nl
			<<  tab << tab << "Temperature = " << TCathodeInit_ << nl
			<<  tab << tab << "Mixture Pressure = " << pMixCathodeInit_ << nl
			<<  tab << tab << "Saturation = " << sLiqCathodeInit_ << nl
			<< endl;	
	}

	// Set the internal and boundary fields in all components

	// APTL
	forAll(APTLRegions, zoneID)
	{
		#include "APTLSetFields.H"
		label patchIdACL = mesh.boundaryMesh().findPatchID("APTL_to_ACL");
		const polyPatch& patchACL = mesh.boundaryMesh()[patchIdACL];

		// Set Cell Centres
		{
			// Electrons
			potElectron.internalField() = potElectronAnodeInit_;
			// Hydrogen
			wH2.internalField() = wH2AnodeInit_;
			// Water Vapour
			wH2OVap.internalField() = wH2OVapAnodeInit_;
			// Temperature
			T.internalField() = TAnodeInit_;
			// Mixture Pressure
			pMix.internalField() = pMixAnodeInit_;
			// Liquid Water Saturation
			sLiq.internalField() = sLiqAnodeInit_;
		}

		// Set Boundary Patch ACL
		forAll(patchACL, faceID)
		{
			// Electrons
			potElectron.boundaryField()[patchIdACL][faceID] = potElectronAnodeInit_;
			// Hydrogen
			wH2.boundaryField()[patchIdACL][faceID] = wH2AnodeInit_;
			// Water Vapour
			wH2OVap.boundaryField()[patchIdACL][faceID] = wH2OVapAnodeInit_;
			// Temperature
			T.boundaryField()[patchIdACL][faceID] = TAnodeInit_;
			// Mixture Pressure
			pMix.boundaryField()[patchIdACL][faceID] = pMixAnodeInit_;
			// Liquid Water Saturation
			sLiq.boundaryField()[patchIdACL][faceID] = sLiqAnodeInit_;
		}

	}

	// ACL
	forAll(ACLRegions, zoneID)
	{
		#include "ACLSetFields.H"
		label patchIdAPTL = mesh.boundaryMesh().findPatchID("ACL_to_APTL");
		const polyPatch& patchAPTL = mesh.boundaryMesh()[patchIdAPTL];
		
		label patchIdPMEM = mesh.boundaryMesh().findPatchID("ACL_to_PMEM");
		const polyPatch& patchPMEM = mesh.boundaryMesh()[patchIdPMEM];

		// Set Cell Centres
		{
			// Electrons
			potElectron.internalField() = potElectronAnodeInit_;
			// Protons
			potProton.internalField() = 0.;
			// Hydrogen
			wH2.internalField() = wH2AnodeInit_;
			// Water Vapour
			wH2OVap.internalField() = wH2OVapAnodeInit_;
			// Temperature
			T.internalField() = TAnodeInit_;
			// Mixture Pressure
			pMix.internalField() = pMixAnodeInit_;
			// Liquid Water Saturation
			sLiq.internalField() = sLiqAnodeInit_;
		}

		// Set Boundary PatchAPTL
		forAll(patchAPTL, faceID)
		{
			// Electrons
			potElectron.boundaryField()[patchIdAPTL][faceID] = potElectronAnodeInit_;
			// Protons
			potProton.boundaryField()[patchIdAPTL][faceID] = 0.;
			// Hydrogen
			wH2.boundaryField()[patchIdAPTL][faceID] = wH2AnodeInit_;
			// Water Vapour
			wH2OVap.boundaryField()[patchIdAPTL][faceID] = wH2OVapAnodeInit_;
			// Temperature
			T.boundaryField()[patchIdAPTL][faceID] = TAnodeInit_;
			// Mixture Pressure
			pMix.boundaryField()[patchIdAPTL][faceID] = pMixAnodeInit_;
			// Liquid Water Saturation
			sLiq.boundaryField()[patchIdAPTL][faceID] = sLiqAnodeInit_;
		}
	
		// Set Boundary PatchPMEM
		forAll(patchPMEM, faceID)
		{
			// Electrons
			potElectron.boundaryField()[patchIdPMEM][faceID] = potElectronAnodeInit_;
			// Protons
			potProton.boundaryField()[patchIdPMEM][faceID] = 0.;
			// Hydrogen
			wH2.boundaryField()[patchIdPMEM][faceID] = wH2AnodeInit_;
			// Water Vapour
			wH2OVap.boundaryField()[patchIdPMEM][faceID] = wH2OVapAnodeInit_;
			// Temperature
			T.boundaryField()[patchIdPMEM][faceID] = TAnodeInit_;
			// Mixture Pressure
			pMix.boundaryField()[patchIdPMEM][faceID] = pMixAnodeInit_;
			// Liquid Water Saturation
			sLiq.boundaryField()[patchIdPMEM][faceID] = sLiqAnodeInit_;
		}
	}

	// PMEM
	forAll(PMEMRegions, zoneID)
	{
		#include "PMEMSetFields.H"
		label patchIdACL = mesh.boundaryMesh().findPatchID("PMEM_to_ACL");
		const polyPatch& patchACL = mesh.boundaryMesh()[patchIdACL];

		label patchIdCCL = mesh.boundaryMesh().findPatchID("PMEM_to_CCL");
		const polyPatch& patchCCL = mesh.boundaryMesh()[patchIdCCL];

		// Set Cell Centres
		{
			// Protons
			potProton.internalField() = 0.;
			// Temperature
			T.internalField() = (TAnodeInit_ + TCathodeInit_)/2.;
		}

		// Set Boundary Patch ACL
		forAll(patchACL, faceID)
		{
			// Protons
			potProton.boundaryField()[patchIdACL][faceID] = 0.;
			// Temperature
			T.boundaryField()[patchIdACL][faceID] = (TAnodeInit_ + TCathodeInit_)/2.;
		}
		
		// Set Boundary Patch CCL
		forAll(patchCCL, faceID)
		{
			// Protons
			potProton.boundaryField()[patchIdCCL][faceID] = 0.;
			// Temperature
			T.boundaryField()[patchIdCCL][faceID] = (TAnodeInit_ + TCathodeInit_)/2.;
		}
	}
	
	// CCL
	forAll(CCLRegions, zoneID)
	{
		#include "CCLSetFields.H"

		label patchIdPMEM = mesh.boundaryMesh().findPatchID("CCL_to_PMEM");
		const polyPatch& patchPMEM = mesh.boundaryMesh()[patchIdPMEM];

		label patchIdCPTL = mesh.boundaryMesh().findPatchID("CCL_to_CPTL");
		const polyPatch& patchCPTL = mesh.boundaryMesh()[patchIdCPTL];

		// Set Cell Centres
		{
			// Electrons
			potElectron.internalField() = potElectronCathodeInit_;
			// Protons
			potProton.internalField() = 0.;
			// Oxygen
			wO2.internalField() = wO2CathodeInit_;
			// Water Vapour
			wH2OVap.internalField() = wH2OVapCathodeInit_;
			// Temperature
			T.internalField() = TCathodeInit_;
			// Mixture Pressure
			pMix.internalField() = pMixCathodeInit_;
			// Liquid Water Saturation
			sLiq.internalField() = sLiqCathodeInit_;
		}

		// Set Boundary Patches
		forAll(patchPMEM, faceID)
		{
			// Electrons
			potElectron.boundaryField()[patchIdPMEM][faceID] = potElectronCathodeInit_;
			// Protons
			potProton.boundaryField()[patchIdPMEM][faceID] = 0.;
			// Oxygen
			wO2.boundaryField()[patchIdPMEM][faceID] = wO2CathodeInit_;
			// Water Vapour
			wH2OVap.boundaryField()[patchIdPMEM][faceID] = wH2OVapCathodeInit_;
			// Temperature
			T.boundaryField()[patchIdPMEM][faceID] = TCathodeInit_;
			// Mixture Pressure
			pMix.boundaryField()[patchIdPMEM][faceID] = pMixCathodeInit_;
			// Liquid Water Saturation
			sLiq.boundaryField()[patchIdPMEM][faceID] = sLiqCathodeInit_;
		}

		forAll(patchCPTL, faceID)
		{
			// Electrons
			potElectron.boundaryField()[patchIdCPTL][faceID] = potElectronCathodeInit_;
			// Protons
			potProton.boundaryField()[patchIdCPTL][faceID] = 0.;
			// Oxygen
			wO2.boundaryField()[patchIdCPTL][faceID] = wO2CathodeInit_;
			// Water Vapour
			wH2OVap.boundaryField()[patchIdCPTL][faceID] = wH2OVapCathodeInit_;
			// Temperature
			T.boundaryField()[patchIdCPTL][faceID] = TCathodeInit_;
			// Mixture Pressure
			pMix.boundaryField()[patchIdCPTL][faceID] = pMixCathodeInit_;
			// Liquid Water Saturation
			sLiq.boundaryField()[patchIdCPTL][faceID] = sLiqCathodeInit_;
		}
	}

	// CPTL
	forAll(CPTLRegions, zoneID)
	{
		#include "CPTLSetFields.H"
		label patchIdCCL = mesh.boundaryMesh().findPatchID("CPTL_to_CCL");
		const polyPatch& patchCCL = mesh.boundaryMesh()[patchIdCCL];

		// Set Cell Centres
		{
			// Electrons
			potElectron.internalField() = potElectronCathodeInit_;
			// Oxygen
			wO2.internalField() = wO2CathodeInit_;
			// Water Vapour
			wH2OVap.internalField() = wH2OVapCathodeInit_;
			// Temperature
			T.internalField() = TCathodeInit_;
			// Mixture Pressure
			pMix.internalField() = pMixCathodeInit_;
			// Liquid Water Saturation
			sLiq.internalField() = sLiqCathodeInit_;
		}

		// Set Boundary Patch CCL
		forAll(patchCCL, faceID)
		{
			// Electrons
			potElectron.boundaryField()[patchIdCCL][faceID] = potElectronCathodeInit_;
			// Oxygen
			wO2.boundaryField()[patchIdCCL][faceID] = wO2CathodeInit_;
			// Water Vapour
			wH2OVap.boundaryField()[patchIdCCL][faceID] = wH2OVapCathodeInit_;
			// Temperature
			T.boundaryField()[patchIdCCL][faceID] = TCathodeInit_;
			// Mixture Pressure
			pMix.boundaryField()[patchIdCCL][faceID] = pMixCathodeInit_;
			// Liquid Water Saturation
			sLiq.boundaryField()[patchIdCCL][faceID] = sLiqCathodeInit_;
		}
	}
	
	// Correct Solved Boundary Conditions
	#include "initCorrectSolvedBoundaryConditions.H"

	// Write the initialization fields
	#include "initWriteSolvedFields.H"
}
