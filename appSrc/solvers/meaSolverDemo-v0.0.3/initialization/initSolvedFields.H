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
	forAll(APTLRegions, zoneID)
	{
		#include "APTLSetFields.H"
		// Electrons
		potElectron.internalField() = potElectronAnodeInit_;
		potElectron.boundaryField() = potElectronAnodeInit_;
		// Hydrogen
		wH2.internalField() = wH2AnodeInit_;
		wH2.boundaryField() = wH2AnodeInit_;
		// Water Vapour
		wH2OVap.internalField() = wH2OVapAnodeInit_;
		wH2OVap.boundaryField() = wH2OVapAnodeInit_;
		// Temperature
		T.internalField() = TAnodeInit_;
		T.boundaryField() = TAnodeInit_;
		// Mixture Pressure
		pMix.internalField() = pMixAnodeInit_;
		pMix.boundaryField() = pMixAnodeInit_;
		// Liquid Water Saturation
		sLiq.internalField() = sLiqAnodeInit_;
		sLiq.boundaryField() = sLiqAnodeInit_;
	}

	forAll(ACLRegions, zoneID)
	{
		#include "ACLSetFields.H"
		// Electrons
		potElectron.internalField() = potElectronAnodeInit_;
		potElectron.boundaryField() = potElectronAnodeInit_;
		// Protons
		potProton.internalField() = 0.;
		potProton.boundaryField() = 0.;
		// Hydrogen
		wH2.internalField() = wH2AnodeInit_;
		wH2.boundaryField() = wH2AnodeInit_;
		// Water Vapour
		wH2OVap.internalField() = wH2OVapAnodeInit_;
		wH2OVap.boundaryField() = wH2OVapAnodeInit_;
		// Temperature
		T.internalField() = TAnodeInit_;
		T.boundaryField() = TAnodeInit_;
		// Mixture Pressure
		pMix.internalField() = pMixAnodeInit_;
		pMix.boundaryField() = pMixAnodeInit_;
		// Liquid Water Saturation
		sLiq.internalField() = sLiqAnodeInit_;
		sLiq.boundaryField() = sLiqAnodeInit_;
	}

	forAll(PMEMRegions, zoneID)
	{
		#include "PMEMSetFields.H"
		label patchIDAnode = mesh.boundaryMesh().findPatchID("PMEM_to_ACL");
		label patchIDCathode = mesh.boundaryMesh().findPatchID("PMEM_to_CCL");
		// Protons
		potProton.internalField() = 0.;
		potProton.boundaryField() = 0.;
		// Temperature
		T.internalField() = (TAnodeInit_ + TCathodeInit_)/2.;
		T.boundaryField() = (TAnodeInit_ + TCathodeInit_)/2.;
	}
	
	forAll(CCLRegions, zoneID)
	{
		#include "CCLSetFields.H"
		// Electrons
		potElectron.internalField() = potElectronCathodeInit_;
		potElectron.boundaryField() = potElectronCathodeInit_;
		// Protons
		potProton.internalField() = 0.;
		potProton.boundaryField() = 0.;
		// Oxygen
		wO2.internalField() = wO2CathodeInit_;
		wO2.boundaryField() = wO2CathodeInit_;
		// Water Vapour
		wH2OVap.internalField() = wH2OVapCathodeInit_;
		wH2OVap.boundaryField() = wH2OVapCathodeInit_;
		// Temperature
		T.internalField() = TCathodeInit_;
		T.boundaryField() = TCathodeInit_;
		// Mixture Pressure
		pMix.internalField() = pMixCathodeInit_;
		pMix.boundaryField() = pMixCathodeInit_;
		// Liquid Water Saturation
		sLiq.internalField() = sLiqCathodeInit_;
		sLiq.boundaryField() = sLiqCathodeInit_;
	}

	forAll(CPTLRegions, zoneID)
	{
		#include "CPTLSetFields.H"
		// Electrons
		potElectron.internalField() = potElectronCathodeInit_;
		potElectron.boundaryField() = potElectronCathodeInit_;
		// Oxygen
		wO2.internalField() = wO2CathodeInit_;
		wO2.boundaryField() = wO2CathodeInit_;
		// Water Vapour
		wH2OVap.internalField() = wH2OVapCathodeInit_;
		wH2OVap.boundaryField() = wH2OVapCathodeInit_;
		// Temperature
		T.internalField() = TCathodeInit_;
		T.boundaryField() = TCathodeInit_;
		// Mixture Pressure
		pMix.internalField() = pMixCathodeInit_;
		pMix.boundaryField() = pMixCathodeInit_;
		// Liquid Water Saturation
		sLiq.internalField() = sLiqCathodeInit_;
		sLiq.boundaryField() = sLiqCathodeInit_;
	}
	
	// Correct the boundary conditions
	forAll(APTLRegions, zoneID)
	{
		#include "APTLSetFields.H"
		// Electrons
		potElectron.correctBoundaryConditions();
		// Hydrogen
		wH2.correctBoundaryConditions();
		// Water Vapour
		wH2OVap.correctBoundaryConditions();
		// Temperature
		T.correctBoundaryConditions();
		// Mixture Pressure
		pMix.correctBoundaryConditions();
		// Liquid Water Saturation
		sLiq.correctBoundaryConditions();
	}
	forAll(ACLRegions, zoneID)
	{
		#include "ACLSetFields.H"
		// Electrons
		potElectron.correctBoundaryConditions();
		// Protons
		potProton.correctBoundaryConditions();
		// Hydrogen
		wH2.correctBoundaryConditions();
		// Water Vapour
		wH2OVap.correctBoundaryConditions();
		// Temperature
		T.correctBoundaryConditions();
		// Mixture Pressure
		pMix.correctBoundaryConditions();
		// Liquid Water Saturation
		sLiq.correctBoundaryConditions();
	}

	forAll(PMEMRegions, zoneID)
	{
		#include "PMEMSetFields.H"
		// Protons
		potProton.correctBoundaryConditions();
		// Temperature
		T.correctBoundaryConditions();
	}
	
	forAll(CCLRegions, zoneID)
	{
		#include "CCLSetFields.H"
		// Electrons
		potElectron.correctBoundaryConditions();
		// Protons
		potProton.correctBoundaryConditions();
		// Oxygen
		wO2.correctBoundaryConditions();
		// Water Vapour
		wH2OVap.correctBoundaryConditions();
		// Temperature
		T.correctBoundaryConditions();
		// Mixture Pressure
		pMix.correctBoundaryConditions();
		// Liquid Water Saturation
		sLiq.correctBoundaryConditions();
	}

	forAll(CPTLRegions, zoneID)
	{
		#include "CPTLSetFields.H"
		// Electrons
		potElectron.correctBoundaryConditions();
		// Oxygen
		wO2.correctBoundaryConditions();
		// Water Vapour
		wH2OVap.correctBoundaryConditions();
		// Temperature
		T.correctBoundaryConditions();
		// Mixture Pressure
		pMix.correctBoundaryConditions();
		// Liquid Water Saturation
		sLiq.correctBoundaryConditions();
	}

	// Write the initialization fields
	#include initWriteSolvedFields.H
}
