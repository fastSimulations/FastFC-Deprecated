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

	// Determine the Open Circuit Voltage
	{

		// Initialization of dimensionedScalars
		dimensionedScalar equilibriumPotential_
		(
			"equilibriumPotential_",
			dimensionSet(OCV.dimensions()),
			1.189
		);

		dimensionedScalar oxideOCV_
		(
			"oxideOCV_",
			dimensionSet(OCV.dimensions()),
			0.066349
		);

		dimensionedScalar nernstOCV_
		(
			"nernstOCV_",
			dimensionSet(OCV.dimensions()),
			VSMALL
		);

		dimensionedScalar satPCoeff1_
		(
			"satPCoeff1_",
			dimensionSet(Fast::Constants::pStd.dimensions()),
			1.26837e-3
		);

		dimensionedScalar satPCoeff2_
		(
			"satPCoeff2_",
			dimensionSet(Fast::Constants::pStd.dimensions()),
			1.49827
		);

		dimensionedScalar satPCoeff3_
		(
			"satPCoeff3_",
			dimensionSet(Fast::Constants::pStd.dimensions()),
			6.70916e2
		);

		dimensionedScalar satPCoeff4_
		(
			"satPCoeff4_",
			dimensionSet(Fast::Constants::pStd.dimensions()),
			1.34832e5
		);
		
		dimensionedScalar satPCoeff5_
		(
			"satPCoeff5_",
			dimensionSet(Fast::Constants::pStd.dimensions()),
			1.02503e7
		);
	
		dimensionedScalar tempTempAnode_
		(
			"tempTempAnode_",
			dimensionSet( 0, 0, 0, 1, 0, 0, 0),
			TAnodeInit_
		);

		dimensionedScalar tempTempCathode_ 
		(
			"tempTempCathode_",
			dimensionSet( 0, 0, 0, 1, 0, 0, 0),
			TCathodeInit_
		);
 
		dimensionedScalar pMixWetAnode_
		(
			"pMixWetAnode_",
			dimensionSet(Fast::Constants::pStd.dimensions()),
			pMixAnodeInit_
		);

		dimensionedScalar pMixWetCathode_
		(
			"pMixWetCathode_",
			dimensionSet(Fast::Constants::pStd.dimensions()),
			pMixCathodeInit_
		);

		dimensionedScalar pSatAnode_ 
		(
			"pSatAnode_",
			dimensionSet(Fast::Constants::pStd.dimensions()),
			VSMALL	
		);

		dimensionedScalar pSatCathode_ 
		(
			"pSatCathode_",
			dimensionSet(Fast::Constants::pStd.dimensions()),
			VSMALL
		);
		
		dimensionedScalar xH2AnodeInit_
		(
			"xH2AnodeInit_",
			dimensionSet( 0, 0, 0, 0, 0, 0, 0),
			VSMALL
		);

		dimensionedScalar xO2CathodeInit_
		(
			"xO2CathodeInit_",
			dimensionSet( 0, 0, 0, 0, 0, 0, 0),
			VSMALL
		);
		
		scalar wN2AnodeInit_ = VSMALL;
		scalar wN2CathodeInit_ = VSMALL;

		// Assign Values to dimensioned and scalar fields
		wN2AnodeInit_ = Foam::max((1. - wH2AnodeInit_ - wH2OVapAnodeInit_), VSMALL);
		wN2CathodeInit_ = Foam::max((1. - wO2CathodeInit_ - wH2OVapCathodeInit_), VSMALL);

		xH2AnodeInit_ = 
			(wH2AnodeInit_/Fast::Constants::molWeightH2.value())/
			(
				(wH2AnodeInit_/Fast::Constants::molWeightH2.value())
				+(wN2AnodeInit_/Fast::Constants::molWeightN2.value())
				+(wH2OVapAnodeInit_/Fast::Constants::molWeightH2O.value())
			);
		
		xO2CathodeInit_ =
			(wO2CathodeInit_/Fast::Constants::molWeightO2.value())
			/(
				(wO2CathodeInit_/Fast::Constants::molWeightO2.value())
				+(wN2CathodeInit_/Fast::Constants::molWeightN2.value())
				+(wH2OVapCathodeInit_/Fast::Constants::molWeightH2O.value())
			);

		pSatAnode_ = 
			(
				satPCoeff1_*Foam::pow(TAnodeInit_,4.)-satPCoeff2_*Foam::pow(TAnodeInit_,3.)
				+satPCoeff3_*Foam::pow(TAnodeInit_,2.)-satPCoeff4_*TAnodeInit_+satPCoeff5_
			);

		pSatCathode_ =
			(
				satPCoeff1_*Foam::pow(TCathodeInit_,4.)-satPCoeff2_*Foam::pow(TCathodeInit_,3.)
				+satPCoeff3_*Foam::pow(TCathodeInit_,2.)-satPCoeff4_*TCathodeInit_+satPCoeff5_
			);

		// Measured Oxide effect on OCV in the DOE Program STC Hardware
		
		nernstOCV_ = equilibriumPotential_ +
			(
				(Fast::Constants::gasConstantR*(tempTempAnode_+tempTempCathode_)/2.)
				/(2.*Fast::Constants::faraday)
			)*Foam::log
			(
				(
					(pMixWetAnode_*xH2AnodeInit_)/Fast::Constants::pStd
				)*
				Foam::pow
				(
					(
						(pMixWetCathode_*xO2CathodeInit_)/Fast::Constants::pStd
					),
					0.5
				)
			);
				

		OCV = (nernstOCV_ - oxideOCV_);

//		Info<< nl
//			<<	tab << "Time = " << runTime.timeName() << nl
//			<<	tab	<< "OCV must equal 1.11012229246918" << nl
//			<<	tab	<< "Open Circuit Voltage = " << OCV.value() << nl
//			<<	tab	<< "Nernst OCV = " << nernstOCV_.value() << nl
//			<<	tab	<< "Nernst Oxide Effect = " << oxideOCV_.value() << nl
//			<<	tab	<< "PmixWetAnode_ = " << pMixWetAnode_ << nl
//			<<	tab	<< "pSatAnode_ = " << pSatAnode_ << nl
//			<<	tab	<< "xH2AnodeInit_ = " << xH2AnodeInit_ << nl
//			<<	tab	<< "pMixWetCathode_ = " << pMixWetCathode_ << nl
//			<<	tab	<< "pSatCathode_ = " << pSatCathode_ << nl
//			<<	tab	<< "xO2CathodeInit_ = " << xO2CathodeInit_ << nl
//			<< endl;
	}

	// Output the Averaged Settings
	{
		Info<< nl
			<<	tab << "Time = " << runTime.timeName() << nl
			<<	tab	<< "Open Circuit Voltage = " << OCV.value() << nl
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
}
