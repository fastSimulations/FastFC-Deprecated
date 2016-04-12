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

