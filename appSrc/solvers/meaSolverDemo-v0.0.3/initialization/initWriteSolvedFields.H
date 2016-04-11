	// Write the initialization fields
	forAll(APTLRegions, zoneID)
	{
		#include "APTLSetFields.H"
		// Electrons
		potElectron.write();
		// Hydrogen
		wH2.write();
		// Water Vapour
		wH2OVap.write();
		// Temperature
		T.write();
		// Mixture Pressure
		pMix.write();
		// Liquid Water Saturation
		sLiq.write();
	}
	forAll(ACLRegions, zoneID)
	{
		#include "ACLSetFields.H"
		// Electrons
		potElectron.write();
		// Protons
		potProton.write();
		// Hydrogen
		wH2.write();
		// Water Vapour
		wH2OVap.write();
		// Temperature
		T.write();
		// Mixture Pressure
		pMix.write();
		// Liquid Water Saturation
		sLiq.write();
	}

	forAll(PMEMRegions, zoneID)
	{
		#include "PMEMSetFields.H"
		// Protons
		potProton.write();
		// Temperature
		T.write();
	}
	
	forAll(CCLRegions, zoneID)
	{
		#include "CCLSetFields.H"
		// Electrons
		potElectron.write();
		// Protons
		potProton.write();
		// Oxygen
		wO2.write();
		// Water Vapour
		wH2OVap.write();
		// Temperature
		T.write();
		// Mixture Pressure
		pMix.write();
		// Liquid Water Saturation
		sLiq.write();
	}

	forAll(CPTLRegions, zoneID)
	{
		#include "CPTLSetFields.H"
		// Electrons
		potElectron.write();
		// Oxygen
		wO2.write();
		// Water Vapour
		wH2OVap.write();
		// Temperature
		T.write();
		// Mixture Pressure
		pMix.write();
		// Liquid Water Saturation
		sLiq.write();
	}

