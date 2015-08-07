{
	// Solve Left 
	forAll(ACLRegions, zoneID)
	{
   		#include <ACLSetFields.H>
	 	#include <ACLReadSIMPLEControls.H>
		#include <protonTransportEquationACL.H>
	}
	// Solve Centre 
	forAll(PMEMRegions, zoneID)
	{
 		#include <PMEMSetFields.H>
		#include <PMEMReadSIMPLEControls.H>
		#include <protonTransportEquationPMEM.H>
	}
	// Solve Right 
	forAll(CCLRegions, zoneID)
	{
   		#include <CCLSetFields.H>
		#include <CCLReadSIMPLEControls.H>
		#include <protonTransportEquationCCL.H>
	}
}


