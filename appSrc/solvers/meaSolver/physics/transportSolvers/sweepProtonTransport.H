{
	// Solve Left 
	forAll(ACLRegions, zoneID)
	{
   		#include <ACLSetFields.H>
	 	#include <ACLReadSIMPLEControls.H>
		#include <protonTransportEquationACL.H>
		#include <ACLGeneralFields.H>
		#include <ACLCoefficientFields.H>
	}
	// Solve Centre 
	forAll(PMEMRegions, zoneID)
	{
 		#include <PMEMSetFields.H>
		#include <PMEMReadSIMPLEControls.H>
		#include <PMEMGetRelativeHumidity.H>
		#include <protonTransportEquationPMEM.H>
		#include <PMEMGeneralFields.H>
		#include <PMEMCoefficientFields.H>
	}
	// Solve Right 
	forAll(CCLRegions, zoneID)
	{
   		#include <CCLSetFields.H>
		#include <CCLReadSIMPLEControls.H>
		#include <protonTransportEquationCCL.H>
		#include <CCLGeneralFields.H>
		#include <CCLCoefficientFields.H>
	}
	// Solve Centre 
	forAll(PMEMRegions, zoneID)
	{
 		#include <PMEMSetFields.H>
		#include <PMEMReadSIMPLEControls.H>
		#include <PMEMGetRelativeHumidity.H>
		#include <protonTransportEquationPMEM.H>
		#include <PMEMGeneralFields.H>
		#include <PMEMCoefficientFields.H>
	}
	// Solve Left 
	forAll(ACLRegions, zoneID)
	{
   		#include <ACLSetFields.H>
	 	#include <ACLReadSIMPLEControls.H>
		#include <protonTransportEquationACL.H>
		#include <ACLGeneralFields.H>
		#include <ACLCoefficientFields.H>
	}
}
