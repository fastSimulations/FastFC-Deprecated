{
	//Solve Anode
	forAll(ACLRegions, zoneID)
	{
		#include <ACLSetFields.H>
 		#include <ACLReadSIMPLEControls.H>
		#include <ACLElectrochemistryFields.H>
		#include <currentSourceACL.H>
	}
	// Solve Cathode  
	forAll(CCLRegions, zoneID)
	{
		#include <CCLSetFields.H>
		#include <CCLReadSIMPLEControls.H>
		#include <CCLElectrochemistryFields.H>
		#include <currentSourceCCL.H>
	}
}
