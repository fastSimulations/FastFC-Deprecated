{
	//Solve Anode
	forAll(APTLRegions, zoneID)
	{
		#include <APTLSetFields.H>
   		#include <APTLReadSIMPLEControls.H>
  		#include <electronTransportEquationAPTL.H>
	}
	forAll(ACLRegions, zoneID)
	{
		#include <ACLSetFields.H>
 		#include <ACLReadSIMPLEControls.H>
		#include <electronTransportEquationACL.H>
	}
	// Solve Cathode  
	forAll(CCLRegions, zoneID)
	{
		#include <CCLSetFields.H>
		#include <CCLReadSIMPLEControls.H>
 		#include <electronTransportEquationCCL.H>
	}
	forAll(CPTLRegions, zoneID)
	{
		#include <CPTLSetFields.H>
		#include <CPTLReadSIMPLEControls.H>
 		#include <electronTransportEquationCPTL.H>
	}
}
