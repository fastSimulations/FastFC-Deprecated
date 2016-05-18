{
	//Solve Anode
	forAll(APTLRegions, zoneID)
	{
		#include <APTLSetFields.H>
   		#include <APTLReadSIMPLEControls.H>
  		#include <electronTransportEquationAPTL.H>
		#include <APTLGeneralFields.H>
		#include <APTLCoefficientFields.H>
	}
	forAll(ACLRegions, zoneID)
	{
		#include <ACLSetFields.H>
 		#include <ACLReadSIMPLEControls.H>
		#include <electronTransportEquationACL.H>
		#include <ACLGeneralFields.H>
		#include <ACLCoefficientFields.H>
	}
	// Solve Cathode  
	forAll(CCLRegions, zoneID)
	{
		#include <CCLSetFields.H>
		#include <CCLReadSIMPLEControls.H>
 		#include <electronTransportEquationCCL.H>
		#include <CCLGeneralFields.H>
		#include <CCLCoefficientFields.H>
	}
	forAll(CPTLRegions, zoneID)
	{
		#include <CPTLSetFields.H>
		#include <CPTLReadSIMPLEControls.H>
 		#include <electronTransportEquationCPTL.H>
		#include <CPTLGeneralFields.H>
		#include <CPTLCoefficientFields.H>
	}
}
