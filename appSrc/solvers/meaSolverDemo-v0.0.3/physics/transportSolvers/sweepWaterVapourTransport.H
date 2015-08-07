{
	//Solve Anode
	forAll(APTLRegions, zoneID)
	{
		#include "APTLSetFields.H"
   		#include "APTLReadSIMPLEControls.H"
  		#include "waterVapourTransportEquationAPTL.H"
	}
	forAll(ACLRegions, zoneID)
	{
		#include "ACLSetFields.H"
 		#include "ACLReadSIMPLEControls.H"
		#include "ACLElectrochemistryFields.H"
 		#include "waterVapourTransportEquationACL.H"
	}
	// Solve Cathode  
	forAll(CCLRegions, zoneID)
	{
		#include "CCLSetFields.H"
		#include "CCLReadSIMPLEControls.H"
		#include "CCLElectrochemistryFields.H"
 		#include "waterVapourTransportEquationCCL.H"
	}
	forAll(CPTLRegions, zoneID)
	{
		#include "CPTLSetFields.H"
		#include "CPTLReadSIMPLEControls.H"
 		#include "waterVapourTransportEquationCPTL.H"
	}
}
