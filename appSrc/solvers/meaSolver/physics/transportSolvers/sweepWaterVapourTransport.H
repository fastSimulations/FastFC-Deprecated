{
	//Solve Anode
	forAll(APTLRegions, zoneID)
	{
		#include "APTLSetFields.H"
   		#include "APTLReadSIMPLEControls.H"
  		#include "waterVapourTransportEquationAPTL.H"
		#include <APTLGeneralFields.H>
		#include <APTLCoefficientFields.H>
	}
	forAll(ACLRegions, zoneID)
	{
		#include "ACLSetFields.H"
 		#include "ACLReadSIMPLEControls.H"
 		#include "waterVapourTransportEquationACL.H"
		#include <ACLGeneralFields.H>
		#include <ACLCoefficientFields.H>
	}
	// Solve Cathode  
	forAll(CCLRegions, zoneID)
	{
		#include "CCLSetFields.H"
		#include "CCLReadSIMPLEControls.H"
 		#include "waterVapourTransportEquationCCL.H"
		#include <CCLGeneralFields.H>
		#include <CCLCoefficientFields.H>
	}
	forAll(CPTLRegions, zoneID)
	{
		#include "CPTLSetFields.H"
		#include "CPTLReadSIMPLEControls.H"
 		#include "waterVapourTransportEquationCPTL.H"
		#include <CPTLGeneralFields.H>
		#include <CPTLCoefficientFields.H>
	}
}
