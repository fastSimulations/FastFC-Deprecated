{
	//Solve Anode
	forAll(APTLRegions, zoneID)
	{
		#include "APTLSetFields.H"
   		#include "APTLReadSIMPLEControls.H"
	  	#include "solverPhysics/transportEquations/energyTransport/energyTransportEquationAPTL.H"
	}
	forAll(ACLRegions, zoneID)
	{
		#include "ACLSetFields.H"
 		#include "ACLReadSIMPLEControls.H"
 		#include "solverPhysics/transportEquations/energyTransportenergyTransportEquationACL.H"
	}
	// Solve membrane
	forAll(PMEMRegions, zoneID)
	{
		#include "PMEMSetFields.H"
 		#include "PMEMReadSIMPLEControls.H"
		#include "solverPhysics/transportEquations/energyTransportenergyTransportEquationPMEM.H"
	}

	// Solve Cathode  
	forAll(CCLRegions, zoneID)
	{
		#include "CCLSetFields.H"
		#include "CCLReadSIMPLEControls.H"
		#include "solverPhysics/transportEquations/energyTransportenergyTransportEquationCCL.H"
	}
	forAll(CPTLRegions, zoneID)
	{
		#include "CPTLSetFields.H"
		#include "CPTLReadSIMPLEControls.H"
		#include "solverPhysics/transportEquations/energyTransportenergyTransportEquationCPTL.H"
	}
}

