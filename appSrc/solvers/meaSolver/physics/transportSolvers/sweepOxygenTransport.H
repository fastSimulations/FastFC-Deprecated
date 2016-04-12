{
	// Solve CPTL
	forAll(CPTLRegions, zoneID)
	{
		#include <CPTLSetFields.H>
		#include <CPTLReadSIMPLEControls.H>
		#include <oxygenTransportEquationCPTL.H>
	}

	// Solve CCL
	forAll(CCLRegions, zoneID)
	{
		#include <CCLSetFields.H>
		#include <CCLReadSIMPLEControls.H>
		#include <oxygenTransportEquationCCL.H>
	}
}
