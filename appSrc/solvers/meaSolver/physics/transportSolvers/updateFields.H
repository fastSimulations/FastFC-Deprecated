forAll(APTLRegions, zoneID)
{
	#include <APTLSetFields.H>
	#include <APTLGeneralFields.H>
	#include <APTLCoefficientFields.H>
}
forAll(ACLRegions, zoneID)
{
	#include <ACLSetFields.H>
	#include <ACLGeneralFields.H>
	#include <ACLCoefficientFields.H>
}
forAll(PMEMRegions, zoneID)
{
	#include <PMEMSetFields.H>
	#include <PMEMGetRelativeHumidity.H>
	#include <PMEMGeneralFields.H>
	#include <PMEMCoefficientFields.H>
}
forAll(CCLRegions, zoneID)
{
	#include <CCLSetFields.H>
	#include <CCLGeneralFields.H>
	#include <CCLCoefficientFields.H>
}
forAll(CPTLRegions, zoneID)
{
	#include <CPTLSetFields.H>
	#include <CPTLGeneralFields.H>
	#include <CPTLCoefficientFields.H>
}
