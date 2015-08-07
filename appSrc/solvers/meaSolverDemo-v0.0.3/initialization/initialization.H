/*---------------------------------------------------------------------------*\
    F ast           	    | FAST-FC: 
	is the		            | The Open Source Analysis and Simulation Toolbox 
    A nalysis and           | for Fuel Cells
	S imulation		        |
	Toolbox for		        | Copyright 2015, David B. Harvey
	F uel                   |         
	C ells                  |                     
-------------------------------------------------------------------------------
License
	FAST-FC and this file are a derivative work of OpenFOAM.

	FAST-FC is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    FAST-FC is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with FAST-FC.  If not, see <http://www.gnu.org/licenses/>.

Application
   meaSolverDemo 

File
	initialization.H

Description
	Master file for setup of the initial variables and solution guess 

Notes
	1/ 	The electronic potential on the cathode uses the initial guess 
		in the change dict file.
	2/	The Transport Coefficients are set to contain "1" which aids in setup of smooth
		solved variable profiles
	3/	Current production is assumed to be zero (eta=0)
	4/	Cathode electrical potential is equal to the set voltage (changeDict)
	5/	Anode electrical potential is set to zero
	6/	Cathode and anode proton potential are derived to keep (eta=0)
	7/	Membrane proton potential is set through a matrix solve

Developers
	David B. Harvey
\*---------------------------------------------------------------------------*/ 
// Determine the active area and electrode thicknesses
#include "initElectrodeGeometry.H"

// Setup the transport coefficients to dummy coefficients (for smooth initial guess)
#include "initTransportCoefficients.H"

// Setup the solved variable profiles and write to the 0/$COMPONENT directory
#include "initSolvedFields.H"	


forAll(APTLRegions, zoneID)
{
	#include "APTLSetFields.H"
	#include "APTLGeneralFields.H"
	#include "APTLCoefficientFields.H"
}
forAll(ACLRegions, zoneID)
{
	#include "ACLSetFields.H"
	#include "ACLGeneralFields.H"
	#include "ACLCoefficientFields.H"
	#include "ACLElectrochemistryFields.H"
}
forAll(CCLRegions,zoneID)
{
	#include "CCLSetFields.H"
	#include "CCLGeneralFields.H"
	#include "CCLCoefficientFields.H"
	#include "CCLElectrochemistryFields.H"
}
forAll(CPTLRegions,zoneID)
{
	#include "CPTLSetFields.H"
	#include "CPTLGeneralFields.H"
	#include "CPTLCoefficientFields.H"
}
forAll(PMEMRegions, zoneID)
{
	#include "PMEMSetFields.H"
	#include "PMEMGetRelativeHumidity.H"
	#include "PMEMGeneralFields.H"
	#include "PMEMCoefficientFields.H"
}
