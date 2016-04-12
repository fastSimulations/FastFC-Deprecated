/*---------------------------------------------------------------------------*\
    F ast           	    | FAST-FC: 
	is the		            | The Open Source Analysis and Simulation Toolbox 
    A nalysis and           | for Fuel Cells
	S imulation		        |
	Toolbox for		        | Copyright 2016, David B. Harvey
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
	initialize the general and coefficient fields based on the initialized
	solution  

Notes

Developers
	David B. Harvey
\*---------------------------------------------------------------------------*/ 

// Initialize the General and Coefficient Fields

// APTL Field Calculations
forAll(APTLRegions, zoneID)
{
	#include "APTLSetFields.H"
	#include "APTLGeneralFields.H"
	#include "APTLCoefficientFields.H"
}

// ACL Field Calculations
forAll(ACLRegions, zoneID)
{
	#include "ACLSetFields.H"
	#include "ACLGeneralFields.H"
	#include "ACLCoefficientFields.H"
	#include "ACLElectrochemistryFields.H"
}

// CCL Field Calculations
forAll(CCLRegions,zoneID)
{
	#include "CCLSetFields.H"
	#include "CCLGeneralFields.H"
	#include "CCLCoefficientFields.H"
	#include "CCLElectrochemistryFields.H"
}

// CPTL Field Calculations
forAll(CPTLRegions,zoneID)
{
	#include "CPTLSetFields.H"
	#include "CPTLGeneralFields.H"
	#include "CPTLCoefficientFields.H"
}

// PMEM Field Calculations
forAll(PMEMRegions, zoneID)
{
	#include "PMEMSetFields.H"
	#include "PMEMGetRelativeHumidity.H"
	#include "PMEMGeneralFields.H"
	#include "PMEMCoefficientFields.H"
}

