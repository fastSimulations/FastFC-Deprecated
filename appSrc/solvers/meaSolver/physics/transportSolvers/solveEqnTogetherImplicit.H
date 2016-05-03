/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    F ast           	    | FAST-FC: 
	is the		            | The Open Source Analysis and Simulation Toolbox 
    A nalysis and           | for Fuel Cells
	S imulation		        |
	Toolbox for		        | Copyright 2015, David B. Harvey
	F uel                   |         
	C ells                  |                     
-------------------------------------------------------------------------------
License
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

	This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Application
    - solveEqnSeparateImplicit.H

Description
	- Solves the transport equations in the separate sub-region loops and 
	uses implicit descriptions for the source terms
:
\*---------------------------------------------------------------------------*/
// Calculated the changing properties and fields
#include <updateFields.H>

forAll(APTLRegions, zoneID)
{
	#include <APTLSetFields.H>
	if(1==solveHydrogen){
		#include <hydrogenTransportEquationAPTL.H>
		#include <APTLGeneralFields.H>
		#include <APTLCoefficientFields.H>
	}
	if(1==solveWaterVapour){
	  	#include <waterVapourTransportEquationAPTL.H>
		#include <APTLGeneralFields.H>
		#include <APTLCoefficientFields.H>
	}
	if(1==solveMass){
		#include <massTransportEquationAPTL.H>
		#include <APTLGeneralFields.H>
		#include <APTLCoefficientFields.H>
	}
	// solve electron transport
	{
  		#include <electronTransportEquationAPTL.H>
		#include <APTLGeneralFields.H>
		#include <APTLCoefficientFields.H>
	}
}
forAll(ACLRegions, zoneID)
{
	#include <ACLSetFields.H>
	#include <ACLReadSIMPLEControls.H>
	#include <ACLElectrochemistryFields.H>
	#include <currentSourceACL.H>

	if(1==solveHydrogen){
		#include <hydrogenTransportEquationACL.H>
		#include <ACLGeneralFields.H>
		#include <ACLCoefficientFields.H>
		#include <ACLElectrochemistryFields.H>
		#include <currentSourceACL.H>
	}

	if(1==solveWaterVapour){
		#include <waterVapourTransportEquationACL.H>
		#include <ACLGeneralFields.H>
		#include <ACLCoefficientFields.H>
		#include <ACLElectrochemistryFields.H>
		#include <currentSourceACL.H>
	}
	if(1==solveMass){
		#include <massTransportEquationACL.H>
		#include <ACLGeneralFields.H>
		#include <ACLCoefficientFields.H>
		#include <ACLElectrochemistryFields.H>
		#include <currentSourceACL.H>
	}
	// solve electron transport
	{
		#include <electronTransportEquationACL.H>
		#include <ACLGeneralFields.H>
		#include <ACLCoefficientFields.H>
		#include <ACLElectrochemistryFields.H>
		#include <currentSourceACL.H>
	}
	// solve proton transport	
	{
		#include <protonTransportEquationACL.H>
		#include <ACLGeneralFields.H>
		#include <ACLCoefficientFields.H>
		#include <ACLElectrochemistryFields.H>
		#include <currentSourceACL.H>
	}
}
forAll(PMEMRegions, zoneID)
{
	#include <PMEMSetFields.H>
	#include <PMEMReadSIMPLEControls.H>
	#include <PMEMGetRelativeHumidity.H>
	#include <protonTransportEquationPMEM.H>
	#include <PMEMGeneralFields.H>
	#include <PMEMCoefficientFields.H>
}
forAll(CCLRegions, zoneID)
{
	#include <CCLSetFields.H>
	#include <CCLReadSIMPLEControls.H>
	#include <CCLElectrochemistryFields.H>
	#include <currentSourceCCL.H>

	if(1==solveOxygen){
		#include <oxygenTransportEquationCCL.H>
		#include <CCLGeneralFields.H>
		#include <CCLCoefficientFields.H>
		#include <CCLElectrochemistryFields.H>
		#include <currentSourceCCL.H>
	}
	if(1==solveWaterVapour){
	 	#include <waterVapourTransportEquationCCL.H>
		#include <CCLGeneralFields.H>
		#include <CCLCoefficientFields.H>
		#include <CCLElectrochemistryFields.H>
		#include <currentSourceCCL.H>
	}
	if(1==solveMass){
		#include <massTransportEquationCCL.H>
		#include <CCLGeneralFields.H>
		#include <CCLCoefficientFields.H>
		#include <CCLElectrochemistryFields.H>
		#include <currentSourceCCL.H>
	}
	// solve electron transport
	{
 		#include <electronTransportEquationCCL.H>
		#include <CCLGeneralFields.H>
		#include <CCLCoefficientFields.H>
		#include <CCLElectrochemistryFields.H>
		#include <currentSourceCCL.H>
	}
	// solve proton transport
	{
		#include <protonTransportEquationCCL.H>
		#include <CCLGeneralFields.H>
		#include <CCLCoefficientFields.H>
		#include <CCLElectrochemistryFields.H>
		#include <currentSourceCCL.H>
	}

}
forAll(CPTLRegions, zoneID)
{
	#include <CPTLSetFields.H>
	if(1==solveOxygen){
		#include <oxygenTransportEquationCPTL.H>
		#include <CPTLGeneralFields.H>
		#include <CPTLCoefficientFields.H>
	}
	if(1==solveWaterVapour){
		#include <waterVapourTransportEquationCPTL.H>
		#include <CPTLGeneralFields.H>
		#include <CPTLCoefficientFields.H>
	}
	if(1==solveMass){
		#include <massTransportEquationCPTL.H>
		#include <CPTLGeneralFields.H>
		#include <CPTLCoefficientFields.H>
	}
	// solve electron transport
	{
		#include <electronTransportEquationCPTL.H>
		#include <CPTLGeneralFields.H>
		#include <CPTLCoefficientFields.H>
	}
}
// ************************************************************************* //
