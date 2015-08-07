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

// Solve Electrochemistry
{
	#include <sweepElectrochemistry.H>
}
// Solve Hydrogen Transport
if(1==solveHydrogen)
{
	#include <sweepHydrogenTransport.H>
}
// Solve Oxygen Transport
if(1==solveOxygen)
{
	#include <sweepOxygenTransport.H>
}
// Solve Water Vapour Transport and Phase Change (if activated)
if(1==solveWaterVapour)
{
	#include <sweepWaterVapourTransport.H>
}	

// Solve Electrons and Protons
{
	// Electrons
	#include <sweepElectronTransport.H>
	// Protons
	#include <sweepProtonTransport.H>
}

// Solve Water Vapour Transport
//if(1==solveWaterVapour)
//{
//	#include "sweepWaterVapourTransport.H"
//}
// Solve Mass Transport
//if(1==solveMass)
//{
//	#include "sweepMassTransport.H"
//}
// Solve Energy Transport
//if(1==solveEnergy)
//{
//	#include "sweepEnergyTransport.H"
//}
// 	Solve Liquid Water Transport
//if(1==solveLiqWater)
//{
//	#include "sweepWaterLiquidTransport.H"
//}
// ************************************************************************* //
