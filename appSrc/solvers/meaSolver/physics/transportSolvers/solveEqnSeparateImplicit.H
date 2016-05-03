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
	Toolbox for		        | Copyright 2016, David B. Harvey
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

// Calculate the fields and coefficients
#include <updateFields.H>

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

// Solve Electron Transport 
#include <sweepElectronTransport.H>
// Solve Proton Transport
#include <sweepProtonTransport.H>
// Electrochemistry
#include <sweepElectrochemistry.H>

// ************************************************************************* //
