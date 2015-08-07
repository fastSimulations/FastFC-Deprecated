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
// Solve Electrons
{
	//Solve Left
	forAll(ACLRegions, zoneID)
	{
		#include "ACLSetFields.H"
   		#include "ACLReadSIMPLEControls.H"
		#include "ACLCalcElectrochemistry.H"
		sourceElectronTotal = reactionRateBVSum;
  		#include "electronTransportEquationA.H"
	}
	// Solve Right 
	forAll(CCLRegions, zoneID)
	{
		#include "CCLSetFields.H"
		#include "CCLReadSIMPLEControls.H"
		#include "CCLCalcElectrochemistry.H"
		sourceElectronTotal = reactionRateBVSum;
   		#include "electronTransportEquationC.H"
	}
}

			
// Solve Protons
{
	// Solve Left 
	forAll(ACLRegions, zoneID)
	{
    	#include "ACLSetFields.H"
	 	#include "ACLReadSIMPLEControls.H"
		#include "ACLCalcElectrochemistry.H"
		sourceProtonTotal = reactionRateBVSum;
		#include "protonTransportEquationA.H"
	}
	// Solve Centre 
	forAll(PMEMRegions, zoneID)
	{
   		#include "PMEMSetFields.H"
		#include "PMEMReadSIMPLEControls.H"
		sourceProtonTotal.internalField() = 0.0;
		sourceProtonTotal.boundaryField() = 0.0;
		#include "protonTransportEquationM.H"
	}
	// Solve Right 
	forAll(CCLRegions, zoneID)
	{
   		#include "CCLSetFields.H"
		#include "CCLReadSIMPLEControls.H"
		#include "CCLCalcElectrochemistry.H"
		sourceProtonTotal = reactionRateBVSum;
		#include "protonTransportEquationC.H"
	}
}
// ************************************************************************* //
