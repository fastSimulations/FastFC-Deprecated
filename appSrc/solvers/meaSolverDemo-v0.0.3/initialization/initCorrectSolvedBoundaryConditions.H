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
	initCorrectSolvedBoundaryConditions.H

Description
	Apply the boundary conditions for the solved variables based on the 
	initialized solution

Notes

Developers
	David B. Harvey
\*---------------------------------------------------------------------------*/ 

// Correct the boundary conditions
forAll(APTLRegions, zoneID)
{
	#include "APTLSetFields.H"
	// Electrons
	potElectron.correctBoundaryConditions();
	// Hydrogen
	wH2.correctBoundaryConditions();
	// Water Vapour
	wH2OVap.correctBoundaryConditions();
	// Temperature
	T.correctBoundaryConditions();
	// Mixture Pressure
	pMix.correctBoundaryConditions();
	// Liquid Water Saturation
	sLiq.correctBoundaryConditions();
}

forAll(ACLRegions, zoneID)
{
	#include "ACLSetFields.H"
	// Electrons
	potElectron.correctBoundaryConditions();
	// Protons
	potProton.correctBoundaryConditions();
	// Hydrogen
	wH2.correctBoundaryConditions();
	// Water Vapour
	wH2OVap.correctBoundaryConditions();
	// Temperature
	T.correctBoundaryConditions();
	// Mixture Pressure
	pMix.correctBoundaryConditions();
	// Liquid Water Saturation
	sLiq.correctBoundaryConditions();
}

forAll(PMEMRegions, zoneID)
{
	#include "PMEMSetFields.H"
	// Protons
	potProton.correctBoundaryConditions();
	// Temperature
	T.correctBoundaryConditions();
}
	
forAll(CCLRegions, zoneID)
{
	#include "CCLSetFields.H"
	// Electrons
	potElectron.correctBoundaryConditions();
	// Protons
	potProton.correctBoundaryConditions();
	// Oxygen
	wO2.correctBoundaryConditions();
	// Water Vapour
	wH2OVap.correctBoundaryConditions();
	// Temperature
	T.correctBoundaryConditions();
	// Mixture Pressure
	pMix.correctBoundaryConditions();
	// Liquid Water Saturation
	sLiq.correctBoundaryConditions();
}

forAll(CPTLRegions, zoneID)
{
	#include "CPTLSetFields.H"
	// Electrons
	potElectron.correctBoundaryConditions();
	// Oxygen
	wO2.correctBoundaryConditions();
	// Water Vapour
	wH2OVap.correctBoundaryConditions();
	// Temperature
	T.correctBoundaryConditions();
	// Mixture Pressure
	pMix.correctBoundaryConditions();
	// Liquid Water Saturation
	sLiq.correctBoundaryConditions();
}

