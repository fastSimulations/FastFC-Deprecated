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
	Master file for setup of the initial variables and solution guess 

Notes

Developers
	David B. Harvey
\*---------------------------------------------------------------------------*/ 
// Determine the active area and electrode thicknesses
#include "initElectrodeGeometry.H"

// Setup the transport coefficients to dummy coefficients (for smooth initial guess)
// This will be redundant on move to the new BC's
#include "initTransportCoefficients.H"

// Setup the solved variable profiles and write to the 0/$COMPONENT directory
#include "initSolvedFields.H"	

// Setup General and Coefficient Fields based on the solved variable init
#include "initGeneralCoefficientFields.H"

// Correct Solved Boundary Conditions
#include "initCorrectSolvedBoundaryConditions.H"

// Write the initialization fields
#include "initWriteSolvedFields.H"
