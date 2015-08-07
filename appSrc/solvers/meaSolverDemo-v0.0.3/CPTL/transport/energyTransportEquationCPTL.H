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
    - waterVapourTransportEquationCPTL.H

Description
	- Solves the transport equation. 
:
\*---------------------------------------------------------------------------*/
{
	//  Protonic Transport
	// Request previous iteration to be stored
	T.storePrevIter();

	// Adapt the transport coefficient 
	transportCoeffEnergy = condThermalEff;	

	// Set the Oxygen Sink/Source Source Term
	#include "CPTLEnergySourceTermFields.H"

	fvScalarMatrix waterVapourTransport
	(
	//	fvm::ddt(densityMix*T)
	    -fvm::laplacian(transportCoeffEnergy,T) 
		==
		sourceEnergyTotal
	);

	// Relax the Equation
	waterVapourTransport.relax();
	// Store the flux (may not be needed, TBD)
//	waterVapourFlux = waterVapourTransport.flux();
	// Solve the Matrix System
	waterVapourTransport.solve();
	// Relax the solved variable
	T.relax();
	// Minimization Check
	{
		dimensionedScalar TMin_
		(
		 	"TMin_",
			dimensionSet(T.dimensions()),
			273.
		);

		T = Foam::max(T,TMin_);
	}	
	// Correct the boundary conditions on the solved variable
	T.correctBoundaryConditions();
}
