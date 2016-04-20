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
    - oxygenTransportEquationCPTL.H

Description
	- Solves the oxygen transport equation 
:
\*---------------------------------------------------------------------------*/
{
	// Request previous iteration to be stored
	//	wO2.storePrevIter();
	volScalarField wO2Old_ = wO2;

	transportCoeffO2 = densityMix*diffO2MixEff;

	fvScalarMatrix oxygenTransport
	(
	//	fvm::ddt(wO2)
	    -fvm::laplacian(transportCoeffO2,wO2) 
	);

	// Relax the Equation
//	oxygenTransport.relax();

	// Solve the Matrix System
	oxygenTransport.solve();

	// Relax the solved variable
	scalar relax_ = relaxOxygen;
	volScalarField deltaO2_ = wO2 - wO2Old_;
	wO2 = wO2Old_ + deltaO2_*relax_;

	// Minimization Check
	{
		dimensionedScalar wO2Min_
		(
		 	"WO2Min_",
			dimensionSet(wO2.dimensions()),
			Foam::VSMALL
		);

		wO2 = Foam::max(wO2,wO2Min_);
	}	

	// Correct the boundary conditions on the solved variable
	wO2.correctBoundaryConditions();
}
