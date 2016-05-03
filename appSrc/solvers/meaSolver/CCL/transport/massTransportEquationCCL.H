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
    - massTransportEquationCCL.H

Description
	- Solves the mass transport equation, recast for mixture pressure,
	with explicit source terms. 
:
\*---------------------------------------------------------------------------*/
{
	// Request previous iteration to be stored
	volScalarField pMixOld_ = pMix;

	dynamicPerm = densityMix*kPermRelGas/viscMix;

	// Set the Mass Sink/Source Source Term
	#include <CCLMassSourceTermFields.H>

	// Assemble the mass transport matrix
	fvScalarMatrix massTransport
	(
	    -fvm::laplacian(dynamicPerm,pMix) 
		==
		sourceMassTotal
	);
	
	// Solve the Matrix System
	massTransport.solve();

	// Relax the solved variable
	scalar relax_ = relaxMass;
	volScalarField pMixDelta_ = pMix - pMixOld_;
	pMix = pMixOld_ + pMixDelta_*relax_;
	
	// Minimization Check
	{
		dimensionedScalar pMixMin_
		(
		 	"pMixMin_",
			dimensionSet(pMix.dimensions()),
			Foam::VSMALL
		);

		pMix = Foam::max(pMix,pMixMin_);
	}	

	// Correct the boundary conditions on the solved variable
	pMix.correctBoundaryConditions();
}
