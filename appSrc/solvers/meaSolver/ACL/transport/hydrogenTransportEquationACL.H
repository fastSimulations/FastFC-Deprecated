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
    - hydrogenTransportEquationACL.H

Description
	- Solves the hydrogen transport equation with explicit/implicit source 
	terms and calculates the slope using a numerically estimated derivative. 
:
\*---------------------------------------------------------------------------*/
{
	// Request previous iteration to be stored
	volScalarField wH2Old_ = wH2;

	// Create fields and steps for the slope of the volumetric current with 
	// hydrogen concentration
	dimensionedScalar stepSize_("stepSize_", dimensionSet(xH2.dimensions()), slopeStep);
	volScalarField xH2Left_ = xH2 - stepSize_;
	volScalarField xH2Right_ = xH2 + stepSize_;
	xH2Left_ = Foam::max(xH2Left_, Foam::VSMALL);

	volScalarField sourceVolHydrogenLeft_(sourceMassH2);
	volScalarField sourceVolHydrogenRight_(sourceMassH2);

	// Compute the numerical first derivative
	#include <hydrogenSourceSlopeACL.H>

	// Determine the slope
	volScalarField dSdPhi_ = (sourceVolHydrogenRight_ - sourceVolHydrogenLeft_)
					/(xH2Right_-xH2Left_);
	
	// Adapt the transport coefficient 
	transportCoeffH2 = densityMix*diffH2MixEff;	

	// Set the Oxygen Sink/Source Source Term
	#include <ACLH2SourceTermFields.H>

	// Calculate the electron transport equation

	fvScalarMatrix hydrogenTransport
	(
	 //   fvm::ddt(wH2)
	    -fvm::laplacian(transportCoeffH2,wH2) //+ fvm::SuSp(dSdPhi_,wH2)
		== 
		sourceMassH2 //+ dSdPhi_*wH2Old_
	);
	
	// Relax the Equation
//	hydrogenTransport.relax();

	// Solve the Matrix System
	hydrogenTransport.solve();
	
	// Relax the solved variable
	scalar relax_ = relaxHydrogen;
	volScalarField deltaH2_ = wH2 - wH2Old_;
	wH2 = wH2Old_ + deltaH2_*relax_;
	
	// Minimization Check
	{
		dimensionedScalar wH2Min_
		(
		 	"WH2Min_",
			dimensionSet(wH2.dimensions()),
			Foam::VSMALL
		);

		wH2 = Foam::max(wH2,wH2Min_);
	}	

	// Correct the boundary conditions on the solved variable
	wH2.correctBoundaryConditions();
}
