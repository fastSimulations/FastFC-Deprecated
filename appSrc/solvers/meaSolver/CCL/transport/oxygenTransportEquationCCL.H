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
    - oxygenTransportEquationCCL.H

Description
	- Solves the oxygen transport equation with explicit/implicit source 
	terms and calculates the slope using a numerically estimated derivative. 
:
\*---------------------------------------------------------------------------*/
{
	// Store value of previous iteration
	volScalarField wO2Old_ = wO2;

	// Adapt the transport coefficient 
	transportCoeffO2 = densityMix*diffO2MixEff;	

	// Create fields and steps for the slope of the volumetric current with proton potential
	dimensionedScalar stepSize_("stepSize_", dimensionSet(xO2.dimensions()), slopeStep);
	volScalarField wO2Left_ = wO2 - stepSize_;
	volScalarField wO2Right_ = wO2 + stepSize_;
	volScalarField wN2Left_ = 1. - wO2Left_ - wH2OVap;
	volScalarField wN2Right_ = 1. - wO2Right_ - wH2OVap;
	volScalarField xO2Left_ = (wO2Left_/Fast::Constants::molWeightO2)/((wO2Left_/Fast::Constants::molWeightO2)
			+(wN2Left_/Fast::Constants::molWeightN2)+(wH2OVap/Fast::Constants::molWeightH2O));
	volScalarField xO2Right_ = (wO2Right_/Fast::Constants::molWeightO2)/((wO2Right_/Fast::Constants::molWeightO2)
			+(wN2Right_/Fast::Constants::molWeightN2)+(wH2OVap/Fast::Constants::molWeightH2O));

	volScalarField sourceVolOxygenLeft_(sourceMassO2);
	volScalarField sourceVolOxygenRight_(sourceMassO2);

	// Compute the current and numerical first derivative
	#include <oxygenSourceSlopeCCL.H>

	// Determine the slope
	volScalarField dSdPhi_ = (sourceVolOxygenRight_ - sourceVolOxygenLeft_)
					/(wO2Right_-wO2Left_);

	// Set the Oxygen Sink/Source Source Term
	#include <CCLOxygenSourceTermFields.H>
	
	// Calculate the oxygen transport equation
	fvScalarMatrix oxygenTransport
	(
	 //   fvm::ddt(wO2)
	    -fvm::laplacian(transportCoeffO2,wO2) + fvm::SuSp(dSdPhi_,wO2)
		== 
		dSdPhi_*wO2Old_ + sourceMassO2
	);
	
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
