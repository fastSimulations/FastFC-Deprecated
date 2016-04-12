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
    - electronTransportEquationACL.H

Description
	- Solves the electron transport equation with explicit/implicit source 
	terms in the Anode Catalyst Layer and calculates the slope using a 
	numerical estimated derivative. 
:
\*---------------------------------------------------------------------------*/
{
	//  Electronic Transport
	// Request previous iteration to be stored
	// potElectron.storePrevIter();
	volScalarField potElectronOld_ = potElectron;

	// Create fields and steps for the slope of the volumetric current with electron potential
	dimensionedScalar stepSize_("stepSize_", dimensionSet(actPot.dimensions()), slopeStep);
	volScalarField actPotLeft_ = actPot - stepSize_;
	volScalarField actPotRight_ = actPot + stepSize_;
	volScalarField sourceVolCurrentLeft_(sourceVolCurrent);
	volScalarField sourceVolCurrentRight_(sourceVolCurrent);
	
	// Compute the current and numerical first derivative	
	#include <electronSourceSlopeACL.H>

	// Determine the slope
	volScalarField dSdPhi_ = (sourceVolCurrentRight_ - sourceVolCurrentLeft_)
					/(actPotRight_-actPotLeft_);

	// Set the Electron Source Term
	#include <ACLElectronSourceTermFields.H>

	// Calculate the electron transport equation

	fvScalarMatrix electronicTransport
	(
	 //   fvm::ddt(phi_e)
	    -fvm::laplacian(condElectronEff,potElectron) + fvm::SuSp(dSdPhi_,potElectron)
		== 
		dSdPhi_*potElectronOld_ - sourceElectronTotal 
	);

	// Relax the Equation
//	electronicTransport.relax();
	// Store the flux (may not be needed, TBD)
//	electronFlux = electronicTransport.flux();
	// Solve the Matrix System
	electronicTransport.solve();
	// Relax the solved variable
	scalar relax_ = relaxElectrons;
	volScalarField deltaElectron_ = potElectron - potElectronOld_;
	potElectron = potElectronOld_ + deltaElectron_*relax_;
//	potElectron.relax();
	// Correct the boundary conditions on the solved variable
	potElectron.correctBoundaryConditions();
}
