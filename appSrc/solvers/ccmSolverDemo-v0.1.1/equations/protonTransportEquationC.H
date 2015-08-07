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
    - protonTransportEquationC.H

Description
	- Solves the proton transport equation with explicit/implicit terms 
	in the Cathode Catalyst Layer and calculates the slope using a numerical 
	estimated derivative dS_dphi. 
:
\*---------------------------------------------------------------------------*/
{
	//  Protonic Transport
	// Request previous iteration to be stored
	potProton.storePrevIter();

	// Calculate the slope
	volScalarField dS_dphi(sourceProtonTotal);
	dS_dphi.dimensions().reset(dimless);
	dS_dphi = 0.;
	dS_dphi.dimensions().reset(sourceProtonTotal.dimensions()/potProton.dimensions());

	// calculate the slope
	{
		// Setup clones of the needed volume fields
		scalar dPhi = slopeStep;
		volScalarField actPotLeft_(actPot);
		volScalarField actPotRight_(actPot);
		volScalarField RRANODICLeft_(sourceProtonTotal);
		volScalarField RRANODICRight_(sourceProtonTotal);
		volScalarField RRCATHODICLeft_(sourceProtonTotal);
		volScalarField RRCATHODICRight_(sourceProtonTotal);
		volScalarField RRLeft_(sourceProtonTotal);
		volScalarField RRRight_(sourceProtonTotal);

		// Reset dims, as required
		RRANODICLeft_.dimensions().reset(dimless);
		RRANODICRight_.dimensions().reset(dimless);
		RRCATHODICLeft_.dimensions().reset(dimless);
		RRCATHODICRight_.dimensions().reset(dimless);

		// Create the change in potProton
		actPotLeft_.internalField() = 
			potElectron.internalField() - potProton.internalField() - 1.23 - dPhi;
		actPotLeft_.boundaryField() = 
			potElectron.boundaryField() - potProton.boundaryField() - 1.23 - dPhi;
		actPotRight_.internalField() = 
			potElectron.internalField() - potProton.internalField() - 1.23 + dPhi;
		actPotRight_.boundaryField() = 
			potElectron.boundaryField() - potProton.boundaryField() - 1.23 + dPhi;

		// Calculate the change in Reaction Rate for potProton +/- SMALL
		RRANODICLeft_.internalField() = exchCurrBV_CCL.value()
			*(Foam::exp(2.*0.5*actPotLeft_.internalField()*96485./(8.314*353.)));
		RRANODICLeft_.boundaryField() =  exchCurrBV_CCL.value()
			*(Foam::exp(2.*0.5*actPotLeft_.boundaryField()*96485./(8.314*353.)));
		RRANODICRight_.internalField() = exchCurrBV_CCL.value()
			*(Foam::exp(2.*0.5*actPotRight_.internalField()*96485./(8.314*353.)));		
		RRANODICRight_.boundaryField() = exchCurrBV_CCL.value()
			*(Foam::exp(2.*0.5*actPotRight_.boundaryField()*96485./(8.314*353.)));

		RRCATHODICLeft_.internalField() = exchCurrBV_CCL.value()
			*(-Foam::exp(-2.*0.5*96485.*actPotLeft_.internalField()/(8.314*353.)));
		RRCATHODICLeft_.boundaryField() = exchCurrBV_CCL.value()
			*(-Foam::exp(-2.*0.5*96485.*actPotLeft_.boundaryField()/(8.314*353.)));
		RRCATHODICRight_.internalField() = exchCurrBV_CCL.value()
			*(-Foam::exp(-2.*0.5*96485.*actPotRight_.internalField()/(8.314*353.)));
		RRCATHODICRight_.boundaryField() = exchCurrBV_CCL.value()
			*(-Foam::exp(-2.*0.5*96485.*actPotRight_.boundaryField()/(8.314*353.)));

		// reset dims on req. variables
		RRANODICLeft_.dimensions().reset(sourceProtonTotal.dimensions());
		RRCATHODICLeft_.dimensions().reset(sourceProtonTotal.dimensions());
		RRANODICRight_.dimensions().reset(sourceProtonTotal.dimensions());
		RRCATHODICRight_.dimensions().reset(sourceProtonTotal.dimensions());

		// Set the value of the RRLeft_ and RRRight_
		RRLeft_ = RRANODICLeft_ + RRCATHODICLeft_;
		RRRight_ = RRANODICRight_ + RRCATHODICRight_;
		dS_dphi = (RRRight_-RRLeft_)/(actPotRight_-actPotLeft_);
	}

	dS_dphi.dimensions().reset(dimless);
	dS_dphi.dimensions().reset(sourceProtonTotal.dimensions()/potProton.dimensions());

	// End calculate the slope

	fvScalarMatrix protonicTransport
	(
	 //   fvm::ddt(potProton)
	    -fvm::laplacian(condProtonEff,potProton) + fvm::SuSp(dS_dphi,potProton)
		== 
		dS_dphi*potProton.prevIter() + sourceProtonTotal
	);

	// Relax the Equation
	protonicTransport.relax();
	// Compute the Flux (may not be needed, TBD)
	protonFlux = protonicTransport.flux();
	// Solve the Equation
	protonicTransport.solve();
	// Relax the Solved Variable
	potProton.relax();
	// Correct the Boundary Conditions
	potProton.correctBoundaryConditions();
}
