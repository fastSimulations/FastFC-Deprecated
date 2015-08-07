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
    - electronTransportEquationA.H

Description
	- Solves the electron transport equation with explicit/implicit terms 
	in the left block and calculates the slope using a numerical 
	estimated derivative dS_dphi. 
:
\*---------------------------------------------------------------------------*/
{
	//  Protonic Transport
	// Request previous iteration to be stored
	potElectron.storePrevIter();

	// Calculate the slope
	volScalarField dS_dphi(sourceElectronTotal);
	dS_dphi.dimensions().reset(dimless);
	dS_dphi = 0.;
	dS_dphi.dimensions().reset(sourceElectronTotal.dimensions()/potElectron.dimensions());

	// calculate the slope
	{
		// Setup clones of the needed volume fields
		scalar dPhi = slopeStep;
		volScalarField actPotLeft_(actPot);
		volScalarField actPotRight_(actPot);
		volScalarField RRANODICLeft_(sourceElectronTotal);
		volScalarField RRANODICRight_(sourceElectronTotal);
		volScalarField RRCATHODICLeft_(sourceElectronTotal);
		volScalarField RRCATHODICRight_(sourceElectronTotal);
		volScalarField RRLeft_(sourceElectronTotal);
		volScalarField RRRight_(sourceElectronTotal);

		// Reset dims as required
		RRANODICLeft_.dimensions().reset(dimless);
		RRANODICRight_.dimensions().reset(dimless);
		RRCATHODICLeft_.dimensions().reset(dimless);
		RRCATHODICRight_.dimensions().reset(dimless);

		// Create the change in potElectron
		actPotLeft_.internalField() = 
			potElectron.internalField() - potProton.internalField() -dPhi;
		actPotLeft_.boundaryField() = 
			potElectron.boundaryField() - potProton.boundaryField() -dPhi;
		actPotRight_.internalField() = 
			potElectron.internalField() - potProton.internalField() +dPhi;
		actPotRight_.boundaryField() = 
			potElectron.boundaryField() - potProton.boundaryField() +dPhi;

		// Calculate the change in Reaction Rate for potElectron +/- SMALL
		RRANODICLeft_.internalField() = exchCurrBV_ACL.value()
			*(Foam::exp(2.*0.5*actPotLeft_.internalField()*96485./(8.314*353.)));
		RRANODICLeft_.boundaryField() =  exchCurrBV_ACL.value()
			*(Foam::exp(2.*0.5*actPotLeft_.boundaryField()*96485./(8.314*353.)));
		RRANODICRight_.internalField() = exchCurrBV_ACL.value()
			*(Foam::exp(2.*0.5*actPotRight_.internalField()*96485./(8.314*353.)));				
		RRANODICRight_.boundaryField() = exchCurrBV_ACL.value()
			*(Foam::exp(2.*0.5*actPotRight_.boundaryField()*96485./(8.314*353.)));

		RRCATHODICLeft_.internalField() = exchCurrBV_ACL.value()
			*(-Foam::exp(-2.*0.5*96485.*actPotLeft_.internalField()/(8.314*353.)));
		RRCATHODICLeft_.boundaryField() = exchCurrBV_ACL.value()
			*(-Foam::exp(-2.*0.5*96485.*actPotLeft_.boundaryField()/(8.314*353.)));
		RRCATHODICRight_.internalField() = exchCurrBV_ACL.value()
			*(-Foam::exp(-2.*0.5*96485.*actPotRight_.internalField()/(8.314*353.)));
		RRCATHODICRight_.boundaryField() = exchCurrBV_ACL.value()
			*(-Foam::exp(-2.*0.5*96485.*actPotRight_.boundaryField()/(8.314*353.)));

		// reset dims on req. variables
		RRANODICLeft_.dimensions().reset(sourceElectronTotal.dimensions());
		RRCATHODICLeft_.dimensions().reset(sourceElectronTotal.dimensions());
		RRANODICRight_.dimensions().reset(sourceElectronTotal.dimensions());
		RRCATHODICRight_.dimensions().reset(sourceElectronTotal.dimensions());

		// Set the value of the RRLeft_ and RRRight_
		RRLeft_ = RRANODICLeft_ + RRCATHODICLeft_;
		RRRight_ = RRANODICRight_ + RRCATHODICRight_;

		dS_dphi = (RRRight_-RRLeft_)/(actPotRight_-actPotLeft_);
	}
	
	dS_dphi.dimensions().reset(dimless);
	dS_dphi.dimensions().reset(sourceElectronTotal.dimensions()/potElectron.dimensions());

	// End calculate the slope

	fvScalarMatrix electronicTransport
	(
	 //   fvm::ddt(potElectron)
	    -fvm::laplacian(condElectronEff,potElectron) + fvm::SuSp(dS_dphi,potElectron)
		== 
		dS_dphi*potElectron.prevIter() - sourceElectronTotal 
	);

	// Relax the Equation
	electronicTransport.relax();
	// Store the flux (may not be needed, TBD)
	electronFlux = electronicTransport.flux();
	// Solve the Matrix System
	electronicTransport.solve();
	// Relax the solved variable
	potElectron.relax();
	// Correct the boundary conditions on the solved variable
	potElectron.correctBoundaryConditions();
}
