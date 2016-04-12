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
    - waterVapourTransportEquationACL.H

Description
	- Solves the water vapour transport equation 
:
\*---------------------------------------------------------------------------*/
{
	// Request previous iteration to be stored
	//	wO2.storePrevIter();
	volScalarField wH2OVapOld_ = wH2OVap;

	// Adapt the transport coefficient 
	transportCoeffH2OVap = densityMix*diffH2OVapMixEff;	

	// Determine Phase change (if two phase is activated)
	#include <ACLPhaseChangeSourceTermFields.H>
	
	// Set the Oxygen Sink/Source Source Term
	#include <ACLWaterVapSourceTermFields.H>
	
	// Calculate the electron transport equation

	fvScalarMatrix waterVapourTransport
	(
	 //   fvm::ddt(wO2)
	    -fvm::laplacian(transportCoeffH2OVap,wH2OVap) //+ fvm::SuSp(dSdPhi_,wH2OVap)
		== 
		sourceMassH2OVap //+ dSdPhi_*wH2OVapOld_
	);
	
	// Relax the Equation
//	waterVapourTransport.relax();

	// Solve the Matrix System
	waterVapourTransport.solve();

	// Relax the solved variable
	scalar relax_ = relaxWaterVapour;
	volScalarField deltaH2OVap_ = wH2OVap - wH2OVapOld_;
	wH2OVap = wH2OVapOld_ + deltaH2OVap_*relax_;
	
	// Minimization Check
	{
		dimensionedScalar wH2OVapMin_
		(
		 	"WH2OVapMin_",
			dimensionSet(wH2OVap.dimensions()),
			Foam::VSMALL
		);

		wH2OVap = Foam::max(wH2OVap,wH2OVapMin_);
	}	

	// Correct the boundary conditions on the solved variable
	wH2OVap.correctBoundaryConditions();
}
