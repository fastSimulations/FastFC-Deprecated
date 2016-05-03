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
    - PMEMGeneralFields.H

Description
	- Calculates the value of the field variables within the polymer membrane
	layer
:
\*---------------------------------------------------------------------------*/

// Calculate the Field Variables
	// Stored Gradients
	{
		potProtonGrad = fvc::grad(potProton);
	}
	// Ionomeric Material Fields
	{

		// Water content within the Ionomer
		{
			scalar lowerRHBound_ = 0.04;
			scalar schroedersSwitch_ = 0.999;
			volScalarField RHTemp_(H2ODissIon*0.);
			RHTemp_ = pmemRH;
			RHTemp_ = Foam::min(RHTemp_,1.0);
			RHTemp_ = Foam::max(RHTemp_, lowerRHBound_);

			H2ODissIon = ((0.043+17.81*RHTemp_-39.85*Foam::pow(RHTemp_,2.)+36*Foam::pow(RHTemp_,3.)));
//Temporarily disabled schroeders paradox

//			dimensionedScalar H2ODissIonLowerBound_
//			(
//			 	"H2ODissIonLowerBound_", 
//				dimensionSet(H2ODissIon.dimensions()), 
//				0.043+17.81*lowerRHBound_-39.85*Foam::pow(lowerRHBound_,2.)+36.*Foam::pow(lowerRHBound_,3.)
//			);
			
//			volScalarField switchH2ODissIon_ = pos(schroedersSwitch_-RHTemp_)-1*neg(schroedersSwitch_-RHTemp_);

//			H2ODissIon = (1+switchH2ODissIon_)
//						*((0.043+17.81*RHTemp_-39.85*Foam::pow(RHTemp_,2.)+36*Foam::pow(RHTemp_,3.))/2.)
//						+(1-switchH2ODissIon_)*
//						((14.)/2.); //Temporarily disabled schroeders paradox

//			H2ODissIon = Foam::max
//						(
//						 	H2ODissIon, 
//							H2ODissIonLowerBound_
//						);
		}

		// Conductivity of Protons within the Ionomer
		{
			dimensionedScalar condProtonCoeff_
			(
			 	"condProtonCoeff_",
			 	dimensionSet(condProtonCoeff_PMEM.dimensions()),
				condProtonCoeff_PMEM.value()
			);

			scalar condProtIonCoeff1_ = -0.00326;
			scalar condProtIonCoeff2_ = 0.005139;
			dimensionedScalar arrheniusFactor_("arrheniusFactor_", dimensionSet(0, 0, 0,1, 0, 0, 0), 1268.);
			dimensionedScalar refTemperature_("refTemperature_", dimensionSet(T.dimensions()), 303.);

			condProtonIon = condProtonCoeff_*(condProtIonCoeff1_+condProtIonCoeff2_*H2ODissIon)
							*Foam::exp(arrheniusFactor_*((1./refTemperature_)-(1./T)));
		}

		// Thermal conductivity of the Ionomer
		{
			condThermalIon = condThermalIonInput_PMEM;
		}
	}		
