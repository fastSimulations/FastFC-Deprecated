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
    tafelHeyrovskyVolmer.H

Description
	Multi-step anode HOR kinetics model

\*---------------------------------------------------------------------------*/

volScalarField tafelHeyrovskyVolmer
(
 	volScalarField returningField_,
	volScalarField pMix_,
	volScalarField xH2_,
	volScalarField T_,
	volScalarField sAC_,
	volScalarField actPot_,
	dimensionedScalar exchCurrTVRef_,
	dimensionedScalar exchCurrHVRef_,
	dimensionedScalar concH2Ref_,
	dimensionedScalar reactionOrderH2_,
	dimensionedScalar kappaTHV_
)
{
	#include "physicalConstants.H"

	// Set the dimensions and size equivalent to the field to be returned	
	volScalarField volCurrentDensity_(returningField_);

	volCurrentDensity_=
	Foam::pow((pMix_*xH2_)/(Fast::Constants::gasConstantR*T_*concH2Ref_),reactionOrderH2_)*sAC_*
	(
		exchCurrTVRef_*
		(
		 	1-Foam::exp
			(
			 	(-2.*Fast::Constants::faraday*actPot_)/(kappaTHV_*Fast::Constants::gasConstantR*T_)
			)
		)+
		exchCurrHVRef_*
		(
		 	Foam::exp
			(
			 	(Fast::Constants::faraday*actPot_)/(2.*Fast::Constants::gasConstantR*T_)
			)-
			Foam::exp
			(
				(-1.*Fast::Constants::faraday*actPot_)/(kappaTHV_*Fast::Constants::gasConstantR*T_)
			)*
			Foam::exp
			(
			 	(-1.*Fast::Constants::faraday*actPot_)/(2.*Fast::Constants::gasConstantR*T_)
			)	
		)
	);
/*
	volCurrentDensity_.internalField() = 1.e10*(Foam::exp(2.*0.5*actPot_.internalField()*96485./(8.314*353.)))
			+1.e10*(-Foam::exp(-2.*0.5*96485.*actPot_.internalField()/(8.314*353.)));
	volCurrentDensity_.boundaryField() = 1.e10*(Foam::exp(2.*0.5*actPot_.boundaryField()*96485./(8.314*353.)))
			+1.e10*(-Foam::exp(-2.*0.5*96485.*actPot_.boundaryField()/(8.314*353.)));
*/
	return volCurrentDensity_;
}
