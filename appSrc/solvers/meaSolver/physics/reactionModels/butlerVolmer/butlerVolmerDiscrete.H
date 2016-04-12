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
    butlerVolmerDiscrete.H

Description
	Discrete-based Model using Butler-Volmer Kinetcs

\*---------------------------------------------------------------------------*/

volScalarField butlerVolmerDiscrete
(
 	volScalarField returningField_,
	volScalarField T_,
	volScalarField sAC_,
	volScalarField actPot_,
	volScalarField pMix_,
	volScalarField xSpecies_,
	volScalarField henrySol_,
	dimensionedScalar nElectrons_,
	dimensionedScalar symmetryFactor_,
	dimensionedScalar exchCurrBVRef_,
	dimensionedScalar actEnergy_,
	dimensionedScalar tempBVRef_,
	dimensionedScalar concBVRef_
)
{
	#include "physicalConstants.H"

	volScalarField volCurrentDensity_(returningField_);

	// Set the dimensions and size equivalent to the field to be returned	
	dimensionedScalar unity_
	(
	 	"unity_",
		dimensionSet( 0, 0, 0, 0, 0, 0, 0),
		1.
	);

	dimensionedScalar chTransCoeffAnodic_
	(
 		"chTransCoeffAnodic_",
   		dimensionSet(nElectrons_.dimensions()),
		nElectrons_.value()*(unity_.value() - symmetryFactor_.value())
	);


	dimensionedScalar chTransCoeffCathodic_
	(
 		"chTransCoeffCathodic_",
   		dimensionSet(nElectrons_.dimensions()),
		nElectrons_.value()*(symmetryFactor_.value())
	);

	volScalarField exchCurrRefCorr_ = exchCurrBVRef_
		*Foam::exp
		(
		 	(-actEnergy_/Fast::Constants::gasConstantR)
			*(
			 	(1./T_)-(1/tempBVRef_)
			)
		);

	volCurrentDensity_.internalField() = 
			exchCurrRefCorr_.internalField()*sAC_.internalField()*(((pMix_.internalField()*xSpecies_.internalField())/henrySol_.internalField())/concBVRef_.value())*(Foam::exp(2.*chTransCoeffAnodic_.value()*actPot_.internalField()*Fast::Constants::faraday.value()/(Fast::Constants::gasConstantR.value()*T_.internalField())))
			+exchCurrRefCorr_.internalField()*sAC_.internalField()*(((pMix_.internalField()*xSpecies_.internalField())/henrySol_.internalField())/concBVRef_.value())*(-Foam::exp(-2.*chTransCoeffCathodic_.value()*Fast::Constants::faraday.value()*actPot_.internalField()/(Fast::Constants::gasConstantR.value()*T_.internalField())));

	volCurrentDensity_.boundaryField() = 
			exchCurrRefCorr_.boundaryField()*sAC_.boundaryField()*(((pMix_.boundaryField()*xSpecies_.boundaryField())/henrySol_.boundaryField())/concBVRef_.value())*(Foam::exp(2.*chTransCoeffAnodic_.value()*actPot_.boundaryField()*Fast::Constants::faraday.value()/(Fast::Constants::gasConstantR.value()*T_.boundaryField())))
			+exchCurrRefCorr_.boundaryField()*sAC_.boundaryField()*(((pMix_.boundaryField()*xSpecies_.boundaryField())/henrySol_.boundaryField())/concBVRef_.value())*(-Foam::exp(-2.*chTransCoeffCathodic_.value()*Fast::Constants::faraday.value()*actPot_.boundaryField()/(Fast::Constants::gasConstantR.value()*T_.boundaryField())));
			 	
	return volCurrentDensity_;
}
