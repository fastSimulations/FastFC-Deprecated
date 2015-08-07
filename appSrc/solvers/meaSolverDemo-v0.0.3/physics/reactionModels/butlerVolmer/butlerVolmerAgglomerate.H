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
    butlerVolmerAggloemrate.H

Description
	Agglomerate-based Model using Butler-Volmer Kinetcs

\*---------------------------------------------------------------------------*/

volScalarField butlerVolmerAgglomerate
(
 	volScalarField returningField_,
	volScalarField T_,
	volScalarField sAC_,
	volScalarField volFPore_,
	volScalarField actPot_,
	volScalarField diffO2Ion_,
	volScalarField pMix_,
	volScalarField xO2_,
	volScalarField henrySolO2_,
	volScalarField tIonFilm_,
	volScalarField sVAgg_,
	volScalarField diffO2IonAgg_,
	dimensionedScalar nElectrons_,
	dimensionedScalar symmetryFactor_,
	dimensionedScalar exchCurrBVRef_,
	dimensionedScalar actEnergyBV_,
	dimensionedScalar tempBVRef_,
	dimensionedScalar concO2BVRef_,
	dimensionedScalar radiusAgg_ 
)
{
	#include "physicalConstants.H"

	// Set the dimensions and size equivalent to the field to be returned	
	volScalarField volCurrentDensity_(returningField_);
	
	dimensionedScalar unity_
	(
	 	"unity_",
		dimensionSet( 0, 0, 0, 0, 0, 0, 0),
		1.
	);

	dimensionedScalar chTransCoeffAnodic_= (nElectrons_*(unity_ - symmetryFactor_));
	dimensionedScalar chTransCoeffCathodic_ = nElectrons_*(symmetryFactor_);

	volScalarField exchCurrRefCorr_ = exchCurrBVRef_
		*Foam::exp
		(
		 	(-1.*actEnergyBV_/Fast::Constants::gasConstantR)
			*(
			 	(1./T_)-(1./tempBVRef_)
			)
		);


	volScalarField kCAgg_ = 
			(
			 	sAC_/
				(
				 	4.*Fast::Constants::faraday
					*(
					 	1-volFPore_
					)
				)
			)
			*(
			 	exchCurrRefCorr_/concO2BVRef_
			)
			*(
			 	Foam::exp
				(
				 	(chTransCoeffAnodic_*Fast::Constants::faraday*actPot_)/
					(Fast::Constants::gasConstantR*T_)
				)-
				Foam::exp
				(
				 	(-chTransCoeffCathodic_*Fast::Constants::faraday*actPot_)/
					(Fast::Constants::gasConstantR*T_)
				)
			);
			
	dimensionedScalar theileModulusCoeff_
	(
		"theileModulusCoeff_",
   		dimensionSet(kCAgg_.dimensions()),
		Foam::VSMALL
	);

	volScalarField theileModulus_ = (radiusAgg_/3.)*Foam::pow((Foam::mag(kCAgg_)/diffO2IonAgg_),0.5);

	theileModulus_ = Foam::max(theileModulus_,1.e-5);


	volScalarField effectiveFactor_ = 
			(
			 	(
				 	(1./theileModulus_)
				)*
			 	(
				 	(1./Foam::tanh(3.*theileModulus_))
					-(1./(3.*theileModulus_))
				)
			);

	dimensionedScalar smallKCAGG_
	(
	 	"smallKCAGG_",
		dimensionSet(kCAgg_.dimensions()),
		Foam::VSMALL
	);

	volCurrentDensity_= -1.*(4.*Fast::Constants::faraday*(pMix_*xO2_/henrySolO2_))/
			(
			 	(
				 	1./
					(
					 	effectiveFactor_*Foam::mag((kCAgg_))*(1.-volFPore_)
					)
				)+
				(
				 	(
					 	(radiusAgg_+tIonFilm_)*tIonFilm_
					)/
					(
					 	sVAgg_*radiusAgg_*diffO2Ion_
					)
				)
			);

	return volCurrentDensity_;
}
