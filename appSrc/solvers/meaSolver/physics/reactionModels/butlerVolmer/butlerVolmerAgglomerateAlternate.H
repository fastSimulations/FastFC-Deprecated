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
    butlerVolmerAggloemrate.H

Description
	Agglomerate-based Model using Butler-Volmer Kinetcs

\*---------------------------------------------------------------------------*/

volScalarField butlerVolmerAgglomerateAlternate
(
 	volScalarField returningField_,
	volScalarField xO2_,
	volScalarField pMix_,
	volScalarField kCAgg_,
	volScalarField effectiveFactor_,
	volScalarField volFPore_,
	volScalarField tIonFilm_,
	volScalarField sVAgg_,
	volScalarField diffO2Ion_,
	volScalarField henrySolO2_,
	dimensionedScalar radiusAgg_ 
)
{
	#include "physicalConstants.H"

	// Set the dimensions and size equivalent to the field to be returned	
	volScalarField volCurrentDensity_(returningField_);

			
	volCurrentDensity_= (Foam::mag(kCAgg_)/kCAgg_)*(4.*Fast::Constants::faraday*(pMix_*xO2_/henrySolO2_))/
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
