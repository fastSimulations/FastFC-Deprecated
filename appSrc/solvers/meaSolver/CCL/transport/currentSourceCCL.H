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
    - currentSourceCCL.H

Description
	- Calls the kinetic function to determine the volumetric current source
	and implements a relaxation method for the stiff kinetics
:
\*---------------------------------------------------------------------------*/

// Calculate the volumetric current production
if(1==cathodeKineticsAgglomerate)
{
	volScalarField sourceVolCurrentOld_ = sourceVolCurrent;
 
 	sourceVolCurrent = 
 		butlerVolmerAgglomerate
		(
			sourceVolCurrent,
			T,
			sAC,
			volFPore,
			actPot,
			diffO2Ion,
			pMix,
			xO2,
			henrySolO2,
			tIonFilm,
			sVAgg,
			diffO2IonAgg,
			nElectronsInput_CCL,
			symmetryFactorInput_CCL,
			exchCurrBVRefInput_CCL,
			actEnergyBVInput_CCL,
			tempBVRefInput_CCL,
			concO2BVRefInput_CCL,
			radiusAggInput_CCL 
		);

	sourceVolCurrent = sourceVolCurrentOld_ + relaxCathodeKinetics*(sourceVolCurrent-sourceVolCurrentOld_);
}
else
{
	sourceVolCurrent = 
		butlerVolmerDiscrete
		(
 			sourceVolCurrent,
			T,
			sAC,
			actPot,
			pMix,
			xO2,
			henrySolO2,
			nElectronsInput_CCL,
			symmetryFactorInput_CCL,
			exchCurrBVRefInput_CCL,
			actEnergyBVInput_CCL,
			tempBVRefInput_CCL,
			concO2BVRefInput_CCL
		);
}
