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
    - currentSourceACL.H

Description
	- Calls the kinetic function to determine the volumetric current source
	and implements a relaxation method for the stiff kinetics
:
\*---------------------------------------------------------------------------*/

// Calculate the volumetric current production
if(1==anodeKineticsTHV)
{
	volScalarField sourceVolCurrentOld_ = sourceVolCurrent;
	sourceVolCurrent = 
 		tafelHeyrovskyVolmer
		(
			sourceVolCurrent,
			pMix,
			xH2,
			T,
			sAC,
			actPot,
			exchCurrTVRefInput_ACL,
			exchCurrHVRefInput_ACL,
			concH2RefInput_ACL,
			reactionOrderH2Input_ACL,
			kappaTHVInput_ACL
		);	
	sourceVolCurrent = sourceVolCurrentOld_ + relaxAnodeKinetics*(sourceVolCurrent-sourceVolCurrentOld_);
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
			xH2,
			henrySolH2,
			nElectronsInput_ACL,
			symmetryFactorBVInput_ACL,
			exchCurrBVRefInput_ACL,
			actEnergyBVInput_ACL,
			tempBVRefInput_ACL,
			concH2RefInput_ACL
		);
}
