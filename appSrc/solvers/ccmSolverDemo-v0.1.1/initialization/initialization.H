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
    initialization.H

Description
	Initialisation of the solution set and calculated variables at 
	the beginning of the simulation.

\*---------------------------------------------------------------------------*/
{
	// Setup the Calculated Fields, i.e. conductivities
	forAll(ACLRegions, zoneID)
	{
		#include "ACLSetFields.H"
		{
			condElectronEff = condElectron_ACL;
			condProtonEff = condProton_ACL;
		}
	}
	forAll(PMEMRegions, zoneID)
	{
		#include "PMEMSetFields.H"
		{
			condProtonEff = condProton_PMEM;
		}
	}
	forAll(CCLRegions, zoneID)
	{
		#include "CCLSetFields.H"
		{
			condElectronEff = condElectron_CCL;
			condProtonEff = condProton_CCL;
		}
	}

	// Write initial guess solution fields
	Info<< nl
		<< "Initialisation of Solution Fields"
		<< endl;

	forAll(ACLRegions, zoneID)
	{
		#include "ACLSetFields.H"
		potElectron.write();
		potProton.write();
		actPot.write();
	}

	forAll(PMEMRegions, zoneID)
	{
		#include "PMEMSetFields.H"
		potProton.write();
	}

	forAll(CCLRegions, zoneID)
	{
		#include "CCLSetFields.H"
		potElectron.write();
		potProton.write();
		actPot.write();
	}
}
