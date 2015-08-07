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
    ccmSolverDemoExtend
Version
	v0.1.1
Description
	Steady State training version of the CCM Steady State 2-3 Eqn Model

\*---------------------------------------------------------------------------*/


#include <iostream>
#include <iomanip>
#include <unistd.h>
#include <fstream>
#include "fvCFD.H"
#include "threeBlockRegionProperties.H"

int main(int argc, char *argv[])
{
	#include "setRootCase.H"
	#include "createTime.H"
	#include "readTimeControls.H"

    threeBlockRegionProperties rp(runTime);

    // Create Component Mesh and Field Objects
    Info<< "Creating component mesh objects at time = " << runTime.timeName() << nl << endl;
	
	// Create the Meshes
   	#include "ACLCreateMesh.H"
	#include "PMEMCreateMesh.H"
 	#include "CCLCreateMesh.H"
	
	// Create the fields
	#include "ACLCreateFields.H"
	#include "PMEMCreateFields.H"
	#include "CCLCreateFields.H"

	// Setup the block identifiers
	#include "fastfcDictionary.H"
	
	Info<<nl<< "Max Iteration set at " << maxSweep << endl;

	
	// Initialization of variables used in looping and convergence tracking - will move these out later
	scalar counter = 0;
	scalar writeCounter = 0;
	scalar anodeCurrentDensity = 0.;
	scalar cathodeCurrentDensity = 0.;
	scalar currentDensityDiff = 1.e15;
	scalar stepTimeOld = 0.;
	scalar totalTime = 0.;

	// Setup OCV -- Should be moved to external calc rather than fixed here
	dimensionedScalar OCV("OCV", dimensionSet( 1, 2, -3, 0, 0, -1, 0 ), 1.18);

	// Determine active area
	#include "activeArea.H"
	// Initialize variables for transport coefficients etc
	#include "initialization.H"

	// Create polarization data file
	{
		std::ofstream polDat("polarizationData.txt");
		polDat << "time" << ", " << "cellVoltage" << ", " << "currentDensity"; 
	}

	while(runTime.run())
	{
		#include "readTimeControls.H"
		runTime++;
		counter = 0;
		currentDensityDiff = 1e15;

		Info<< "Polarization Point #" << runTime.timeName() << nl << endl;

    	while
		(
			(counter<int(maxSweep))
			&&
			(
			 	(currentDensityDiff>currentTolerance)
				||
				(anodeCurrentDensity/cathodeCurrentDensity>0.)
			)
		)
    	{
			counter++;

			// Calculated the changing properties and fields
			forAll(ACLRegions, zoneID)
			{
		 		#include "ACLSetFields.H"
				#include "ACLCalcFields.H"
			}
			forAll(CCLRegions, zoneID)
			{
				#include "CCLSetFields.H"
				#include "CCLCalcFields.H"
			}
			forAll(PMEMRegions, zoneID)
			{
		   		#include "PMEMSetFields.H"
				#include "PMEMCalcFields.H"
			}
		
			// Solve Equations	
			#include "solveEqnSeparateImplicit.H"

			// Check Current and Charge Convergence
			#include "convergenceOutput.H"

			// Output Solution
			//runTime.write();
    	}

		forAll(ACLRegions, zoneID)
		{
			#include "ACLSetFields.H"
			potElectron.write();
			potProton.write();
			actPot.write();
			sourceElectronTotal.write();
			reactionRateBVSum.write();
			condProtonEff.write();
			condElectronEff.write();
		}

		forAll(PMEMRegions, zoneID)
		{
			#include "PMEMSetFields.H"
			potProton.write();
			condProtonEff.write();
		}

		forAll(CCLRegions, zoneID)
		{
			#include "CCLSetFields.H"
			potElectron.write();
			potProton.write();
			actPot.write();
			sourceElectronTotal.write();
			reactionRateBVSum.write();
			condProtonEff.write();
			condElectronEff.write();
		}

	//	#include "voltageCurrent.H"

	}
		
 	Info<< "End\n" << endl;

 	return 0;
}


// ************************************************************************* //
