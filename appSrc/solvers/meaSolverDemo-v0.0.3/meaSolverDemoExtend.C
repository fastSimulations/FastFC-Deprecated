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
    meaSolverDemoExtend

Description
	Steady State version of the MEA Steady State Performance Model

\*---------------------------------------------------------------------------*/

#include <iostream>
#include <iomanip>
#include <unistd.h>
#include <fstream>
#include <directMappedPatchBase.H>
#include <fvCFD.H>
#include <fastRegionProperties.H>
#include <physicalConstants.H>
#include <tafelHeyrovskyVolmer.H>
#include <butlerVolmerAgglomerate.H>
#include <butlerVolmerAgglomerateAlternate.H>
#include <butlerVolmerDiscrete.H>

int main(int argc, char *argv[])
{

	#include <setRootCase.H>
	#include <licenseStatement.H>

	#include <createTime.H>
	#include <readTimeControls.H>

    fastRegionProperties rp(runTime);

    // Create Component Mesh and Field Objects
    Info<< "Creating component mesh objects at time = " << runTime.timeName() << nl << endl;
	
	// Create the Meshes
	#include <APTLCreateMesh.H>
   	#include <ACLCreateMesh.H>
	#include <PMEMCreateMesh.H>
 	#include <CCLCreateMesh.H>
	#include <CPTLCreateMesh.H>
	
	// Create the fields
	#include <APTLCreateFields.H>
	#include <ACLCreateFields.H>
	#include <PMEMCreateFields.H>
	#include <CCLCreateFields.H>
	#include <CPTLCreateFields.H>

	// Setup the block identifiers
	#include <fastfcDictionary.H>

	Info<<nl<< "Max Iteration set at " << maxSweep << endl;
	
	// Initialization of variables used in looping and convergence tracking - will move these out later
	#include <initSolverLoopingVariables.H>
	
	// Setup OCV -- Should be moved to external calc rather than fixed here
	dimensionedScalar OCV("OCV", dimensionSet( 1, 2, -3, 0, 0, -1, 0 ), 1.1101222924691811);

	// Initialize variables for transport coefficients etc
	#include <initialization.H>

	// Create polarization data file
	{
		std::ofstream polDat("polarizationData.txt");
		polDat << "time" << ", " << "currentDensity" << ", " << "cellVoltage"; 
	}
	
	while ( runTime.run())
	{
		#include <readTimeControls.H>

		runTime++;
		counter = 0;
		currentDensityDiff = 1e15;

		Info<< "Polarization Point #" << runTime.timeName() << nl << endl;

    	while 	(
					(counter<int(maxSweep))
					&&
					(
					 	(currentDensityDiff>currentTolerance)
					 	||
					 	((anodeCurrentDensity/cathodeCurrentDensity)>0.)
					)
				)
		{
			counter++;
	
			// Solve Equations	
//			#include <solveEqnSeperateImplicit.H>
			#include <solveEqnTogetherImplicit.H>

			// Check Current and Charge Convergence
			#include <convergenceOutput.H>
		
	    }

		// Output Solution
//		forAll(CCLRegions, zoneID)
//		{
//			#include <CCLSetFields.H>
//			sourceVolCurrent.write();
//		}
		runTime.write();	
	
		#include <voltageCurrent.H>
	}

 	Info<< "End\n" << endl;

 	return 0;
}


// ************************************************************************* //
