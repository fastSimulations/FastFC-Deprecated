Info<< nl << nl << nl << nl << nl << nl << nl << nl << endl;
forAll(ACLRegions, zoneID)
{
	#include "ACLSetFields.H"
	const scalar anodeCurrent = gSum(sourceProtonTotal.internalField()*mesh.V());
	anodeCurrentDensity = (anodeCurrent)/(geomActiveArea.value()*10.);
	scalar reactionRateCurrentDensity = 0.;
	scalar reactionRateMax = 0.;
	scalar reactionRateMin = 0.;
	scalar reactionRateAve = 0.;

	reactionRateCurrentDensity = gSum(sourceVolCurrent.internalField()
								*mesh.V())/(geomActiveArea.value()*10.);
	reactionRateMax = gMax(sourceVolCurrent);
	reactionRateMin = gMin(sourceVolCurrent);
	reactionRateAve = gAverage(sourceVolCurrent);

	
	cout<< "\n"
		<< "***************************** Anode Electrode *****************************"
		<< "\n"
		<< std::left
		<< std::setw(30) << "Total Curent Density" 
		<< "\n"	
		<< std::setw(30) << "Source Term" 
		<< std::setw(30) << "Reaction Rate HOR"
		<< "\n"
		<< std::left << std::setprecision(2) << std::fixed
		<< std::setw(30) << anodeCurrentDensity
		<< std::setw(30) << reactionRateCurrentDensity
		<< "\n"	
		<< resetiosflags(ios_base::floatfield)
		<< std::setprecision(6)
		<< std::setw(30) << "Max Reaction Rate HOR"
		<< std::setw(30) << "Min Reaction Rate HOR"
		<< std::setw(30) << "Ave Reaction Rate HOR"
		<< "\n"
		<< std::setw(30) << reactionRateMax
		<< std::setw(30) << reactionRateMin
		<< std::setw(30) << reactionRateAve
		<< "\n"
		<< std::setw(30) << "Average eta"
		<< std::setw(30) << "Max eta"
		<< std::setw(30) << "Min eta"
		<< "\n"
		<< std::setw(30) << gAverage(actPot)
		<< std::setw(30) << gMax(actPot)
		<< std::setw(30) << gMin(actPot)
		<< "\n"
		<< std::setw(30) << "Average potElectron"
		<< std::setw(30) << "Max potElectron"
		<< std::setw(30) << "Min potElectron"
		<< "\n"
		<< std::setw(30) << gAverage(potElectron)
		<< std::setw(30) << gMax(potElectron)
		<< std::setw(30) << gMin(potElectron)
		<< "\n"
		<< std::setw(30) << "Average potProton"
		<< std::setw(30) << "Max potProton"
		<< std::setw(30) << "Min potProton"
		<< "\n"
		<< std::setw(30) << gAverage(potProton)
		<< std::setw(30) << gMax(potProton)
		<< std::setw(30) << gMin(potProton)
		<< "\n";

}
forAll(CCLRegions, zoneID)
{
	
	#include "CCLSetFields.H"
	const scalar cathodeCurrent = gSum(sourceProtonTotal.internalField()*mesh.V());
	cathodeCurrentDensity = (cathodeCurrent)/(geomActiveArea.value()*10.);

	scalar reactionRateCurrentDensity = gSum(sourceVolCurrent.internalField()
										*mesh.V())/(geomActiveArea.value()*10.);
	scalar reactionRateMax = gMax(sourceVolCurrent);
	scalar reactionRateMin = gMin(sourceVolCurrent);
	scalar reactionRateAve = gAverage(sourceVolCurrent);

	cout<< "\n"
		<< "**************************** Cathode Electrode *****************************"
		<< "\n"
		<< std::left
		<< std::setw(30) << "Total Curent Density" 
		<< "\n"	
		<< std::setw(30) << "Source Term" 
		<< std::setw(30) << "Reaction Rate OOR"
		<< "\n"
		<< std::left << std::setprecision(2) << std::fixed
		<< std::setw(30) << cathodeCurrentDensity
		<< std::setw(30) << reactionRateCurrentDensity
		<< "\n"	
		<< resetiosflags(ios_base::floatfield)
		<< std::setprecision(6)
		<< std::setw(30) << "Max Reaction Rate ORR"
		<< std::setw(30) << "Min Reaction Rate ORR"
		<< std::setw(30) << "Ave Reaction Rate ORR"
		<< "\n"
		<< std::setw(30) << reactionRateMax
		<< std::setw(30) << reactionRateMin
		<< std::setw(30) << reactionRateAve
		<< "\n"
		<< std::setw(30) << "Average eta"
		<< std::setw(30) << "Max eta"
		<< std::setw(30) << "Min eta"
		<< "\n"
		<< std::setw(30) << gAverage(actPot)
		<< std::setw(30) << gMax(actPot)
		<< std::setw(30) << gMin(actPot)
		<< "\n"
		<< std::setw(30) << "Average potElectron"
		<< std::setw(30) << "Max potElectron"
		<< std::setw(30) << "Min potElectron"
		<< "\n"
		<< std::setw(30) << gAverage(potElectron)
		<< std::setw(30) << gMax(potElectron)
		<< std::setw(30) << gMin(potElectron)
		<< "\n"
		<< std::setw(30) << "Average potProton"
		<< std::setw(30) << "Max potProton"
		<< std::setw(30) << "Min potProton"
		<< "\n"
		<< std::setw(30) << gAverage(potProton)
		<< std::setw(30) << gMax(potProton)
		<< std::setw(30) << gMin(potProton)
		<< "\n";

}	
// Check the current imbalance
{
	if(counter>10)
	{
		currentDensityDiff = 
				fabs(fabs(anodeCurrentDensity) - fabs(cathodeCurrentDensity));	
	}

	cout<< "\n" 
		<< "**************************** Iteration Summary ****************************"
		<< "\n"
		<< std::left
		<< std::setw(30) << "anodeCurrentDensity" 
		<< std::setw(30) << "cathodeCurrentDensity"
		<< "\n" 
		<< std::setw(30) << anodeCurrentDensity
		<< std::setw(30) << cathodeCurrentDensity
		<< "\n";
	const scalar stepTime = runTime.elapsedCpuTime();
	const scalar deltaTime = stepTime - stepTimeOld;
	stepTimeOld = stepTime;
	totalTime += deltaTime;
	
	cout<< "\n"
   		<< std::left
		<< std::setw(12) << " "
		<< std::setw(16) << "Convergence Checking"
		<< "\n"
		<< std::setw(12) << "Iteration" 
		<< std::setw(16) << "Current Consv."
		<< std::setw(16) << "Iteration Time" 
		<< std::setw(16) << "Total Time"
		<< "\n"
		<< std::setw(12) << counter
		<< std::setw(16) << std::setprecision(4)
		<< std::setw(16) << currentDensityDiff
		<< std::setw(16) << deltaTime
		<< std::setw(16) << totalTime
		<< "\n";
}
cout<< resetiosflags(ios_base::floatfield);
