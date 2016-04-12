// Electrochemistry Field Variables
{
	actPot = potElectron - potProton - OCV;

	// Set the dimensions and size equivalent to the field to be returned	
/*	dimensionedScalar unity_
	(
	 	"unity_",
		dimensionSet( 0, 0, 0, 0, 0, 0, 0),
		1.
	);

	chTransCoeffAnodic = nElectronsInput_CCL*(unity_ - symmetryFactorInput_CCL);
	chTransCoeffCathodic = nElectronsInput_CCL*(symmetryFactorInput_CCL);

	exchCurrRefCorr = exchCurrBVRefInput_CCL*
		Foam::exp
		(
		 	(-1.*actEnergyBVInput_CCL/Fast::Constants::gasConstantR)*
			((1./T)-(1./tempBVRefInput_CCL))
		);

	kCAgg =	
		(
			sAC/
			(
			 	4.*Fast::Constants::faraday
				*(
				 	1-volFPore
				)
			)
		 )*
		(
		 	exchCurrRefCorr/concO2BVRefInput_CCL
		)*
		(
		 	Foam::exp
			(
			 	(chTransCoeffAnodic*Fast::Constants::faraday*actPot)/
				(Fast::Constants::gasConstantR*T)
			)-
			Foam::exp
			(
			 	(-chTransCoeffCathodic*Fast::Constants::faraday*actPot)/
				(Fast::Constants::gasConstantR*T)
			)
		);

	theileModulus = (radiusAggInput_CCL/3.)*Foam::pow((Foam::mag(kCAgg)/diffO2IonAgg),0.5); 

	theileModulus = Foam::max(theileModulus,1.e-5);
	
	effectiveFactor =
		(
		 	(
			 	(1./theileModulus)
			)*
		 	(
			 	(1./Foam::tanh(3.*theileModulus))
				-(1./(3.*theileModulus))
			)
		);	
		*/
}
