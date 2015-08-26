// Start of New file
//
//
// Calculate the Field Variables
	// Stored Gradients
	{
		pMixGrad = fvc::grad(pMix);
		potElectronGrad = fvc::grad(potElectron);
		potProtonGrad = fvc::grad(potProton);
	}
	// Morphology Fields
	{
		ratioPtSupp = ratioPtCInput_ACL;
		ratioCatIon = ratioCatIonInput_ACL;

		massSupp = ((1.-ratioPtSupp)/ratioPtSupp)*massPtInput_ACL;
		massIon = (ratioCatIon/(1.-ratioCatIon))*(massSupp+massPtInput_ACL);

		if(1==anodeInputECSA)
		{			
			sAC = ((experimentEPSAInput_ACL*geomActiveArea)/volumeACL);
		}
		else
		{
			Info<< nl << "anode Particle Size Calculation Selected - this is not available yet"
				<< nl << " DEFAULTING TO EXPERIMENTAL INPUT"
				<< nl << endl;
			sAC = ((experimentEPSAInput_ACL*geomActiveArea)/volumeACL);
		}

		volFPt = (massPtInput_ACL*geomActiveArea/densityPtInput_ACL)/(volumeACL);
		volFSupp = (massSupp*geomActiveArea/densitySupportInput_ACL)/(volumeACL);
		volFIon = (massIon*geomActiveArea/densityIonomerInput_ACL)/(volumeACL);
		volFPore = 1.-volFPt-volFSupp-volFIon;
		volFPoreEff = volFPore*(1-sLiq);
		tIonFilm = (volFIon*(volumeACL))/(sAC*volumeACL);
	
		//sVAgg 
		//numAgg = 
		//volFSuppAgg = 
		//volFPtAgg = 
		//volFIonAgg = 
	}
	// Material Property Fields
	{
		kPermDry = dryPermYInput_ACL;
		kPermRelLiq = kPermDry*pow(sLiq, sCorrOrderPermInput_ACL);
		kPermRelGas = kPermDry*pow(volFPoreEff, sCorrOrderPermInput_ACL); // Need to revisit the order of this
		
		{
			dimensionedScalar henryCoeff_
			(
				"henryCoeff_",
				dimensionSet( 1, 2, -2, 0,-1, 0, 0 ),
				0.255e5
			);

			dimensionedScalar expCoeff_
			(
				"expCoeff_", 
				dimensionSet(T.dimensions()), 
				-170.
			);

			henrySolH2 =henryCoeff_*exp(expCoeff_/T);
		}
		//diffH2Ion =
		//diffH2IonAgg =
	}
	// Nitrogen Transport Fields
	{
		wN2 = 1. - wH2 - wH2OVap;
		Foam::max
		(
			wN2, 
			dimensionedScalar
			(
				"wN2LowerLimit", 
				dimensionSet(wN2.dimensions()), 
				Foam::doubleScalarVSMALL
			)
		);
	}
	// Mole Fraction Fields
	{
		xH2 = (wH2/Fast::Constants::molWeightH2)/((wH2/Fast::Constants::molWeightH2)
			+(wN2/Fast::Constants::molWeightN2)+(wH2OVap/Fast::Constants::molWeightH2O));
		xN2 = (wN2/Fast::Constants::molWeightN2)/((wH2/Fast::Constants::molWeightH2)
			+(wN2/Fast::Constants::molWeightN2)+(wH2OVap/Fast::Constants::molWeightH2O));
		xH2OVap =(wH2OVap/Fast::Constants::molWeightH2O)/((wH2/Fast::Constants::molWeightH2)
				+(wN2/Fast::Constants::molWeightN2)+(wH2OVap/Fast::Constants::molWeightH2O));
	}
	// Fluid Property Fields
	{
		viscH2 = Fast::Constants::viscH2Ref;
		viscN2 = Fast::Constants::viscN2Ref;
		viscH2OVap = Fast::Constants::viscH2OVapRef;
	}
	// Binary Diffusivity Fields
	{
		dimensionedScalar diffCoeff_("diffCoeff_", dimensionSet( 0.5, 4, -1,-1.75,-0.5, 0, 0), 1.e-7);
		
		diffH2N2 = diffCoeff_*Foam::pow(T,1.75)*Foam::pow(((1.e-3/Fast::Constants::molWeightH2)
					+(1.e-3/Fast::Constants::molWeightN2)),0.5)/((pMix/Fast::Constants::pStd)
					*Foam::pow((Foam::pow(Fast::Constants::diffVolH2,(1./3.))
					+Foam::pow(Fast::Constants::diffVolN2,(1./3.))),2.));

		diffH2H2O = diffCoeff_*Foam::pow(T,1.75)*Foam::pow(((1.e-3/Fast::Constants::molWeightH2)
					+(1.e-3/Fast::Constants::molWeightH2O)),0.5)/((pMix/Fast::Constants::pStd)
					*Foam::pow((Foam::pow(Fast::Constants::diffVolH2,(1./3.))
					+Foam::pow(Fast::Constants::diffVolH2OVap,(1./3.))),2.));

		diffN2H2O = diffCoeff_*Foam::pow(T,1.75)*Foam::pow(((1.e-3/Fast::Constants::molWeightN2)
					+(1.e-3/Fast::Constants::molWeightH2O)),0.5)/((pMix/Fast::Constants::pStd)
					*Foam::pow((Foam::pow(Fast::Constants::diffVolN2,(1./3.))
					+Foam::pow(Fast::Constants::diffVolH2OVap,(1./3.))),2.));
	}
	// Mixture Property Fields
	{
		molWeightMix = xH2*Fast::Constants::molWeightH2+xN2*Fast::Constants::molWeightN2
					+xH2OVap*Fast::Constants::molWeightH2O;
		densityMix = pMix*molWeightMix/(Fast::Constants::gasConstantR*T);
		viscMix = xH2*viscH2+xN2*viscN2+xH2OVap*viscH2OVap;
		velAvg = -kPermRelGas*pMixGrad.component(1)/viscMix;
	}
	// Local Condition Fields
	{
		pH2 = xH2*pMix;
		pH2OVap = xH2OVap*pMix;

		dimensionedScalar satPCoeff1_("satPCoeff1_", dimensionSet(pH2OSat.dimensions()), 1.26837e-3);
		dimensionedScalar satPCoeff2_("satPCoeff2_", dimensionSet(pH2OSat.dimensions()), 1.49827);
		dimensionedScalar satPCoeff3_("satPCoeff3_", dimensionSet(pH2OSat.dimensions()), 6.70916e2);
		dimensionedScalar satPCoeff4_("satPCoeff4_", dimensionSet(pH2OSat.dimensions()), 1.34832e5);
		dimensionedScalar satPCoeff5_("satPCoeff5_", dimensionSet(pH2OSat.dimensions()), 1.02503e7);

		volScalarField Ttemp_(T);
		Ttemp_.dimensions().reset(dimless);

		pH2OSat = (satPCoeff1_*Foam::pow(Ttemp_,4.)-satPCoeff2_*Foam::pow(Ttemp_,3.)
				+satPCoeff3_*Foam::pow(Ttemp_,2.)-satPCoeff4_*Ttemp_+satPCoeff5_);

		RH = (pMix*xH2OVap/pH2OSat);
	}
	// Ionomeric Material Fields
	{
		// Water content within the Ionomer
		{
			scalar lowerRHBound_ = 0.04;
			scalar schroedersSwitch_ = 0.999;
			volScalarField RHTemp_ = RH;
			RHTemp_ = Foam::min(RHTemp_, 1.0);

			dimensionedScalar H2ODissIonLowerBound_
			(
			 	"H2ODissIonLowerBound_", 
				dimensionSet(H2ODissIon.dimensions()), 
				0.043+17.81*lowerRHBound_-39.85*Foam::pow(lowerRHBound_,2.)+36.*Foam::pow(lowerRHBound_,3.)
			);
			
			volScalarField switchH2ODissIon_ = pos(schroedersSwitch_-RHTemp_)-1*neg(schroedersSwitch_-RHTemp_);

			H2ODissIon = (1+switchH2ODissIon_)*((0.043+17.81*RHTemp_-39.85*Foam::pow(RHTemp_,2.)
						+36*Foam::pow(RHTemp_,3.))/2.)+(1-switchH2ODissIon_)*((22.)/2.);

			H2ODissIon = Foam::max
						(
						 	H2ODissIon, 
							H2ODissIonLowerBound_
						);
		}
		// Conductivity of Protons within the Ionomer
		{	
			dimensionedScalar condProtonCoeff_
			(
			 	"condProtonCoeff_",
			 	dimensionSet(condProtonCoeff_ACL.dimensions()),
				condProtonCoeff_ACL.value()
			);

			scalar condProtIonCoeff1_ = -0.00326;
			scalar condProtIonCoeff2_ = 0.005139;
			dimensionedScalar arrheniusFactor_("arrheniusFactor_", dimensionSet(0, 0, 0,1, 0, 0, 0), 1268.);
			dimensionedScalar refTemperature_("refTemperature_", dimensionSet(T.dimensions()), 303.);

			condProtonIon = condProtonCoeff_*(condProtIonCoeff1_+condProtIonCoeff2_*H2ODissIon)
							*Foam::exp(arrheniusFactor_*((1./refTemperature_)-(1./T)));
		}
	}		

	// Liquid Water Fields
	{
		jS = (1.417*sLiq-2.120*Foam::pow(sLiq,2)+1.263*Foam::pow(sLiq,3.));
		pCap = Foam::pow((volFPore/kPermDry),0.5)*Fast::Constants::surfTensH2OLiq
			*Foam::cos(contactAngleInput_ACL)*jS;
		{
			volScalarField conditionMin_ = pos(sLiq - 0.);
			volScalarField conditionBelow_ = pos(0.5595 - sLiq);
			
			scalar sLiqNegFlag = Foam::min
			(
			 	Foam::min(conditionMin_.internalField()),
				Foam::min(conditionMin_.boundaryField())
			);
		
			if(0 == sLiqNegFlag)
			{
				Info<< "Liquid water in the ACL has a negative fraction" << nl
					<< "Exiting ...." << nl
					<< endl;
				exit(1);
			}
		
			djSDs = (conditionBelow_)*(1.417-4.24*sLiq+3.789*Foam::pow(sLiq,2.))
				 +(1-conditionBelow_)*(0.2308295065);
		}
		
		dpCapDs = Foam::pow((volFPore/kPermDry),0.5)*Fast::Constants::surfTensH2OLiq
				*Foam::cos(contactAngleInput_ACL)*djSDs;	
		{ 
			volScalarField conditionBelow_ = pos(1.0 - RH);
			volScalarField conditionAbove_ = pos(RH - 1.0);
			swPhaseChangeLiqVap = (1)*(conditionAbove_)+(-1)*(conditionBelow_);
		}
	}
