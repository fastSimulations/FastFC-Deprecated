// Start of New file
//
//
// Calculate the Field Variables
	// Stored Gradients
	{
		potProtonGrad = fvc::grad(potProton);
	}
	// Ionomeric Material Fields
	{

		// Water content within the Ionomer
		{
			scalar lowerRHBound_ = 0.04;
			scalar schroedersSwitch_ = 0.999;
			volScalarField RHTemp_(H2ODissIon*0.);
			RHTemp_ = pmemRH;
			RHTemp_ = Foam::min(RHTemp_,1.0);


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
			 	dimensionSet(condProtonCoeff_PMEM.dimensions()),
				condProtonCoeff_PMEM.value()
			);

			scalar condProtIonCoeff1_ = -0.00326;
			scalar condProtIonCoeff2_ = 0.005139;
			dimensionedScalar arrheniusFactor_("arrheniusFactor_", dimensionSet(0, 0, 0,1, 0, 0, 0), 1268.);
			dimensionedScalar refTemperature_("refTemperature_", dimensionSet(T.dimensions()), 303.);

			condProtonIon = condProtonCoeff_*(condProtIonCoeff1_+condProtIonCoeff2_*H2ODissIon)
							*Foam::exp(arrheniusFactor_*((1./refTemperature_)-(1./T)));
		}

		// Thermal conductivity of the Ionomer
		{
			condThermalIon = condThermalIonInput_PMEM;
		}
	}		
