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
    - CPTLCoefficientFields.H

Description
	- Calculate the transport coefficient terms for the Cathode Porous 
	Trasport layer. 
:
\*---------------------------------------------------------------------------*/
{
	// Electron Transport Fields
	condElectronEff = condElectronCarbFibInput_CPTL;
				// = condElectronCarbFibInput_CPTL*Foam::pow((1.-volFPore), torCarbFibXInput_CPTL);
	
	// Oxygen Transport Fields
	diffO2MixEff = ((1.-xO2)/((xN2/diffO2N2)+(xH2OVap/diffO2H2O)))
					*Foam::pow(volFPoreEff,torPoreYInput_CPTL);
		
	//dimensionedScalar diffO2MixEffTemp("diffO2MixEffTemp", dimensionSet(diffO2MixEff.dimensions()), 1e-6);
	//diffO2MixEff = diffO2MixEffTemp;
	
	// Water Vapour Transport Fields
	
	diffH2OVapMixEff = (1.-xH2OVap)/((xN2/diffN2H2O)+(xO2/diffO2H2O))
						*Foam::pow(volFPoreEff,torPoreYInput_CPTL);
	
	// Energy Transport Fields
	condThermalEff = condThermalCarbFibInput_CPTL;
	// = condThermalCarbFibInput_CPTL*Foam::pow((1-volFPore),torCarbFibYInput_CPTL);
	
	// Volume Averaged thermal conductivity - test after implementation of main code
	//((volFPt*condThermalPtInput_CPTL)+(volFSupp*condThermalSuppInput_CPTL)+(volFIon*condThermalIonInput_CPTL))*Foam::pow((volFPt+volFSupp+volFIon),(torSuppYInput_CPTL*(condThermalSuppInput_CPTL/(condThermalSuppInput_CPTL+condThermalIonInput_CPTL)))+(torIonYInput_CPTL*(condThermalIonInput_CPTL/(condThermalSuppInput_CPTL+condThermalIonInput_CPTL))));
  		
	// Mass Transport (Pressure) Fields
	dynamicPerm = densityMix*(kPermRelGas/viscMix);
	
 	// Liquid Water Transport Fields
	diffH2OLiq = (Fast::Constants::densityH2OLiq*kPermRelLiq/Fast::Constants::viscH2OLiqRef);
}
