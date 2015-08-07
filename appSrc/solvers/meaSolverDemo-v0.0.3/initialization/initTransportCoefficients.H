/*---------------------------------------------------------------------------*\
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

Application
   meaSolverDemo 

File
	initTransportCoefficients.H

Description
	Setup of the transport ceofficients using a "dummy" coefficient

Notes

Developers
	David B. Harvey
\*---------------------------------------------------------------------------*/ 

// Setup a dummy set of transport coefficients for the domains
{
	forAll(APTLRegions, zoneID)
	{
		#include "APTLSetFields.H"
		// Electrons
		dimensionedScalar condElectronEff_("CondElectronEff_", dimensionSet(condElectronEff.dimensions()), 1.);
		condElectronEff = condElectronEff_;
		// Hydrogen, Water Vapour
		dimensionedScalar diffEff_("diffEff_", dimensionSet(diffH2MixEff.dimensions()), 1.);
		dimensionedScalar density_("density_", dimensionSet(densityMix.dimensions()), 1.);
		transportCoeffH2.dimensions().reset(density_.dimensions()*diffEff_.dimensions());
		transportCoeffH2 = diffEff_*density_;
		diffH2MixEff = diffEff_;
		transportCoeffH2OVap.dimensions().reset(density_.dimensions()*diffEff_.dimensions());
		transportCoeffH2OVap = diffEff_*density_;
		diffH2OVapMixEff = diffEff_;
		// Energy
		dimensionedScalar condThermalEff_("condThermalEff_", dimensionSet(condThermalEff.dimensions()), 1.);
		condThermalEff = condThermalEff_;
		transportCoeffEnergy.dimensions().reset(condThermalEff_.dimensions());
		transportCoeffEnergy = condThermalEff_;
		// Mass
		dimensionedScalar dynamicPerm_("dynamicPerm_", dimensionSet(dynamicPerm.dimensions()), 1.);
		dynamicPerm = dynamicPerm_;
		// Liquid Water 
		dimensionedScalar diffH2OLiq_("diffH2OLiq_", dimensionSet(diffH2OLiq.dimensions()), 1.);
		diffH2OLiq = diffH2OLiq_;
	}
	forAll(ACLRegions, zoneID)
	{
		#include "ACLSetFields.H"
		// Electrons
		dimensionedScalar condElectronEff_("condElectronEff_", dimensionSet(condElectronEff.dimensions()), 1.);
		condElectronEff = condElectronEff_;
		// Protons
		dimensionedScalar condProtonEff_("condProtonEff_", dimensionSet(condProtonEff.dimensions()), 1.);
		condProtonEff = condProtonEff_;
		// Hydrogen, Water Vapour
		dimensionedScalar diffEff_("diffEff_", dimensionSet(diffH2MixEff.dimensions()), 1.);
		dimensionedScalar density_("density_", dimensionSet(densityMix.dimensions()), 1.);
		transportCoeffH2.dimensions().reset(density_.dimensions()*diffEff_.dimensions());
		transportCoeffH2 = diffEff_*density_;
		diffH2MixEff = diffEff_;
		transportCoeffH2OVap.dimensions().reset(density_.dimensions()*diffEff_.dimensions());
		transportCoeffH2OVap = diffEff_*density_;
		diffH2OVapMixEff = diffEff_;
		// Energy
		dimensionedScalar condThermalEff_("condThermalEff_", dimensionSet(condThermalEff.dimensions()), 1.);
		condThermalEff = condThermalEff_;
		transportCoeffEnergy.dimensions().reset(condThermalEff_.dimensions());
		transportCoeffEnergy = condThermalEff_;
		// Mass
		dimensionedScalar dynamicPerm_("dynamicPerm_", dimensionSet(dynamicPerm.dimensions()), 1.);
		dynamicPerm = dynamicPerm_;
		// Liquid Water 
		dimensionedScalar diffH2OLiq_("diffH2OLiq_", dimensionSet(diffH2OLiq.dimensions()), 1.);
		diffH2OLiq = diffH2OLiq_;
	}
	forAll(PMEMRegions, zoneID)
	{
		#include "PMEMSetFields.H"
		// Protons
		dimensionedScalar condProtonEff_("condProtonEff_", dimensionSet(condProtonEff.dimensions()), 1.);
		condProtonEff = condProtonEff_;
		// Energy
		dimensionedScalar condThermalEff_("condThermalEff_", dimensionSet(condThermalEff.dimensions()), 1.);
		condThermalEff = condThermalEff_;
		transportCoeffEnergy.dimensions().reset(condThermalEff_.dimensions());
		transportCoeffEnergy = condThermalEff_;
	}
	forAll(CCLRegions, zoneID)
	{
		#include "CCLSetFields.H"
		// Electrons
		dimensionedScalar condElectronEff_("condElectronEff_", dimensionSet(condElectronEff.dimensions()), 1.);
		condElectronEff = condElectronEff_;
		// Protons
		dimensionedScalar condProtonEff_("condProtonEff_", dimensionSet(condProtonEff.dimensions()), 1.);
		condProtonEff = condProtonEff_;
		// Oxygen, Water Vapour
		dimensionedScalar diffEff_("diffEff_", dimensionSet(diffO2MixEff.dimensions()), 1.);
		dimensionedScalar density_("density_", dimensionSet(densityMix.dimensions()), 1.);
		diffO2MixEff = diffEff_;
		transportCoeffO2.dimensions().reset(density_.dimensions()*diffEff_.dimensions());
		transportCoeffO2 = diffEff_*density_;
		transportCoeffH2OVap.dimensions().reset(density_.dimensions()*diffEff_.dimensions());
		transportCoeffH2OVap = diffEff_*density_;
		diffH2OVapMixEff = diffEff_;
		// Energy
		dimensionedScalar condThermalEff_("condThermalEff_", dimensionSet(condThermalEff.dimensions()), 1.);
		condThermalEff = condThermalEff_;
		transportCoeffEnergy.dimensions().reset(condThermalEff_.dimensions());
		transportCoeffEnergy = condThermalEff_;

		// Mass
		dimensionedScalar dynamicPerm_("dynamicPerm_", dimensionSet(dynamicPerm.dimensions()), 1.);
		dynamicPerm = dynamicPerm_;
		// Liquid Water 
		dimensionedScalar diffH2OLiq_("diffH2OLiq_", dimensionSet(diffH2OLiq.dimensions()), 1.);
		diffH2OLiq = diffH2OLiq_;
	}
	forAll(CPTLRegions, zoneID)
	{
		#include "CPTLSetFields.H"
		// Electrons
		dimensionedScalar condElectronEff_("condElectronEff_", dimensionSet(condElectronEff.dimensions()), 1.);
		condElectronEff = condElectronEff_;
		// Oxygen, Water Vapour
		dimensionedScalar diffEff_("diffEff_", dimensionSet(diffO2MixEff.dimensions()), 1.);
		dimensionedScalar density_("density_", dimensionSet(densityMix.dimensions()), 1.);
		diffO2MixEff = diffEff_;
		transportCoeffO2.dimensions().reset(density_.dimensions()*diffEff_.dimensions());
		transportCoeffO2 = diffEff_*density_;
		transportCoeffH2OVap.dimensions().reset(density_.dimensions()*diffEff_.dimensions());
		transportCoeffH2OVap = diffEff_*density_;
		diffH2OVapMixEff = diffEff_;
		// Energy
		dimensionedScalar condThermalEff_("condThermalEff_", dimensionSet(condThermalEff.dimensions()), 1.);
		condThermalEff = condThermalEff_;
		transportCoeffEnergy.dimensions().reset(condThermalEff_.dimensions());
		transportCoeffEnergy = condThermalEff_;
		// Mass
		dimensionedScalar dynamicPerm_("dynamicPerm_", dimensionSet(dynamicPerm.dimensions()), 1.);
		dynamicPerm = dynamicPerm_;
		// Liquid Water 
		dimensionedScalar diffH2OLiq_("diffH2OLiq_", dimensionSet(diffH2OLiq.dimensions()), 1.);
		diffH2OLiq = diffH2OLiq_;
	}
} // End setting up dummy coefficients
