/****************************/
//    constantsLibrary.H
/***************************/
#ifndef constantsLibrary_H
#define constantsLibrary_H

#include "dimensionedTypes.H"

// Nomenclature
// [ |kg|  m|  s|  K|mol|  A |cd|]

namespace Foam
{
	namespace Fast
	{
		namespace Constants	// Constants Library Namespace
		{

			const dimensionedScalar faraday
			(
				"faraday",
				dimensionSet( 0, 0,1, 0,-1, 1, 0 ),
				96485.          // Faraday's Constant (C/mol e-)
			);
			const dimensionedScalar gasConstantR
			(
				"gasConstantR",
				dimensionSet(  1, 2,-2,-1,-1, 0, 0 ), 
				8.314           // Universal Gas Constant (J/(mol*K))
			);
			const dimensionedScalar molWeightH2
			(
				"molWeightH2",
				dimensionSet(  1, 0, 0, 0,-1, 0, 0 ),
				2.002e-3        // Molecular Weight of H2 (kg/mol)
			);
			const dimensionedScalar molWeightO2
			(
				"molWeightO2",
				dimensionSet(  1, 0, 0, 0,-1, 0, 0 ),
				32.00e-3        // Molecular Weight of O2 (kg/mol)
			);
			const dimensionedScalar molWeightN2
			(
				"molWeightN2",
				dimensionSet(  1, 0, 0, 0,-1, 0, 0 ),
				28.01e-3        // Molecular Weight of N2 (kg/mol)
			);
			const dimensionedScalar	molWeightH2O
			(
				"molWeightH2O",
				dimensionSet(  1, 0, 0, 0,-1, 0, 0 ),
				18.02e-3        // Molecular Weight of H2O (kg/mol)
			);
			const dimensionedScalar	viscO2Ref
			(
				"viscO2Ref",
				dimensionSet(  1,-1,-1, 0, 0, 0, 0 ),
				2.4403e-5       // Reference Viscosity for O2 (kg/(m*s))
			);
			const dimensionedScalar viscH2Ref
			(
				"viscH2Ref",
				dimensionSet(  1,-1,-1, 0, 0, 0, 0 ),
				7.99842e-5      // Reference Viscosity for H2 (kg/(m*s))
			);
			const dimensionedScalar viscN2Ref
			(
				"viscN2Ref",
				dimensionSet(  1,-1,-1, 0, 0, 0, 0 ),
				2.1045e-5       // Reference Viscosity for N2 (kg/(m*s))
			);
			const dimensionedScalar viscH2OVapRef
			(
				"viscH2OVapRef",
				dimensionSet(  1,-1,-1, 0, 0, 0, 0 ),
				1.25e-5         // Reference Viscosity for H2O vapour  (kg/(m*s))
			);
			const dimensionedScalar viscH2OLiqRef
			(
				"viscH2OLiqRef",
				dimensionSet(  1,-1,-1, 0, 0, 0, 0 ),
				3.5e-4          // Reference Viscosity for H2O liquid (kg/(m*s))
			);
			const dimensionedScalar diffVolO2
			(
				"diffVolO2",
				dimensionSet(  0, 3, 0, 0, 0, 0, 0 ),
				16.6            // Atomic diffusion volume for O2 (reference: Cussler) (m^3)
			);
			const dimensionedScalar diffVolN2
			(
				"diffVolN2",
				dimensionSet(  0, 3, 0, 0, 0, 0, 0 ),
				17.9            // Atomic diffusion volume for N2 (reference: Cussler) (m^3)
			);
			const dimensionedScalar diffVolH2OVap
			(
				"diffVolH2OVap",
				dimensionSet(  0, 3, 0, 0, 0, 0, 0 ),
				12.7            // Atomic diffusion volume for H2O vapour (reference: Cussler) (m^3)
			);
			const dimensionedScalar diffVolH2
			(
				"diffVolH2",
				dimensionSet(  0, 3, 0, 0, 0, 0, 0 ),
				7.07            // Atomic diffusion volume for H2 (reference: Cussler) (m^3)
			);
			const dimensionedScalar surfTensH2OLiq
			(
				"surfTensH2OLiq",
				dimensionSet(  1, 0,-2, 0, 0, 0, 0 ),
				6.25e-2         // Surface Tension of Liquid water at ~75 C
			);
			const dimensionedScalar pStd
			(
				"pStd",
				dimensionSet(  1,-1,-2, 0, 0, 0, 0 ),
				1.01325e5
			);
			const dimensionedScalar densityH2OLiq
			(
				"rho_H2O_l",
				dimensionSet (  1,-3, 0, 0, 0, 0, 0 ),
				1000.
			);
			const dimensionedScalar latentHeatH2OVap
			(
				"latentHeatH2OVap",
				dimensionSet ( 0, 2,-2, 0, 0, 0, 0 ),
				2.26e6
			);
		}; // end namespace constants
	}; // end namespace Fast
} // end namespace Foam
#endif
