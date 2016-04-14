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
	This file is part of FAST-FC.

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
\*---------------------------------------------------------------------------*/

#include "interiorAdvectionDiffusionScalarFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "directMappedPatchBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Fastfc
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

interiorAdvectionDiffusionScalarFvPatchScalarField::
interiorAdvectionDiffusionScalarFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    interiorCoupledGenericFvPatchScalarField(p, iF),
    diffusivityFieldName_("undefined-diffusivityField")
{
    this->nbrCellVal() = 0.0;
    this->nbrCoeffVal() = 1.0;
	this->ownerCoeffVal() = 0.0;
}


interiorAdvectionDiffusionScalarFvPatchScalarField::
interiorAdvectionDiffusionScalarFvPatchScalarField
(
    const interiorAdvectionDiffusionScalarFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    interiorCoupledGenericFvPatchScalarField(ptf, p, iF, mapper),
    diffusivityFieldName_(ptf.diffusivityFieldName_)
{}


interiorAdvectionDiffusionScalarFvPatchScalarField::
interiorAdvectionDiffusionScalarFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    interiorCoupledGenericFvPatchScalarField(p, iF),
    diffusivityFieldName_(dict.lookup("diffusivityField"))
{
    if (!isA<directMappedPatchBase>(this->patch().patch()))
    {
        FatalErrorIn
        (
            "interiorAdvectionDiffusionScalarFvPatchScalarField::"
            "interiorAdvectionDiffusionScalarFvPatchScalarField\n"
            "(\n"
            "    const fvPatch& p,\n"
            "    const DimensionedField<scalar, volMesh>& iF,\n"
            "    const dictionary& dict\n"
            ")\n"
        )   << "\n    patch type '" << p.type()
            << "' not type '" << directMappedPatchBase::typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << dimensionedInternalField().name()
            << " in file " << dimensionedInternalField().objectPath()
            << exit(FatalError);
    }

    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    if (dict.found("nbrCellVal"))
    {
        // Full restart
        nbrCellVal() = scalarField("nbrCellVal", dict, p.size());
		nbrCoeffVal() = scalarField("nbrCoeffVal", dict, p.size());
		ownerCoeffVal() = scalarField("ownerCoeffVal", dict, p.size());
    }
    else
    {
        // Start from user entered data. Assume fixedValue.
        nbrCellVal() = *this;
		nbrCoeffVal() = 1.0;	
		ownerCoeffVal() = 0.0;
    }
}


interiorAdvectionDiffusionScalarFvPatchScalarField::
interiorAdvectionDiffusionScalarFvPatchScalarField
(
    const interiorAdvectionDiffusionScalarFvPatchScalarField& wtcsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    interiorCoupledGenericFvPatchScalarField(wtcsf, iF),
    diffusivityFieldName_(wtcsf.diffusivityFieldName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void interiorAdvectionDiffusionScalarFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Get the coupling information from the directMappedPatchBase
    const directMappedPatchBase& mpp = 
			refCast<const directMappedPatchBase>
    		(
        		patch().patch()
    		);
	
    const polyMesh& neighbourMesh = mpp.sampleMesh();
    const fvPatch& neighbourPatch = 
			refCast<const fvMesh>
    		(
        		neighbourMesh
			).boundary()[mpp.samplePolyPatch().index()];

    // Force recalculation of mapping and schedule
    const mapDistribute& distMap = mpp.map();


	// Determine the unit normal vectors for the principal and neighbour patches
	vectorField pUnitNormal(this->patch().Sf()/this->patch().magSf());
	vectorField nUnitNormal(neighbourPatch.Sf()/neighbourPatch.magSf());

	mapDistribute::distribute
    (
        Pstream::defaultCommsType,
        distMap.schedule(),
        distMap.constructSize(),
        distMap.subMap(),           // what to send
        distMap.constructMap(),     // what to receive
        nUnitNormal
    );
	
	// Lookup and assign internal solved field values for the principal and neighbour patches
	const fvPatchScalarField& pFieldPhi =
			refCast<const fvPatchScalarField>
			(
				patch().lookupPatchField<volScalarField, scalar>
				(
				 	dimensionedInternalField().name()
				)
			);
	const fvPatchScalarField& nFieldPhi = 
			refCast<const fvPatchScalarField>
			(
				neighbourPatch.lookupPatchField<volScalarField, scalar>
				(
				 	dimensionedInternalField().name()
				)
			);

	// Lookup and assign diffusivity coefficient field values for the principal and neighbour patches
	const fvPatchScalarField& pFieldDiffusivity = 
			refCast<const fvPatchScalarField>
			(
				patch().lookupPatchField<volScalarField, scalar>
				(
			 		diffusivityFieldName_
				)
			);

	const fvPatchScalarField& nFieldDiffusivity = 
			refCast<const fvPatchScalarField>
			(
				neighbourPatch.lookupPatchField<volScalarField, scalar>
				(
				 	diffusivityFieldName_
				)
			);

	// Lookup and assign mixture density field values for the principal and neighbour patches
	const fvPatchScalarField& pFieldDensityMix = 
			refCast<const fvPatchScalarField>
			(
				patch().lookupPatchField<volScalarField, scalar>
				(
					"densityMix"
				)
			);

	const fvPatchScalarField& nFieldDensityMix =
			refCast<const fvPatchScalarField>
			(
				patch().lookupPatchField<volScalarField, scalar>
				(
					"densityMix"
				)
			);

	// Lookup and assign mixture velocity field values for the principal and neighbour patches
	const fvPatchScalarField& pFieldVelocityAvgMix =
			refCast<const fvPatchScalarField>
			(
				patch().lookupPatchField<volScalarField, scalar>
				(
					"velAvg"
				)
			);

	const fvPatchScalarField& nFieldVelocityAvgMix = 
			refCast<const fvPatchScalarField>
			(
				patch().lookupPatchField<volScalarField, scalar>
				(
					"velAvg"
				)
			);

	// Swap to obtain full local values for the solved field
	scalarField pIntFldPhi(pFieldPhi.patchInternalField());
	scalarField nIntFldPhi(nFieldPhi.patchInternalField());

    mapDistribute::distribute
    (
        Pstream::defaultCommsType,
        distMap.schedule(),
        distMap.constructSize(),
        distMap.subMap(),           // what to send
        distMap.constructMap(),     // what to receive
        nIntFldPhi
    );

	// Swap to obtain full local values for the diffusivity coefficient
	scalarField pIntFldDiffusivity(pFieldDiffusivity.patchInternalField());
	scalarField nIntFldDiffusivity(nFieldDiffusivity.patchInternalField());

    mapDistribute::distribute
    (
        Pstream::defaultCommsType,
        distMap.schedule(),
        distMap.constructSize(),
        distMap.subMap(),           // what to send
        distMap.constructMap(),     // what to receive
        nIntFldDiffusivity
    );

	// Swap to obtain full local values for the mixture density
	scalarField pIntFldDensityMix(pFieldDensityMix.patchInternalField());
	scalarField nIntFldDensityMix(nFieldDensityMix.patchInternalField());

    mapDistribute::distribute
    (
        Pstream::defaultCommsType,
        distMap.schedule(),
        distMap.constructSize(),
        distMap.subMap(),           // what to send
        distMap.constructMap(),     // what to receive
        nIntFldDensityMix
    );

	// Swap to obtain full local values for the mixture average velocity
	scalarField pIntFldVelocityAvgMix(pFieldVelocityAvgMix.patchInternalField());
	scalarField nIntFldVelocityAvgMix(nFieldVelocityAvgMix.patchInternalField());

    mapDistribute::distribute
    (
        Pstream::defaultCommsType,
        distMap.schedule(),
        distMap.constructSize(),
        distMap.subMap(),           // what to send
        distMap.constructMap(),     // what to receive
        nIntFldVelocityAvgMix
    );

	// Inverse cell spacing
	scalarField pDelta(patch().deltaCoeffs());
	scalarField nDelta(neighbourPatch.deltaCoeffs());

    mapDistribute::distribute
    (
        Pstream::defaultCommsType,
        distMap.schedule(),
        distMap.constructSize(),
        distMap.subMap(),           // what to send
        distMap.constructMap(),     // what to receive
        nDelta
    );

	// Calculate Coefficients
	scalarField tauP = pIntFldDensityMix*pIntFldDiffusivity*pDelta;
	scalarField tauN = nIntFldDensityMix*nIntFldDiffusivity*nDelta;
	scalarField betaP = pIntFldDensityMix*pIntFldVelocityAvgMix;
	scalarField betaN = nIntFldDensityMix*nIntFldVelocityAvgMix;

	// Check for zero or non-zero coefficient and variables
	// Diffusivity Check
	tauP = Foam::max(tauP, VSMALL);
	tauN = Foam::max(tauN, VSMALL);

	// Determine the flux at the boundary patch based on a weighted
	// average of the principal and neighbour cells
	this->nbrCellVal() = nIntFldPhi;

	this->ownerCoeffVal() = (tauP + betaP)  / (tauP + tauN);

	this->nbrCoeffVal() = (tauN + betaN)  / (tauP + tauN);

	this->nbrDeltaCoeffVal() = nDelta;

	interiorCoupledGenericFvPatchScalarField::updateCoeffs();
}


void interiorAdvectionDiffusionScalarFvPatchScalarField::write
(
    Ostream& os
) const
{
    interiorCoupledGenericFvPatchScalarField::write(os);
	os.writeKeyword("diffusivityField") << diffusivityFieldName_ 
		<< token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    interiorAdvectionDiffusionScalarFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Fastfc
} // End namespace Foam


// ************************************************************************* //
