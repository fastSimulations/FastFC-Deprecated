/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
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

\*---------------------------------------------------------------------------*/

#include "interiorScalarMixedFvPatchScalarField.H"
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

interiorScalarMixedFvPatchScalarField::
interiorScalarMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    interiorCoupledMixedFvPatchScalarField(p, iF),
    transportCoeffName_("undefined-transportCoeffName")
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 1.0;
}


interiorScalarMixedFvPatchScalarField::
interiorScalarMixedFvPatchScalarField
(
    const interiorScalarMixedFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    interiorCoupledMixedFvPatchScalarField(ptf, p, iF, mapper),
    transportCoeffName_(ptf.transportCoeffName_)
{}


interiorScalarMixedFvPatchScalarField::
interiorScalarMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    interiorCoupledMixedFvPatchScalarField(p, iF),
    transportCoeffName_(dict.lookup("transportCoeff"))
{
    if (!isA<directMappedPatchBase>(this->patch().patch()))
    {
        FatalErrorIn
        (
            "interiorScalarMixedFvPatchScalarField::"
            "interiorScalarMixedFvPatchScalarField\n"
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

    if (dict.found("refValue"))
    {
        // Full restart
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        // Start from user entered data. Assume fixedValue.
        refValue() = *this;
        refGrad() = 0.0;
        valueFraction() = 1.0;
    }
}


interiorScalarMixedFvPatchScalarField::
interiorScalarMixedFvPatchScalarField
(
    const interiorScalarMixedFvPatchScalarField& wtcsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    interiorCoupledMixedFvPatchScalarField(wtcsf, iF),
    transportCoeffName_(wtcsf.transportCoeffName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void interiorScalarMixedFvPatchScalarField::updateCoeffs()
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
	
	// Lookup and assign internal field values for the principal and neighbour patches
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

	// Lookup and assign transport coefficient field values for the principal and neighbour patches
	const fvPatchScalarField& pFieldK = 
			refCast<const fvPatchScalarField>
			(
				patch().lookupPatchField<volScalarField, scalar>
				(
			 		transportCoeffName_
				)
			);

	const fvPatchScalarField& nFieldK = 
			refCast<const fvPatchScalarField>
			(
				neighbourPatch.lookupPatchField<volScalarField, scalar>
				(
				 	transportCoeffName_
				)
			);

	// Swap to obtain full local values for the neighbouring field values
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

	// transport coefficient
	scalarField pIntFldK(pFieldK.patchInternalField());
	scalarField nIntFldK(nFieldK.patchInternalField());

    mapDistribute::distribute
    (
        Pstream::defaultCommsType,
        distMap.schedule(),
        distMap.constructSize(),
        distMap.subMap(),           // what to send
        distMap.constructMap(),     // what to receive
        nIntFldK
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

	// Apply a cell spacing based weighting on the conductivities
	scalarField pKDelta(pIntFldK*pDelta);
	scalarField nKDelta(nIntFldK*nDelta);

	// Determine the flux at the boundary patch based on a weighted
	// average of the principal and neighbour cells
	this->refValue() = nIntFldPhi;

	this->refGrad() = 0.;

	this->valueFraction() = (nKDelta  / (nKDelta + pKDelta));


	interiorCoupledMixedFvPatchScalarField::updateCoeffs();
}


void interiorScalarMixedFvPatchScalarField::write
(
    Ostream& os
) const
{
    interiorCoupledMixedFvPatchScalarField::write(os);
	os.writeKeyword("transportCoeff") << transportCoeffName_ 
		<< token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    interiorScalarMixedFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Fastfc
} // End namespace Foam


// ************************************************************************* //
