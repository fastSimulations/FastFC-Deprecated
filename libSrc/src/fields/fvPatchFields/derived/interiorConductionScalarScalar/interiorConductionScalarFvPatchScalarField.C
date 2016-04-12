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

#include "interiorConductionScalarFvPatchScalarField.H"
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

interiorConductionScalarFvPatchScalarField::
interiorConductionScalarFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    interiorCoupledGenericFvPatchScalarField(p, iF),
    conductivityFieldName_("undefined-ConductivityFieldName")
{
    this->nbrCellVal() = 0.0;
    this->nbrCoeffVal() = 1.0;
	this->ownerCoeffVal() = 0.0;
}


interiorConductionScalarFvPatchScalarField::
interiorConductionScalarFvPatchScalarField
(
    const interiorConductionScalarFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    interiorCoupledGenericFvPatchScalarField(ptf, p, iF, mapper),
    conductivityFieldName_(ptf.conductivityFieldName_)
{}


interiorConductionScalarFvPatchScalarField::
interiorConductionScalarFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    interiorCoupledGenericFvPatchScalarField(p, iF),
    conductivityFieldName_(dict.lookup("conductivityField"))
{
    if (!isA<directMappedPatchBase>(this->patch().patch()))
    {
        FatalErrorIn
        (
            "interiorConductionScalarFvPatchScalarField::"
            "interiorConductionScalarFvPatchScalarField\n"
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


interiorConductionScalarFvPatchScalarField::
interiorConductionScalarFvPatchScalarField
(
    const interiorConductionScalarFvPatchScalarField& wtcsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    interiorCoupledGenericFvPatchScalarField(wtcsf, iF),
    conductivityFieldName_(wtcsf.conductivityFieldName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void interiorConductionScalarFvPatchScalarField::updateCoeffs()
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
			 		conductivityFieldName_
				)
			);

	const fvPatchScalarField& nFieldK = 
			refCast<const fvPatchScalarField>
			(
				neighbourPatch.lookupPatchField<volScalarField, scalar>
				(
				 	conductivityFieldName_
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

	// Check for zero or non-zero transport coefficient
	pIntFldK = Foam::max(pIntFldK, SMALL);
	nIntFldK = Foam::max(nIntFldK, SMALL);

	// Apply a cell spacing based weighting on the conductivities
	scalarField pKDelta(pIntFldK*pDelta);
	scalarField nKDelta(nIntFldK*nDelta);

	// Determine the flux at the boundary patch based on a weighted
	// average of the principal and neighbour cells
	this->nbrCellVal() = nIntFldPhi;

	this->nbrCoeffVal() = (nKDelta  / (nKDelta + pKDelta));

	this->ownerCoeffVal() = (pKDelta  / (nKDelta + pKDelta));

	this->nbrDeltaCoeffVal() = nDelta;

	interiorCoupledGenericFvPatchScalarField::updateCoeffs();
}


void interiorConductionScalarFvPatchScalarField::write
(
    Ostream& os
) const
{
    interiorCoupledGenericFvPatchScalarField::write(os);
	os.writeKeyword("conductivityField") << conductivityFieldName_ 
		<< token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    interiorConductionScalarFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Fastfc
} // End namespace Foam


// ************************************************************************* //
