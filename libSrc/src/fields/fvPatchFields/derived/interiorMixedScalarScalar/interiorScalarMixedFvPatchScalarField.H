/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2013-01-13 David B. Harvey:

Notes:
	- Re-wrote turbulentTemperatureCoupleBaffleMixed to apply an interior
	coupling for a given scalar field and tensor transport property
	- Re-wrote modified boundary condition to work in parallel
	- Moved modified boundary condition to operate in FAST namespace
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

Class
	Fastfc::interiorScalarMixedFvPatchScalarField

Description
    This is a mixed boundary condition intended for a scalar field with
	a tensor transport ceofficient.

	This follows the convention of the mixed boundary conditions as
	written for OpenFOAM's chtMultiRegionHeater tutorial

    Specifies gradient and potential to equate the flux on both sides
	- refGradient = zero gradient
    - refValue = neighbour value
    - mixFraction = nKDelta / (nKDelta + pKDelta)

    where pKDelta and nKDelta are the transport coefficient * inverse cell spacing

    Example usage:
        pInterfaceName_to_nInterfaceName
        {
            type        		FCAPOLLO::interiorMixedScalarScalar;
            transportCoeff		transportCoefficientName;
			value				initialization value for the face;
        }

    Needs to be on underlying mapped(Wall)FvPatch.

SourceFiles
    interiorScalarMixedFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef interiorScalarMixedFvPatchScalarField_H
#define interiorScalarMixedFvPatchScalarField_H

#include "interiorCoupledMixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Fastfc
{

/*---------------------------------------------------------------------------*\
   Class interiorScalarMixedFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class interiorScalarMixedFvPatchScalarField
:
    public interiorCoupledMixedFvPatchScalarField
{
    // Private data

        //- Name of field on the neighbour region
        const word transportCoeffName_;

public:

    //- Runtime type information
    TypeName("Fastfc::interiorMixedScalarScalar");


    // Constructors

        //- Construct from patch and internal field
        interiorScalarMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        interiorScalarMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  interiorScalarMixedFvPatchScalarField onto a
        //  new patch
        interiorScalarMixedFvPatchScalarField
        (
            const interiorScalarMixedFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new interiorScalarMixedFvPatchScalarField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        interiorScalarMixedFvPatchScalarField
        (
            const interiorScalarMixedFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new interiorScalarMixedFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Fastfc
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
