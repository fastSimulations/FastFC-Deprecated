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

Class
	Fastfc::interiorConductionScalarFvPatchScalarField

Description
    This is a boundary condition that takes an external boundary and 
	changes it to an interior boundary condition for multi-domain 
	simulations.  It is written for use with a scalar transport coefficient
	(i.e. conductivity) and a scalar solved variable (i.e. scalar field). 

    This boundary condition sets the value of the solved field on the boundary 
	face based on flux conservation across the boundary between the two physical
	domains.  Both domains must be set to use the same boundary condition for the 
	system of equations to operate correctly.

	Parameters:
	- nIntValue = Value of the solved field at the neighbour cell centre
	- pCoeffValue = coefficient term for the primary cell centre
	- nCoeffValue = coefficient term for the neighbour cell centre 
	- nDeltaCoeffValue = Delta values for the neighbour cell centre to the boundary 

    Example usage:
        pInterfaceName_to_nInterfaceName
        {
            type        		Fastfc::interiorConductionScalarScalar;
            ConductivityField	Name of Conductivity Field;
        }

    This boundary condition needs to be applied on a underlying mapped(Wall)FvPatch.

SourceFiles
    interiorConductionScalarFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef interiorConductionScalarFvPatchScalarField_H
#define interiorConductionScalarFvPatchScalarField_H

#include "interiorCoupledGenericFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Fastfc
{

/*---------------------------------------------------------------------------*\
   Class interiorConductionScalarFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class interiorConductionScalarFvPatchScalarField
:
    public interiorCoupledGenericFvPatchScalarField
{
    // Private data

        //- Name of field on the neighbour region
        const word conductivityFieldName_;

public:

    //- Runtime type information
    TypeName("Fastfc::interiorConductionScalarScalar");


    // Constructors

        //- Construct from patch and internal field
        interiorConductionScalarFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        interiorConductionScalarFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  interiorConductionScalarFvPatchScalarField onto a
        //  new patch
        interiorConductionScalarFvPatchScalarField
        (
            const interiorConductionScalarFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new interiorConductionScalarFvPatchScalarField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        interiorConductionScalarFvPatchScalarField
        (
            const interiorConductionScalarFvPatchScalarField&,
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
                new interiorConductionScalarFvPatchScalarField
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
