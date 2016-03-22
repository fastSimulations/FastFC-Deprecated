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
    Foam::interiorCoupledGenericFvPatchField

Group
    grpGenericBoundaryConditions

Description
    This boundary condition provides a base class for 'interiorCoupledGeneric' 
	type boundary conditions, i.e. conditions that couple two domains as though 
	the interface was an interior patch.  This base clase is derived based on a 
	a generic transport equation  with a control volume drawn around the patch.

    The respective equation for the value at the interior patch is determined by:

        \f[
            x_f = C_p x_p + C_n x_n
        \f]

    where
    \vartable
        x_f   | patch value
        x_p   | patch primary (internal) cell values
		x_n   | patch neighbour (internal) cell values		
        C_p   | Coefficient for the primary variable
		C_n   | Coefficient for the neighbour variable
    \endvartable

    \heading Patch usage

    \table
        Property      	  | Description            	| Required  | Default value
   		nbrCellVal     	  | neighbour cell value  	| yes       |
		ownerCoeffVal	  | weight field           	| yes       |
        nbrCoeffVal 	  | weight field			| yes		|
		nbrDeltaCoeffVal  | neighbour cell spacing	| yes		|
    \endtable

Note
    This condition is not usually applied directly; instead, use a derived
    interiorCoupledGeneric condition such as \c interiorAdvDiffScalarScalar 

SeeAlso
    Foam::fastfc::interiorAdvDiffScalarScalar

SourceFiles
    interiorCoupledGenericFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef interiorCoupledGenericFvPatchField_H
#define interiorCoupledGenericFvPatchField_H

#include "fvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class mixedFvPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class interiorCoupledGenericFvPatchField
:
    public fvPatchField<Type>
{
    // Private data

        //- Value field
//        Field<Type> refValue_;

        //- Normal gradient field
//        Field<Type> refGrad_;

        //- Fraction (0-1) of value used for boundary condition
//        scalarField valueFraction_;

	// New Variables
	
		//- Interior cell centre value of the neighbour field
		Field<Type> nbrCellVal_;

		//- Primary Cell Coefficient Value
		scalarField ownerCoeffVal_;
		
		//- Neighbour Cell Coefficient Value
		scalarField nbrCoeffVal_;

		//- Neighbour Cell Spacing Coefficients
		scalarField nbrDeltaCoeffVal_;

public:

    //- Runtime type information
    TypeName("interiorCoupledGeneric");


    // Constructors

        //- Construct from patch and internal field
        interiorCoupledGenericFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        interiorCoupledGenericFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping the given interiorCoupledGenericFvPatchField onto a new patch
        interiorCoupledGenericFvPatchField
        (
            const interiorCoupledGenericFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        interiorCoupledGenericFvPatchField
        (
            const interiorCoupledGenericFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
            (
                new interiorCoupledGenericFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        interiorCoupledGenericFvPatchField
        (
            const interiorCoupledGenericFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type> > clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type> >
            (
                new interiorCoupledGenericFvPatchField<Type>(*this, iF)
            );
        }

    // Member functions

        // Access

            //- Return true if this patch field fixes a value.
            //  Needed to check if a level has to be specified while solving
            //  Poissons equations.
            virtual bool fixesValue() const
            {
                return true;
            }


        // Return defining fields

//            virtual Field<Type>& refValue()
//            {
//                return refValue_;
//            }

//            virtual const Field<Type>& refValue() const
//            {
//                return refValue_;
//            }

//            virtual Field<Type>& refGrad()
//            {
//                return refGrad_;
//            }

//            virtual const Field<Type>& refGrad() const
//            {
//               return refGrad_;
//            }

//            virtual scalarField& valueFraction()
//            {
//                return valueFraction_;
//            }

//            virtual const scalarField& valueFraction() const
//            {
//                return valueFraction_;
//            }

		// New defining fields

			virtual Field<Type>& nbrCellVal()
			{
				return nbrCellVal_;
			}

			virtual const Field<Type>& nbrCellVal() const
			{
				return nbrCellVal_;
			}

			virtual scalarField& ownerCoeffVal()
			{
				return ownerCoeffVal_;
			}

			virtual const scalarField& ownerCoeffVal() const
			{
				return ownerCoeffVal_;
			}

			virtual scalarField& nbrCoeffVal()
			{
				return nbrCoeffVal_;
			}

			virtual const scalarField& nbrCoeffVal() const
			{
				return nbrCoeffVal_;
			}

			virtual scalarField& nbrDeltaCoeffVal()
			{
				return nbrDeltaCoeffVal_;
			}

			virtual const scalarField& nbrDeltaCoeffVal() const
			{
				return nbrDeltaCoeffVal_;
			}

		// Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<Type>&,
                const labelList&
            );


        // Evaluation functions

            //- Return gradient at boundary
            virtual tmp<Field<Type> > snGrad() const;

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType=Pstream::blocking
            );

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<Type> > valueInternalCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the value of this patchField with given weights
            virtual tmp<Field<Type> > valueBoundaryCoeffs
            (
                const tmp<scalarField>&
            ) const;

            //- Return the matrix diagonal coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<Type> > gradientInternalCoeffs() const;

            //- Return the matrix source coefficients corresponding to the
            //  evaluation of the gradient of this patchField
            virtual tmp<Field<Type> > gradientBoundaryCoeffs() const;


        //- Write
        virtual void write(Ostream&) const;


    // Member operators

        virtual void operator=(const UList<Type>&) {}

        virtual void operator=(const fvPatchField<Type>&) {}
        virtual void operator+=(const fvPatchField<Type>&) {}
        virtual void operator-=(const fvPatchField<Type>&) {}
        virtual void operator*=(const fvPatchField<scalar>&) {}
        virtual void operator/=(const fvPatchField<scalar>&) {}

        virtual void operator+=(const Field<Type>&) {}
        virtual void operator-=(const Field<Type>&) {}

        virtual void operator*=(const Field<scalar>&) {}
        virtual void operator/=(const Field<scalar>&) {}

        virtual void operator=(const Type&) {}
        virtual void operator+=(const Type&) {}
        virtual void operator-=(const Type&) {}
        virtual void operator*=(const scalar) {}
        virtual void operator/=(const scalar) {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "interiorCoupledGenericFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
