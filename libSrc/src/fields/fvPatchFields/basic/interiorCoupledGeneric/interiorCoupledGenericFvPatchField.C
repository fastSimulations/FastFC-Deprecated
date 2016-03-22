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

#include "interiorCoupledGenericFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
interiorCoupledGenericFvPatchField<Type>::interiorCoupledGenericFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(p, iF),
//    refValue_(p.size()),
//    refGrad_(p.size()),
//    valueFraction_(p.size()),
	nbrCellVal_(p.size()),
	ownerCoeffVal_(p.size()),
	nbrCoeffVal_(p.size()),
	nbrDeltaCoeffVal_(p.size())
{}


template<class Type>
interiorCoupledGenericFvPatchField<Type>::interiorCoupledGenericFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fvPatchField<Type>(p, iF, dict),
//    refValue_("refValue", dict, p.size()),
//    refGrad_("refGradient", dict, p.size()),
//    valueFraction_("valueFraction", dict, p.size()),
	nbrCellVal_("nbrCellVal", dict, p.size()),
	ownerCoeffVal_("ownerCoeffVal", dict, p.size()),
	nbrCoeffVal_("nbrCoeffVal", dict, p.size()),
	nbrDeltaCoeffVal_("nbrDeltaCoeffVal", dict, p.size())
{
    evaluate();
}


template<class Type>
interiorCoupledGenericFvPatchField<Type>::interiorCoupledGenericFvPatchField
(
    const interiorCoupledGenericFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fvPatchField<Type>(ptf, p, iF, mapper),
//    refValue_(ptf.refValue_, mapper),
//    refGrad_(ptf.refGrad_, mapper),
//    valueFraction_(ptf.valueFraction_, mapper),
	nbrCellVal_(ptf.nbrCellVal_, mapper),
	ownerCoeffVal_(ptf.ownerCoeffVal_, mapper),
	nbrCoeffVal_(ptf.nbrCoeffVal_, mapper),
	nbrDeltaCoeffVal_(ptf.nbrDeltaCoeffVal_, mapper)
{}


template<class Type>
interiorCoupledGenericFvPatchField<Type>::interiorCoupledGenericFvPatchField
(
    const interiorCoupledGenericFvPatchField<Type>& ptf
)
:
    fvPatchField<Type>(ptf),
//    refValue_(ptf.refValue_),
//    refGrad_(ptf.refGrad_),
//    valueFraction_(ptf.valueFraction_),
	nbrCellVal_(ptf.nbrCellVal_),
	ownerCoeffVal_(ptf.ownerCoeffVal_),
	nbrCoeffVal_(ptf.nbrCoeffVal_),
	nbrDeltaCoeffVal_(ptf.nbrDeltaCoeffVal_)
{}


template<class Type>
interiorCoupledGenericFvPatchField<Type>::interiorCoupledGenericFvPatchField
(
    const interiorCoupledGenericFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(ptf, iF),
//    refValue_(ptf.refValue_),
//    refGrad_(ptf.refGrad_),
//    valueFraction_(ptf.valueFraction_),
	nbrCellVal_(ptf.nbrCellVal_),
	ownerCoeffVal_(ptf.ownerCoeffVal_),
	nbrCoeffVal_(ptf.nbrCoeffVal_),
	nbrDeltaCoeffVal_(ptf.nbrDeltaCoeffVal_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void interiorCoupledGenericFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fvPatchField<Type>::autoMap(m);
//    refValue_.autoMap(m);
//    refGrad_.autoMap(m);
//    valueFraction_.autoMap(m);
	nbrCellVal_.autoMap(m);
	ownerCoeffVal_.autoMap(m);
	nbrCoeffVal_.autoMap(m);
	nbrDeltaCoeffVal_.autoMap(m);
}


template<class Type>
void interiorCoupledGenericFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fvPatchField<Type>::rmap(ptf, addr);

    const interiorCoupledGenericFvPatchField<Type>& mptf =
        refCast<const interiorCoupledGenericFvPatchField<Type> >(ptf);

//    refValue_.rmap(mptf.refValue_, addr);
//    refGrad_.rmap(mptf.refGrad_, addr);
//    valueFraction_.rmap(mptf.valueFraction_, addr);
	nbrCellVal_.rmap(mptf.nbrCellVal_, addr);
	ownerCoeffVal_.rmap(mptf.ownerCoeffVal_, addr);
	nbrCoeffVal_.rmap(mptf.nbrCoeffVal_, addr);
	nbrDeltaCoeffVal_.rmap(mptf.nbrDeltaCoeffVal_, addr);
}


template<class Type>
void interiorCoupledGenericFvPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    Field<Type>::operator=
    (
        nbrCoeffVal_*nbrCellVal_
      +
        ownerCoeffVal_*this->patchInternalField()
    );

    fvPatchField<Type>::evaluate();
}


template<class Type>
tmp<Field<Type> > interiorCoupledGenericFvPatchField<Type>::snGrad() const
{
    return
       (nbrCellVal_ - this->patchInternalField())
       *(this->patch().deltaCoeffs() + nbrDeltaCoeffVal_);
//       valueFraction_
//      *(refValue_ - this->patchInternalField())
//      *this->patch().deltaCoeffs()
//     + (1.0 - valueFraction_)*refGrad_;
}


template<class Type>
tmp<Field<Type> > interiorCoupledGenericFvPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    return pTraits<Type>::one*ownerCoeffVal_;
//    return pTraits<Type>::one*(1.0 - valueFraction_);

}


template<class Type>
tmp<Field<Type> > interiorCoupledGenericFvPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    return
         nbrCoeffVal_*nbrCellVal_;
//    return
//         valueFraction_*refValue_
//       + (1.0 - valueFraction_)*refGrad_/this->patch().deltaCoeffs();
}


template<class Type>
tmp<Field<Type> > interiorCoupledGenericFvPatchField<Type>::gradientInternalCoeffs() const
{
    return 
		-pTraits<Type>::one*nbrCoeffVal_*this->patch().deltaCoeffs();
//    return -pTraits<Type>::one*valueFraction_*this->patch().deltaCoeffs();
}


template<class Type>
tmp<Field<Type> > interiorCoupledGenericFvPatchField<Type>::gradientBoundaryCoeffs() const
{
    return
        nbrCoeffVal_*this->patch().deltaCoeffs()*nbrCellVal_;
//    return
//        valueFraction_*this->patch().deltaCoeffs()*refValue_
//      + (1.0 - valueFraction_)*refGrad_;
}


template<class Type>
void interiorCoupledGenericFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
	nbrCellVal_.writeEntry("nbrCellVal", os);
	ownerCoeffVal_.writeEntry("ownerCoeffVal", os);
	nbrCoeffVal_.writeEntry("nbrCoeffVal", os);
	nbrDeltaCoeffVal_.writeEntry("nbrDeltaCoeffVal", os);
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
