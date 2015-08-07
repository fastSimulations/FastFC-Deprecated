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

Class
    threeBlockRegionProperties

Description
    Simple class to hold region information for coupled region simulations
    based on a physics naming convention.

SourceFiles
    threeBlockRegionProperties.C

\*---------------------------------------------------------------------------*/

#ifndef threeBlockRegionProperties_H
#define threeBlockRegionProperties_H

#include "IOdictionary.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class regionProperties Declaration
\*---------------------------------------------------------------------------*/

class threeBlockRegionProperties
:
    public IOdictionary
{
    // Private data

        //- List of the region names
        List<word> ACLRegionNames_;
		List<word> PMEMRegionNames_;
		List<word> CCLRegionNames_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        threeBlockRegionProperties(const threeBlockRegionProperties&);

        //- Disallow default bitwise assignment
        void operator=(const threeBlockRegionProperties&);

public:

    // Constructors

        //- Construct from components
        threeBlockRegionProperties(const Time& runTime);


    //- Destructor
    ~threeBlockRegionProperties();


    // Member Functions

        // Access

            //- Return const reference to the list of region names
            const List<word>& ACLRegionNames() const;
	   		const List<word>& PMEMRegionNames() const;
		    const List<word>& CCLRegionNames() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
