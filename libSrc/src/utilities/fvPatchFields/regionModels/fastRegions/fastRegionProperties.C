/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "objectRegistry.H"
#include "Time.H"
#include "fastRegionProperties.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fastRegionProperties::fastRegionProperties(const Time& runTime)
:
    IOdictionary
    (
        IOobject
        (
            "fastRegionsList",
            runTime.time().constant(),
            runTime.db(),
			IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    APTLRegionNames_(lookup("APTLRegionNames")),
	ACLRegionNames_(lookup("ACLRegionNames")),
    PMEMRegionNames_(lookup("PMEMRegionNames")),
    CCLRegionNames_(lookup("CCLRegionNames")),
	CPTLRegionNames_(lookup("CPTLRegionNames"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fastRegionProperties::~fastRegionProperties()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::List<Foam::word>& Foam::fastRegionProperties::APTLRegionNames() const
{
    return APTLRegionNames_;
}

const Foam::List<Foam::word>& Foam::fastRegionProperties::ACLRegionNames() const
{
    return ACLRegionNames_;
}

const Foam::List<Foam::word>& Foam::fastRegionProperties::PMEMRegionNames() const
{
    return PMEMRegionNames_;
}

const Foam::List<Foam::word>& Foam::fastRegionProperties::CCLRegionNames() const
{
    return CCLRegionNames_;
}

const Foam::List<Foam::word>& Foam::fastRegionProperties::CPTLRegionNames() const
{
	return CPTLRegionNames_;
}

// ************************************************************************* //
