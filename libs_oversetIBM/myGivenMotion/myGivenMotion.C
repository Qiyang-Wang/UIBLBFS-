/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "myGivenMotion.H"
#include "addToRunTimeSelectionTable.H"
#include "vectorIOList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(myGivenMotion, 0);
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        myGivenMotion,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::myGivenMotion::myGivenMotion
(
    const dictionary& SBMFCoeffs,
    const Time& runTime
)
:
    solidBodyMotionFunction(SBMFCoeffs, runTime)
{
    read(SBMFCoeffs);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::myGivenMotion::
~myGivenMotion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::myGivenMotion::transformation() const
{
    scalar t = time_.value();

    //const vector displacement = amplitude_*sin(omega_*t);
    vectorIOList test = time_.lookupObject<vectorIOList>("oversetMotion");
    const vector displacement(test[0].x(), test[0].y(), test[0].z());

    quaternion R(1);
    septernion TR(septernion(-displacement)*R);
    
    // 6DOF motion, an API would be helpful  
    // // vector trans(0, 0.5*sin(2*3.141592653*t/0.5), 0);
    // vector trans(0.5*sin(2*3.141592653*t/0.5), 0, 0);
    // vector rot(0, 0, 0.5*3.141592653*t);
    // vector CofG_(0, 0, 0);

    // quaternion R(quaternion::XYZ, rot);
    // septernion TR(septernion(-CofG_ + -trans)*R*septernion(CofG_));

    DebugInFunction << "Time = " << t << " transformation: " << TR << endl;

    return TR;
}


bool Foam::solidBodyMotionFunctions::myGivenMotion::read
(
    const dictionary& SBMFCoeffs
)
{
    solidBodyMotionFunction::read(SBMFCoeffs);

    return true;
}


// ************************************************************************* //
