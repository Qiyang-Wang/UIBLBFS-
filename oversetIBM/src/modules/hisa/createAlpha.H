/*---------------------------------------------------------------------------*\

    HiSA: High Speed Aerodynamic solver

    Copyright (C) 2014-2018 Johan Heyns - CSIR, South Africa
    Copyright (C) 2014-2018 Oliver Oxtoby - CSIR, South Africa
    Copyright (C) 2011-2013 OpenFOAM Foundation

-------------------------------------------------------------------------------
License
    This file is part of HiSA.

    HiSA is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    HiSA is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with HiSA.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

Info<< "Reading porosity properties\n" << endl;

const IOdictionary dictPor
(
    IOobject
    (
        "porosityProperties",
        mesh_().time().caseConstant(),
        mesh_(),
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false // Do not register
    )
);

epsilon_ = dictPor.lookupOrDefault<scalar>("epsilon", 1);//如果不存在，则返回默认值1
eta_ = dictPor.lookupOrDefault<scalar>("eta", 30);
scalar KK = dictPor.lookupOrDefault<scalar>("K", 1e-4);
K_ = dimensionedScalar("K", dimensionSet(0,2,0,0,0,0,0), KK);
C_F = 1.75/(sqrt(150*epsilon_*epsilon_*epsilon_));


Info << "hahahahaha " << epsilon_  << ", " << eta_ << ", " << K_ << endl;

alpha_.set
(
    new volScalarField
    (
    IOobject
    (
        "alpha",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("myZero",dimless, 0.0),
    zeroGradientFvPatchScalarField::typeName
    )
);
volScalarField& alpha = alpha_();


GAMMA_.set
(
    new volScalarField
    (
    IOobject
    (
        "GAMMA",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("myZero",dimless, 0.0),
    zeroGradientFvPatchScalarField::typeName
    )
);
volScalarField& GAMMA = GAMMA_();


// initialise alpha field 
const fvMesh& mesh_whole = mesh_();
const labelIOList& zoneID_tmp = (mesh_whole.lookupObject<labelIOList>("zoneIDList"));

oversetZoneID_ = 1;
scalar r, z;

for(int j = 0; j<mesh.C().size(); j++)
{
    if(zoneID_tmp[j]!=oversetZoneID_)
    {
        alpha[j] = 1;

        GAMMA[j] = 0;
    }
    else
    {
        r = mag(mesh.C()[j] - vector(0, 0, 0));
        z = r - 0.5;

        alpha[j] = (1 + epsilon_)/2 + (1 - epsilon_)/2* tanh(eta_*z);

        if(z < 0)
        {
            GAMMA[j] = 1;
        }


    }

    // if (alpha[j] < 1 && zoneID_tmp[j]!=oversetZoneID_)
    // {
    //     Info << "hahahahaha " << zoneID_tmp[j] << ", " << alpha[j] << endl;
    // }

}

alpha.correctBoundaryConditions();
GAMMA.correctBoundaryConditions();

// creat grad alpha field
gradAlpha_.set
(
    new volVectorField
    (
    IOobject
    (
        "grad_alpha",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::grad(alpha,"Gauss linear")/alpha // divided by alpha for convenience, see Equations in CH's thesis 
    )
);
volVectorField& gradAlpha = gradAlpha_();

// Darcy-Forchheimer force
F_r_.set
(
    new volVectorField
    (
    IOobject
    (
        "Fr",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    rho/sqrt(K_)*U * mag(U),
    zeroGradientFvPatchVectorField::typeName
    )
);
volVectorField& F_r = F_r_();

