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

Info<< "Reading thermophysical properties\n" << endl;
pThermo_ =
(
    psiThermo::New(mesh)
);
Info << endl;
psiThermo& thermo = pThermo_();

Cs_macro_.set(new dimensionedScalar("Cs", average(sqrt(1/thermo.psi()))));

Info << "Cs: " << Cs_macro_() << nl << endl;

volScalarField& p = thermo.p();

const volScalarField& mu = thermo.mu();

if (gMax(mu.internalField()) > 0.0)
{
    inviscid_ = false;
    Info << "Viscous analysis detected" << nl << endl;
}
else
{
    inviscid_ = true;
    Info << "Inviscid analysis detected" << nl << endl;
}

Info<< "Reading field U" << nl << endl;
U_.set
(
    new volVectorField
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    )
);
volVectorField& U = U_();

U_vol_.set  //新加的  指针内容实体化
(
    new volVectorField
    (
        IOobject
        (
            "U_vol",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    )
);
volVectorField& U_vol = U_vol_();

volTensorField gradU
(
   IOobject
   (
       "gradU",
       runTime.timeName(),
       mesh,
       IOobject::NO_READ,
       IOobject::NO_WRITE
   ),
   fvc::grad(U,"Gauss linear")
);

scalarVars_.setSize(1);
vectorVars_.setSize(1);

#include "rhoBoundaryTypes.H"
scalarVars_.set
(
    0,
    new volScalarField
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        thermo.rho(),
        rhoBoundaryTypes // NOTE: fixRho does a psi*p update (Check if working with nozzle BCs)
    )
);
volScalarField& rho = scalarVars_[0];

vectorVars_.set
(
    0,
    new volVectorField
    (
        IOobject
        (
            "rhoU",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        rho*U
    )
);
volVectorField& rhoU = vectorVars_[0];

phi_.set
(
    new surfaceScalarField
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh.Sf() & fvc::interpolate(rhoU)
    )
);
surfaceScalarField& phi = phi_();

Info<< "Creating turbulence model\n" << endl;
turbulence_ =
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);
compressibleTurbulenceModel& turbulence = turbulence_();
Info << endl;

Info<< "Creating flux scheme\n" << endl;
flux_ =
(
    fluxScheme_incomp::New(mesh.schemesDict(), thermo, rho, U, rhoU, turbulence) // from lyg: using new incompressible flux class
);
Info << endl;

rPseudoDeltaT_.set
(
    new volScalarField
    (
        IOobject
        (
            "rPseudoDeltaT",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("rPseudoDeltaT", dimless/dimTime, 0.0),
        zeroGradientFvPatchScalarField::typeName
    )
);

ddtCoeff_.set(new scalarField);

scalarResiduals_.setSize(1);
vectorResiduals_.setSize(1);

// Continuity residual
scalarResiduals_.set
(
    0,
    new volScalarField
    (
        IOobject
        (
            "rhoR",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", rho.dimensions()/dimTime, 0.0)
    )
);


// Mom Residual
vectorResiduals_.set
(
    0,
    new volVectorField
    (
        IOobject
        (
            "rhoUR",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", rhoU.dimensions()/dimTime, vector::zero)
    )
);


if (!inviscid_)
{
    // tauMC is stored for the benefit of the maxwellSlip BC - maybe better
    // to modify the BC to recalculate it rather
    tauMC_.set
    (
        new volTensorField
        (
            IOobject
            (
                "tauMC",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            turbulence_->muEff()*dev2(Foam::T(fvc::grad(U)))
        )
    );
}

rho.correctBoundaryConditions();
