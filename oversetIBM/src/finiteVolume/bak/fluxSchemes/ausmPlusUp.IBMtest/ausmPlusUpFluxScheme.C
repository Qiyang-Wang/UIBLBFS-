/*---------------------------------------------------------------------------*\

    HiSA: High Speed Aerodynamic solver

    Copyright (C) 2014-2018 Oliver Oxtoby - CSIR, South Africa
    Copyright (C) 2014-2018 Johan Heyns - CSIR, South Africa
    Copyright (C) 1991-2008 OpenCFD Ltd.

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

#include "ausmPlusUpFluxScheme.H"
#include "addToRunTimeSelectionTable.H"
#include "bound.H"
#include "fvcSurfaceReconstruct.H"
#include "cellFaceFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(ausmPlusUpFluxScheme, 0);
addToRunTimeSelectionTable(fluxScheme, ausmPlusUpFluxScheme, dictionary);


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

ausmPlusUpFluxScheme::ausmPlusUpFluxScheme
(
    const dictionary& dict,
    const psiThermo& thermo,
    const volScalarField& rho,
    const volVectorField& U,
    const volVectorField& rhoU,
    const volScalarField& rhoE,
    const compressibleTurbulenceModel& turbulence
)
:
    fluxScheme(typeName, dict),
    mesh_(U.mesh()),
    thermo_(thermo),
    rho_(rho),
    U_(U),
    rhoU_(rhoU),
    rhoE_(rhoE),
    dict_(dict),
    turbulence_(turbulence)
{}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

ausmPlusUpFluxScheme::~ausmPlusUpFluxScheme()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::ausmPlusUpFluxScheme::calcFlux(surfaceScalarField& phi, surfaceVectorField& phiUp, surfaceScalarField& phiEp, surfaceVectorField& Up)
{
    bool lowMach =
    dict_.lookupOrDefault<Switch>
    (
        "lowMachAusm",
        true
    );

    const volScalarField& p = thermo_.p();

    tmp<surfaceVectorField> U_L, U_R;
    fvc::surfaceReconstruct(U_, U_L, U_R, "reconstruct(U)");
    oversetInterpolateV(U_, U_L, U_R);

    tmp<surfaceScalarField> phi_L = U_L()&mesh_.Sf();
    tmp<surfaceScalarField> phi_R = U_R()&mesh_.Sf();

    // Flux relative to mesh movement
    if (mesh_.moving())
    {
        fvc::makeRelative(phi_L.ref(), U_);
        fvc::makeRelative(phi_R.ref(), U_);
    }
    surfaceScalarField un_L = phi_L/mesh_.magSf();
    surfaceScalarField un_R = phi_R/mesh_.magSf();

    // Critical acoustic velocity (Liou 2006)
    tmp< volScalarField > gamma = thermo_.gamma();
    //tmp< volScalarField > H((rhoE_+p)/rho_);
    tmp< volScalarField > H
    (
        (max(rhoE_/rho_,dimensionedScalar("0", rhoE_.dimensions()/rho_.dimensions(), SMALL)) +
         max(p/rho_,dimensionedScalar("0", p.dimensions()/rho_.dimensions(), SMALL)))
    );
    H->rename("H");

    if (mesh_.moving())
    {
        H.ref() -= 0.5*(U_&U_);
        volVectorField Urel(U_);
        Urel -= fvc::reconstruct(fvc::meshPhi(U_));
        H.ref() += 0.5*(Urel&Urel);
    }
//    H().rename("H");
//    bound(H(), dimensionedScalar("H0", H().dimensions(), 0.0));

    tmp< volScalarField > c_tmp = sqrt(2.0*(gamma()-1.0)/(gamma()+1.0)*H());
    c_tmp->rename("c");
    const volScalarField& c = c_tmp();
    gamma.clear();

    tmp<surfaceScalarField> c_L, c_R;
    fvc::surfaceReconstruct(c, c_L, c_R, "reconstruct(T)");
    oversetInterpolate(c, c_L, c_R);
    
    c_L = sqr(c_L())/max(c_L(), un_L);
    c_R = sqr(c_R())/max(c_R(),-un_R);
    tmp< surfaceScalarField > c_face(min(c_L(),c_R()));
    c_L.clear();
    c_R.clear();

    // Critical Mach number
    tmp< surfaceScalarField > Mach_L(un_L/c_face());

    // Split Mach numbers
    tmp<surfaceScalarField> Mach_plus_L = 
        calcForEachFace
        (
            [](const scalar& MLf)
            {  
                if (mag(MLf) < 1.0)
                {
                    scalar ML2p =  0.25*sqr(MLf+1);
                    scalar ML2m = -0.25*sqr(MLf-1);
                    //return ML2p;             // beta = 0
                    return ML2p*(1 - 2*ML2m);  // beta = 1/8
                }
                else
                {
                    return max(MLf, 0);
                }
            },
            Mach_L()
        );

    // Pressure flux
    tmp<surfaceScalarField> p_plus_L =
        calcForEachFace
        (
            [](const scalar& MLf)
            {
                if (mag(MLf) < 1.0)
                {
                    scalar ML2p =  0.25*sqr(MLf+1);
                    scalar ML2m = -0.25*sqr(MLf-1);
                    return ML2p*(2 - MLf - 3*MLf*ML2m);  //alpha = 3/16
                }
                else
                {
                    return (MLf > 0 ? 1.0 : 0.0);
                }
            },
            Mach_L()
        );
    Mach_L.clear();

    tmp< surfaceScalarField > Mach_R(un_R/c_face());

    // Split Mach numbers
    tmp<surfaceScalarField> Mach_minus_R = 
        calcForEachFace
        (
            [](const scalar& MRf)
            {  
                if (mag(MRf) < 1.0)
                {
                    scalar MR2m = -0.25*sqr(MRf-1);
                    scalar MR2p =  0.25*sqr(MRf+1);
                    //return MR2m;             // beta = 0
                    return MR2m*(1 + 2*MR2p);  // beta = 1/8
                }
                else
                {
                    return min(MRf, 0);
                }
            },
            Mach_R()
        );

    // Pressure flux
    tmp<surfaceScalarField> p_minus_R =
        calcForEachFace
        (
            [](const scalar& MRf)
            {
                if (mag(MRf) < 1.0)
                {
                    scalar MR2m = -0.25*sqr(MRf-1);
                    scalar MR2p =  0.25*sqr(MRf+1);
                    return MR2m*(-2 - MRf + 3*MRf*MR2p);  //alpha = 3/16
                }
                else
                {
                    return (MRf < 0 ? 1.0 : 0.0);
                }
            },
            Mach_R()
        );
    Mach_R.clear();

    tmp<surfaceScalarField> p_L, p_R;
    fvc::surfaceReconstruct(p, p_L, p_R, "reconstruct(rho)");
    oversetInterpolate(p, p_L, p_R);

    tmp<surfaceScalarField> rho_L, rho_R;
    fvc::surfaceReconstruct(rho_, rho_L, rho_R, "reconstruct(rho)");
    oversetInterpolate(rho_, rho_L, rho_R);

    tmp< surfaceScalarField > Mach_1_2 = Mach_plus_L() + Mach_minus_R();
    Mach_plus_L.clear();
    Mach_minus_R.clear();

    surfaceScalarField p_1_2 = p_plus_L()*p_L() + p_minus_R()*p_R();

    // Low Mach number diffusive term
    tmp< surfaceScalarField > M_mean = 0.5*(sqr(un_L) + sqr(un_R))/sqr(c_face());   // Mean local Mach number
    tmp< surfaceScalarField > MDiff = -0.25*max((1.0-M_mean()),0.0)*(p_R() - p_L())/(0.5*(rho_L()+rho_R())*sqr(c_face()));    // 0 < K_p < 1, Liou suggest 0.25
    M_mean.clear();

    // Mach_1_2.ref() += MDiff();
    Mach_1_2.ref() +=
        calcForEachFace
        (
            [](const scalar& Mach_1_2, const scalar& MDiff)
            {
                if 
                (
                    (Mach_1_2 > 0.0 && Mach_1_2 + MDiff <= 0.0) ||
                    (Mach_1_2 < 0.0 && Mach_1_2 + MDiff >= 0.0)
                )
                {
                    return 0.2*MDiff;
                }
                else
                {
                    return MDiff;
                }
            },
            Mach_1_2(),
            MDiff()
        );

    MDiff.clear();

    // Low Mach number diffusive term
    if(lowMach)
    {
        tmp< surfaceScalarField > pDiff = -0.25*p_plus_L()*p_minus_R()*(rho_L()+rho_R())*c_face()*(un_R-un_L); // 0 < Ku < 1; Liou suggests 0.75
        #ifdef OPENFOAM_PLUS
            pDiff->setOriented(false);
        #endif
        p_1_2 +=  pDiff();
        pDiff.clear();
    }

    #ifdef OPENFOAM_PLUS
        Mach_1_2->setOriented(true);
    #endif

    p_L.clear();
    p_R.clear();
    p_plus_L.clear();
    p_minus_R.clear();

    tmp<surfaceVectorField> U_f = surfaceFieldSelect(U_L, U_R, Mach_1_2(), 0);

    surfaceScalarField rhoa_LR  = Mach_1_2()*c_face()*surfaceFieldSelect(rho_L, rho_R, Mach_1_2(), 0);
    surfaceVectorField rhoaU_LR = rhoa_LR*U_f();
    // volScalarField ee("ee",rhoE_/rho_-0.5*magSqr(U_));
    // surfaceScalarField rhoah_LR = rhoa_LR*(fvc::surfaceReconstruct(ee, Mach_1_2(), "reconstruct(T)") + 0.5*magSqr(U_f())) + p_1_2*Mach_1_2()*c_face();
    // NOTE: According to Liou enthalpy should be interpolated.
    //surfaceScalarField rhoah_LR = rhoa_LR*(fvc::surfaceReconstruct(H(), Mach_1_2(), "reconstruct(T)"));

//    volScalarField& Href = H();
    surfaceScalarField tmp12 = (fvc::surfaceReconstruct(H(), Mach_1_2(), "reconstruct(T)"));
    oversetInterpolate(H, tmp12);
    surfaceScalarField rhoah_LR = rhoa_LR*tmp12;

    // Face velocity for sigmaDotU (turbulence term)
    Up = U_f()*mesh_.magSf();
    U_f.clear();
    H.clear();
//    c_face.clear();

    phi = rhoa_LR*mesh_.magSf();
    phiUp = rhoaU_LR*mesh_.magSf() + p_1_2*mesh_.Sf();
    phiEp = rhoah_LR*mesh_.magSf();

    if (mesh_.moving())
    {
        //phiEp += fvc::meshPhi(U_)*fvc::surfaceReconstruct(p, Mach_1_2(), "reconstruct(T)");
        //phiEp += fvc::meshPhi(U_)*fvc::surfaceReconstruct(rho_, Mach_1_2(), "reconstruct(rho)")*fvc::surfaceReconstruct(p/rho_, Mach_1_2(), "reconstruct(T)"); //Ensure consistent interpolation with pressure term above
        //Not neccessary if using second last line
        phiEp += p_1_2 * fvc::meshPhi(U_); //Ensure consistent interpolation with pressure term above
    }
}

void Foam::ausmPlusUpFluxScheme::oversetInterpolate
(
    const volScalarField& vvl,
    tmp<surfaceScalarField>& vsf_L_tmp,
    tmp<surfaceScalarField>& vsf_R_tmp
)
{
    surfaceScalarField& vsf_L = vsf_L_tmp.ref();
    surfaceScalarField& vsf_R = vsf_R_tmp.ref();

    oversetInterpolate(vvl, vsf_L, vsf_R);
}

void Foam::ausmPlusUpFluxScheme::oversetInterpolate
(
    const volScalarField& vvl,
    surfaceScalarField& vsf_L,
    surfaceScalarField& vsf_R
)
{
    const surfaceScalarField &faceType = mesh_.lookupObject<surfaceScalarField>("faceType");
    const surfaceScalarField &faceDonor = mesh_.lookupObject<surfaceScalarField>("faceDonor");

    //Info << faceType << endl;
    //cin.get();
    //Info << faceDonor << endl;
    //cin.get();

    for (label fi=0; fi<mesh_.nInternalFaces(); fi++)
    {
        if (faceType[fi] == 1)
        {
            vsf_L[fi] = vvl[faceDonor[fi]];
            vsf_R[fi] = vsf_L[fi];
        }
    }

    forAll(mesh_.boundary(), pi)
    {
        fvsPatchField<scalar>& vsfp_L = vsf_L.boundaryFieldRef()[pi];
        fvsPatchField<scalar>& vsfp_R = vsf_R.boundaryFieldRef()[pi];
        const fvsPatchField<scalar>& faceDonorp = faceDonor.boundaryField()[pi];

        forAll(mesh_.boundary()[pi], fi)
        {
            if (faceType.boundaryField()[pi][fi] == 1)
            {
                vsfp_L[fi] = vvl[faceDonorp[fi]];
                vsfp_R[fi] = vsfp_L[fi];
            }
        }
    }
}

void Foam::ausmPlusUpFluxScheme::oversetInterpolateV
(
    const volVectorField& vvl,
    tmp<surfaceVectorField>& vsf_L_tmp,
    tmp<surfaceVectorField>& vsf_R_tmp
)
{
    surfaceVectorField& vsf_L = vsf_L_tmp.ref();
    surfaceVectorField& vsf_R = vsf_R_tmp.ref();

    oversetInterpolateV(vvl, vsf_L, vsf_R);
}

void Foam::ausmPlusUpFluxScheme::oversetInterpolateV
(
    const volVectorField& vvl,
    surfaceVectorField& vsf_L,
    surfaceVectorField& vsf_R
)
{
    const surfaceScalarField &faceType = mesh_.lookupObject<surfaceScalarField>("faceType");
    const surfaceScalarField &faceDonor = mesh_.lookupObject<surfaceScalarField>("faceDonor");

    //Info << faceType << endl;
    //cin.get();
    //Info << faceDonor << endl;
    //cin.get();

    for (label fi=0; fi<mesh_.nInternalFaces(); fi++)
    {
        if (faceType[fi] == 1)
        {
            vsf_L[fi] = vvl[faceDonor[fi]];
            vsf_R[fi] = vsf_L[fi];
        }
    }

    forAll(mesh_.boundary(), pi)
    {
        fvsPatchField<vector>& vsfp_L = vsf_L.boundaryFieldRef()[pi];
        fvsPatchField<vector>& vsfp_R = vsf_R.boundaryFieldRef()[pi];
        const fvsPatchField<scalar>& faceDonorp = faceDonor.boundaryField()[pi];

        forAll(mesh_.boundary()[pi], fi)
        {
            if (faceType.boundaryField()[pi][fi] == 1)
            {
                vsfp_L[fi] = vvl[faceDonorp[fi]];
                vsfp_R[fi] = vsfp_L[fi];
            }
        }
    }
}

void Foam::ausmPlusUpFluxScheme::oversetInterpolate
(
    const volScalarField& vvl,
    tmp<surfaceScalarField>& vsf_tmp
)
{
    surfaceScalarField& vsf = vsf_tmp.ref();

    oversetInterpolate(vvl, vsf);
}

void Foam::ausmPlusUpFluxScheme::oversetInterpolate
(
    const volScalarField& vvl,
    surfaceScalarField& vsf
)
{
    const surfaceScalarField &faceType = mesh_.lookupObject<surfaceScalarField>("faceType");
    const surfaceScalarField &faceDonor = mesh_.lookupObject<surfaceScalarField>("faceDonor");

    for (label fi=0; fi<mesh_.nInternalFaces(); fi++)
    {
        if (faceType[fi] == 1)
        {
            vsf[fi] = vvl[faceDonor[fi]];
        }
    }

    forAll(mesh_.boundary(), pi)
    {
        fvsPatchField<scalar>& vsfp = vsf.boundaryFieldRef()[pi];
        const fvsPatchField<scalar>& faceDonorp = faceDonor.boundaryField()[pi];

        forAll(mesh_.boundary()[pi], fi)
        {
            if (faceType.boundaryField()[pi][fi] == 1)
            {
                vsfp[fi] = vvl[faceDonorp[fi]];
		Info << mesh_.Cf().boundaryField()[pi][fi] << endl;
		Info << mesh_.C()[faceDonorp[fi]] << endl;
		cin.get();
            }
        }
    }
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
