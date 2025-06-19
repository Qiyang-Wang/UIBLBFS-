/*---------------------------------------------------------------------------*\

    HiSA: High Speed Aerodynamic solver

    Copyright (C) 2014-2018 Oliver Oxtoby - CSIR, South Africa
    Copyright (C) 2014-2018 Johan Heyns - CSIR, South Africa

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

#include "discreteLaxFriedrich.H"
#include "emptyFvPatch.H"
#include "transformFvPatchFields.H"
#include "fvjOperators.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(discreteLaxFriedrich, 0);
addToRunTimeSelectionTable(jacobianInviscid, discreteLaxFriedrich, dictionary);


// * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * * * * //

discreteLaxFriedrich::discreteLaxFriedrich
(
    const dictionary& dict,
    const fvMesh& mesh,
    const volScalarField& rho,
    const volVectorField& rhoU,
    const volScalarField& rhoE,
    const psiThermo& thermo,
    const scalarField& ddtCoeff
)
 :
    jacobianInviscid(typeName, dict),
    mesh_(mesh),
    rho_(rho),
    rhoU_(rhoU),
    rhoE_(rhoE),
    thermo_(thermo),
    ddtCoeff_(ddtCoeff)
{}


// * * * * * * * * *  Private Member Functions  * * * * * * * * * * * * * * * //

void discreteLaxFriedrichJacobian::computeFluxJacobian()
void discreteLaxFriedrich::addFluxTerms(compressibleJacobianMatrix& jacobian)
{
    const surfaceScalarField& w(mesh_.weights());
    const volScalarField& rho(rho_);
    const volVectorField& rhoU(rhoU_);
    const volScalarField& rhoE(rhoE_);
    const volScalarField& gamma(gamma_());
    const surfaceScalarField& lambdaConv(lambdaConv_());
    const surfaceScalarField& phi =
        mesh_.lookupObject<surfaceScalarField>("phi");

    volVectorField U = rhoU/rho;

    // Convective part

    fvjMatrix<scalar>& dContByRho = jacobian.dSByS(0,0);
    fvjMatrix<vector>& dContByRhoU = jacobian.dSByV(0,0);
    //fvjMatrix<scalar>& dContByRhoE = jacobian.dSByS(0,1);
    fvjMatrix<vector>& dMomByRho = jacobian.dVByS(0,0);
    fvjMatrix<tensor>& dMomByRhoU = jacobian.dVByV(0,0);
    fvjMatrix<vector>& dMomByRhoE = jacobian.dVByS(0,1);
    fvjMatrix<scalar>& dEnergyByRho = jacobian.dSByS(1,0);
    fvjMatrix<vector>& dEnergyByRhoU = jacobian.dSByV(1,0);
    fvjMatrix<scalar>& dEnergyByRhoE = jacobian.dSByS(1,1);

    dContByRhoU += fvj::grad(w, geometricOneField());

    tmp<volScalarField> aStar = gamma-1;
    tmp<volScalarField> theta = 0.5*aStar()*magSqr(U);
    tmp<volTensorField> UU = U*U;
    surfaceScalarField wu = surfaceInterpolationScheme<vector>::New
    (
        mesh_,
        phi,
        mesh_.interpolationScheme("reconstruct(U)")
    )->weights(U);
    surfaceVectorField Uf = surfaceInterpolationScheme<vector>::interpolate(U, wu);

//    dMomByRho = fvj::div(w, -UU()) + fvj::grad(w, theta());
    dMomByRho = fvj::divMult(wu, phi, -U/rho) + fvj::grad(w, theta());
//    dMomByRhoU = fvj::div(w, U)*tensor::I + fvj::grad(w, U) + fvj::gradT(w, -aStar()*U);
    dMomByRhoU = fvj::divMult(wu, phi, 1.0/rho)*tensor::I + fvj::grad(w, Uf) + fvj::gradT(w, -aStar()*U);
//    dMomByRhoE = fvj::grad(w, aStar());

    tmp<volScalarField> gammaE = gamma*rhoE/rho;
    surfaceScalarField we = surfaceInterpolationScheme<scalar>::New
    (
        mesh_,
        phi,
        mesh_.interpolationScheme("reconstruct(T)")
    )->weights(rhoE/rho-0.5*magSqr(U));
    surfaceScalarField ef = surfaceInterpolationScheme<scalar>::interpolate(rhoE/rho-0.5*magSqr(U), we);
    surfaceVectorField UbyRhof = surfaceInterpolationScheme<vector>::interpolate(U/rho, we);
    surfaceScalarField wrho = surfaceInterpolationScheme<scalar>::New
    (
        mesh_,
        phi,
        mesh_.interpolationScheme("reconstruct(rho)")
    )->weights(rho);
    surfaceScalarField rhof = surfaceInterpolationScheme<scalar>::interpolate(rho, wrho);

//    dEnergyByRho = fvj::div(w, -gammaE()*U+2*theta()*U);
    dEnergyByRho = fvj::divMult(we, phi, -rhoE/sqr(rho)+magSqr(U)/rho)
                 + fvj::divDot(wu, phi*Uf, -U/rho)
                 + fvj::divMult(wrho, -phi/sqr(rhof), aStar()*rhoE-rho*theta())
                 + fvj::divMult(w, phi/rhof, theta());
//    dEnergyByRhoU = fvj::grad(w, gammaE()-theta()) + fvj::div(w, -aStar()*UU());
    dEnergyByRhoU = fvj::grad(w, ef+0.5*magSqr(Uf))
                  + fvj::divMult(we, phi, -U/rho)
                  + fvj::divMult(wu, phi*Uf, 1/rho)
                  + fvj::grad(w, 1/rhof, aStar()*rhoE-rho*theta())
                  + fvj::divMult(w, phi/rhof, -aStar()*U);
//    dEnergyByRhoE = fvj::div(w, gamma*U);
    dEnergyByRhoE = fvj::divMult(we, phi, 1/rho)
                  + fvj::divMult(w, phi/rhof, aStar());

    // Moving mesh part
    if (meshPhi_.valid())
    {
        const surfaceScalarField& meshPhi(meshPhi_()());
        tmp<fvjMatrix<scalar> > divMeshPhi = fvj::div(w, meshPhi);

        dContByRho -= divMeshPhi();
        dMomByRhoU -= divMeshPhi()*tensor::I;
        dEnergyByRhoE -= divMeshPhi;
    }

/*
        */

        // Preconditioning can't cope with arbitrary weights; need to use pure
        // Lax-Friedrich

        preContByRho = tmp<fvjMatrix<scalar> >(new fvjMatrix<scalar>(mesh_));
        preContByRhoU = fvj::grad(w, geometricOneField());
        preContByRhoE = tmp<fvjMatrix<scalar> >(new fvjMatrix<scalar>(mesh_));

        preMomByRho = fvj::div(w, -UU()) + fvj::grad(w, theta());
        preMomByRhoU = fvj::div(w, U)*tensor::I + fvj::grad(w, U) + fvj::gradT(w, -aStar()*U);
        preMomByRhoE = fvj::grad(w, aStar());

        preEnergyByRho = fvj::div(w, -gammaE()*U+2*theta()*U);
        preEnergyByRhoU = fvj::grad(w, gammaE-theta) + fvj::div(w, -aStar()*UU);
        preEnergyByRhoE = fvj::div(w, gamma*U);

        // Moving mesh part
        if (meshPhi_.valid())
        {
            const surfaceScalarField& meshPhi(meshPhi_()());
            tmp<fvjMatrix<scalar> > divMeshPhi = fvj::div(w, meshPhi);

            preContByRho() -= divMeshPhi();
            preMomByRhoU() -= divMeshPhi()*tensor::I;
            preEnergyByRhoE() -= divMeshPhi;
        }

    // Add Lax-Friedrich stabilisation term
    tmp<fvjMatrix<scalar> > stab = fvj::laplacian(0.5*lambdaConv, geometricOneField(), true);

    dContByRho -= stab();
    dMomByRhoU -= stab()*tensor::I;
    dEnergyByRhoE -= stab;
}


void discreteLaxFriedrich::boundaryJacobian
(
    label patchi,
    tmp<scalarField>& dContFluxdp, tmp<vectorField>& dContFluxdU, tmp<scalarField>& dContFluxdT,
    tmp<vectorField>& dMomFluxdp, tmp<tensorField>& dMomFluxdU, tmp<vectorField>& dMomFluxdT,
    tmp<scalarField>& dEnergyFluxdp, tmp<vectorField>& dEnergyFluxdU, tmp<scalarField>& dEnergyFluxdT
)
{
    // Compute approximate and full Jacobian diagonal matrices

    const fvMesh& mesh(mesh_);
    const volScalarField& rho(rho_);
    const volVectorField& rhoU(rhoU_);
    const volScalarField& rhoE(rhoE_);
    const volScalarField& p(thermo_.p());
    const volScalarField& T(thermo_.T());
    tmp< volScalarField > tcv = thermo_.Cv();
    const volScalarField& cv = tcv();

    const volScalarField::Boundary& pbf = p.boundaryField();
    const volVectorField::Boundary& ubf =

        mesh.lookupObject<volVectorField>("U").boundaryField();
    const volScalarField::Boundary& tbf = T.boundaryField();

    const scalarField& wP = mesh.weights().boundaryField()[patchi];

    // ValueInternalCoeffs
    // fixedValue   -> 0
    // zeroGradient -> 1
    tmp < vectorField > tuVIC = ubf[patchi].valueInternalCoeffs(wP);
    const vectorField& uVIC = tuVIC();
    tmp < scalarField > tpVIC = pbf[patchi].valueInternalCoeffs(wP);
    scalarField& pVIC = tpVIC.ref();
    tmp < scalarField > ttVIC = tbf[patchi].valueInternalCoeffs(wP);
    scalarField& tVIC = ttVIC.ref();

    const vectorField& SfB = mesh_.Sf().boundaryField()[patchi];
    const scalarField& rhoB = rho.boundaryField()[patchi];
    tmp<vectorField> UB = rhoU.boundaryField()[patchi]/rhoB;
    tmp<scalarField> UrelBdotSf = UB() & SfB;
    if (meshPhi_.valid())
    {
        UrelBdotSf.ref() -= meshPhi_()->boundaryField()[patchi];
    }
    const scalarField& pB = p.boundaryField()[patchi];
    const scalarField& TB = T.boundaryField()[patchi];
    const scalarField& rhoEB = rhoE.boundaryField()[patchi];
    const scalarField& cvB = cv.boundaryField()[patchi];

    dContFluxdp = rhoB/pB * UrelBdotSf() * pVIC;
    dContFluxdU = rhoB*cmptMultiply(SfB,uVIC);
    dContFluxdT = -rhoB/TB * UrelBdotSf() * tVIC;

    dMomFluxdp = (rhoB/pB*UB() * UrelBdotSf() + SfB) * pVIC;
    dMomFluxdU = rhoB * UB() * cmptMultiply(SfB,uVIC);
    tmp<vectorField> dMomFluxdUDiag = rhoB * UrelBdotSf() * uVIC;
    dMomFluxdT = -rhoB/TB*UB() * UrelBdotSf() * tVIC;

    dEnergyFluxdp = (rhoEB/pB*UrelBdotSf() + (UB() & SfB)) * pVIC;
    dEnergyFluxdU = cmptMultiply(SfB, uVIC)*(rhoEB+pB) + rhoB*UrelBdotSf()*cmptMultiply(UB(), uVIC);
    dEnergyFluxdT = UrelBdotSf*(rhoB*cvB-rhoEB/TB)*tVIC;

    dMomFluxdU.ref().replace(tensor::XX, dMomFluxdU().component(tensor::XX)+dMomFluxdUDiag().component(vector::X));
    dMomFluxdU.ref().replace(tensor::YY, dMomFluxdU().component(tensor::YY)+dMomFluxdUDiag().component(vector::Y));
    dMomFluxdU.ref().replace(tensor::ZZ, dMomFluxdU().component(tensor::ZZ)+dMomFluxdUDiag().component(vector::Z));
}


// Compute the temporal contribution to the approximate
// and full Jacobian diagonal matrices
void laxFriedrich::addTemporalTerms(compressibleJacobianMatrix& jacobian)
{
    const fvMesh& mesh(mesh_);
    const volScalarField& rho(rho_);
    const volVectorField& rhoU(rhoU_);

    volVectorField U = rhoU/rho;

    fvjMatrix<scalar>& dContByRho = jacobian.dSByS(0,0);
    fvjMatrix<tensor>& dMomByRhoU = jacobian.dVByV(0,0);
    fvjMatrix<scalar>& dEnergyByRhoE = jacobian.dSByS(1,1);

    // Temporal contribution (Weak form)
    scalarField diagCoeff = ddtCoeff_*mesh.V();
    dContByRho.diag() += diagCoeff;
    dMomByRhoU.diag() += diagCoeff*tensor::I;
    dEnergyByRhoE.diag() += diagCoeff;
}


// * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * * //

void discreteLaxFriedrich::addInviscidJacobian(compressibleJacobianMatrix& jacobian)
{
    gamma_ = thermo_.gamma();
    const volScalarField& gamma = gamma_();

    const volVectorField& U = mesh_.lookupObject<volVectorField>("U");

    if (mesh_.moving())
    {
        meshPhi_.reset(new tmp<surfaceScalarField>(fvc::meshPhi(U)));
    }
    else
    {
        meshPhi_.clear();
    }

    // Wave speed: Lax-Friedrich flux approximation of left-hand side Jacobian
    tmp< volScalarField > c = sqrt(gamma/thermo_.psi());
    if (mesh_.moving())
    {
        lambdaConv_ = fvc::interpolate(c) + mag((fvc::interpolate(U)&mesh_.Sf())-fvc::meshPhi(U))/mesh_.magSf();
    }
    else
    {
        lambdaConv_ = fvc::interpolate(c) + mag(fvc::interpolate(U)&mesh_.Sf()/mesh_.magSf());
    }
    // LU-SGS sweeps

    addFluxTerms(jacobian);
    addTemporalTerms(jacobian);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
