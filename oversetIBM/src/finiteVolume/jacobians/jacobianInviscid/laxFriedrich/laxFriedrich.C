/*---------------------------------------------------------------------------*\

    HiSA: High Speed Aerodynamic solver

    Copyright (C) 2014-2018 Johan Heyns - CSIR, South Africa
    Copyright (C) 2014-2018 Oliver Oxtoby - CSIR, South Africa

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

#include "laxFriedrich.H"
#include "emptyFvPatch.H"
#include "transformFvPatchFields.H"
#include "fvjOperators.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(laxFriedrich, 0);
addToRunTimeSelectionTable(jacobianInviscid, laxFriedrich, dictionary);


// * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * * * * //

laxFriedrich::laxFriedrich
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

// Assemble full matrix with diagonal and off-diagonal contribution from the
// flux terms.
void laxFriedrich::addFluxTerms(compressibleJacobianMatrix& jacobian)
{
    #ifdef OPENFOAM_PLUS
        surfaceScalarField w(mesh_.weights());
        w.setOriented(false);  // Workaround; it should not be oriented
    #else
        const surfaceScalarField& w(mesh_.weights());
    #endif
    const volScalarField& rho(rho_);
    const volVectorField& rhoU(rhoU_);
    const volScalarField& rhoE(rhoE_);
    const volScalarField& gamma(gamma_());
    const surfaceScalarField& lambdaConv(lambdaConv_());

    volVectorField U = rhoU/rho;

    // Intermediate variables
    tmp<volScalarField> aStar = gamma-1;
    tmp<volScalarField> theta = 0.5*aStar()*(U&U);
    tmp<volTensorField> UU = U*U;
    tmp<volScalarField> gammaE = gamma*rhoE/rho;

    fvjMatrix<scalar>& dContByRho = jacobian.dSByS(0,0);
    fvjMatrix<vector>& dContByRhoU = jacobian.dSByV(0,0);
    //fvjMatrix<scalar>& dContByRhoE = jacobian.dSByS(0,1);
    fvjMatrix<vector>& dMomByRho = jacobian.dVByS(0,0);
    fvjMatrix<tensor>& dMomByRhoU = jacobian.dVByV(0,0);
    fvjMatrix<vector>& dMomByRhoE = jacobian.dVByS(0,1);
    fvjMatrix<scalar>& dEnergyByRho = jacobian.dSByS(1,0);
    fvjMatrix<vector>& dEnergyByRhoU = jacobian.dSByV(1,0);
    fvjMatrix<scalar>& dEnergyByRhoE = jacobian.dSByS(1,1);

    // Diagonal and off-diagonal contribution of convective part
    dContByRhoU += fvj::grad(w, geometricOneField());

    dMomByRho += fvj::div(w, -UU()) + fvj::grad(w, theta());
    dMomByRhoU += fvj::div(w, U)*tensor::I + fvj::grad(w, U) + fvj::gradT(w, -aStar()*U);
    dMomByRhoE += fvj::grad(w, aStar());

    dEnergyByRho += fvj::div(w, -gammaE()*U+2*theta()*U);
    dEnergyByRhoU += fvj::grad(w, gammaE-theta) + fvj::div(w, -aStar*UU);
    dEnergyByRhoE += fvj::div(w, gamma*U);

    // Moving mesh part
    if (meshPhi_.valid())
    {
        const surfaceScalarField& meshPhi(meshPhi_()());
        tmp<fvjMatrix<scalar> > divMeshPhi = fvj::div(w, meshPhi);

        dContByRho -= divMeshPhi();
        dMomByRhoU -= divMeshPhi()*tensor::I;
        dEnergyByRhoE -= divMeshPhi;
    }

    // Add Lax-Friedrich stabilisation term
    tmp<fvjMatrix<scalar> > stab = fvj::laplacian(0.5*lambdaConv, geometricOneField(), true);

    dContByRho -= stab();
    dMomByRhoU -= stab()*tensor::I;
    dEnergyByRhoE -= stab;
}


void laxFriedrich::boundaryJacobian
(
    label patchi,
    tmp<scalarField>& dContFluxdp, tmp<vectorField>& dContFluxdU, tmp<scalarField>& dContFluxdT,
    tmp<vectorField>& dMomFluxdp, tmp<tensorField>& dMomFluxdU, tmp<vectorField>& dMomFluxdT,
    tmp<scalarField>& dEnergyFluxdp, tmp<vectorField>& dEnergyFluxdU, tmp<scalarField>& dEnergyFluxdT
)
{
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
    const scalarField& pVIC = tpVIC();
    tmp < scalarField > ttVIC = tbf[patchi].valueInternalCoeffs(wP);
    const scalarField& tVIC = ttVIC();

    const vectorField& SfB = mesh_.Sf().boundaryField()[patchi];
    const scalarField& rhoB = rho.boundaryField()[patchi];
    tmp<vectorField> UB = rhoU.boundaryField()[patchi]/rhoB;
    tmp<scalarField> UrelBdotSf = UB() & SfB;
    if (meshPhi_.valid())
    {
        const surfaceScalarField& meshPhi(meshPhi_()());
        UrelBdotSf.ref() -= meshPhi.boundaryField()[patchi];
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

void laxFriedrich::addInviscidJacobian(compressibleJacobianMatrix& jacobian)
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

    addFluxTerms(jacobian);
    addTemporalTerms(jacobian);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
