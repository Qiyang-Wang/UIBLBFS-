/*---------------------------------------------------------------------------*\

    HiSA: High Speed Aerodynamic solver

    Copyright (C) 2018 Johan Heyns - CSIR, South Africa
    Copyright (C) 2018 Oliver Oxtoby - CSIR, South Africa

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

#include "laplacian.H"
#include "emptyFvPatch.H"
#include "transformFvPatchFields.H"
#include "fvjOperators.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(laplacian, 0);
addToRunTimeSelectionTable(jacobianViscous, laplacian, dictionary);


// * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * * * * //

laplacian::laplacian
(
    const dictionary& dict,
    const fvMesh& mesh,
    const volScalarField& rho,
    const volVectorField& rhoU,
    const volScalarField& rhoE,
    const psiThermo& thermo,
    const compressible::turbulenceModel& turbulence
)
 :
    jacobianViscous(typeName, dict),
    mesh_(mesh),
    rho_(rho),
    rhoU_(rhoU),
    rhoE_(rhoE),
    thermo_(thermo),
    turbulence_(turbulence)
{}


// * * * * * * * * *  Private Member Functions  * * * * * * * * * * * * * * * //

// Assemble full matrix with diagonal and off-diagonal contribution from the
// flux terms.
void laplacian::addFluxTerms(compressibleJacobianMatrix& jacobian)
{
    fvjMatrix<scalar>& dContByRho = jacobian.dSByS(0,0);
    fvjMatrix<tensor>& dMomByRhoU = jacobian.dVByV(0,0);
    fvjMatrix<scalar>& dEnergyByRhoE = jacobian.dSByS(1,1);

    // Viscous part
    // Add laplacian term similar to Lax-Friedrich stabilisation term
    tmp<fvjMatrix<scalar> > stab = fvj::laplacian(0.5*lambdaVisc_, geometricOneField(), true);

    dContByRho -= stab();
    dMomByRhoU -= stab()*tensor::I;
    dEnergyByRhoE -= stab;
}

void laplacian::boundaryJacobian
(
    label patchi,
    tmp<scalarField>& dContFluxdp, tmp<vectorField>& dContFluxdU, tmp<scalarField>& dContFluxdT,
    tmp<vectorField>& dMomFluxdp, tmp<tensorField>& dMomFluxdU, tmp<vectorField>& dMomFluxdT,
    tmp<scalarField>& dEnergyFluxdp, tmp<vectorField>& dEnergyFluxdU, tmp<scalarField>& dEnergyFluxdT
)
{
    dContFluxdp = tmp<scalarField>(new scalarField(mesh_.boundary()[patchi].patch().size(), Zero));
    dContFluxdU = tmp<vectorField>(new vectorField(mesh_.boundary()[patchi].patch().size(), Zero));
    dContFluxdT = tmp<scalarField>(new scalarField(mesh_.boundary()[patchi].patch().size(), Zero));
    dMomFluxdp = tmp<vectorField>(new vectorField(mesh_.boundary()[patchi].patch().size(), Zero));
    dMomFluxdU = tmp<tensorField>(new tensorField(mesh_.boundary()[patchi].patch().size(), Zero));
    dMomFluxdT = tmp<vectorField>(new vectorField(mesh_.boundary()[patchi].patch().size(), Zero));
    dEnergyFluxdp = tmp<scalarField>(new scalarField(mesh_.boundary()[patchi].patch().size(), Zero));
    dEnergyFluxdU = tmp<vectorField>(new vectorField(mesh_.boundary()[patchi].patch().size(), Zero));
    dEnergyFluxdT = tmp<scalarField>(new scalarField(mesh_.boundary()[patchi].patch().size(), Zero));
}


// * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * * //

void laplacian::addViscousJacobian(compressibleJacobianMatrix& jacobian)
{
    alphaEff_ = fvc::interpolate(turbulence_.alphaEff());
    muEff_ = fvc::interpolate(turbulence_.muEff());
    const surfaceScalarField& alphaEff = alphaEff_();
    const surfaceScalarField& muEff = muEff_();
    gamma_ = thermo_.gamma();

    // Viscous wave speed (Luo et al, 2001)
    // CHECK: Blazek proposes an alternative approximation of viscous spectrial
    // radius
    lambdaVisc_ = (muEff+alphaEff)/fvc::interpolate(rho_)*mesh_.deltaCoeffs();
//    lambda_() += (muEff)/fvc::interpolate(rho_)*mesh_.deltaCoeffs(); // Luo et al. (2000) do not include alphaEff

    addFluxTerms(jacobian);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
