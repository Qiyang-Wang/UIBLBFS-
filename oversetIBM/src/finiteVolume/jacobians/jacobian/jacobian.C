/*---------------------------------------------------------------------------*\

    HiSA: High Speed Aerodynamic solver

    Copyright (C) 2017 Oliver Oxtoby - CSIR, South Africa
    Copyright (C) 2017 Johan Heyns - CSIR, South Africa

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

#include "jacobian.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(jacobian, 0);

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

void jacobian::addBoundaryTerms()
{
    const fvMesh& mesh(mesh_);
    const volScalarField& rho(rho_);
    const volVectorField& rhoU(rhoU_);
    const volScalarField& rhoE(rhoE_);
    tmp< volScalarField > tcv = thermo_.Cv();
    const volScalarField& cv = tcv();
    tmp< volScalarField > tgamma(thermo_.gamma());
    const volScalarField& gamma = tgamma();

    volVectorField U = rhoU/rho;

    fvjMatrix<scalar>& dContByRho = jacobianMatrix_.dSByS(0,0);
    fvjMatrix<vector>& dContByRhoU = jacobianMatrix_.dSByV(0,0);
    fvjMatrix<scalar>& dContByRhoE = jacobianMatrix_.dSByS(0,1);
    fvjMatrix<vector>& dMomByRho = jacobianMatrix_.dVByS(0,0);
    fvjMatrix<tensor>& dMomByRhoU = jacobianMatrix_.dVByV(0,0);
    fvjMatrix<vector>& dMomByRhoE = jacobianMatrix_.dVByS(0,1);
    fvjMatrix<scalar>& dEnergyByRho = jacobianMatrix_.dSByS(1,0);
    fvjMatrix<vector>& dEnergyByRhoU = jacobianMatrix_.dSByV(1,0);
    fvjMatrix<scalar>& dEnergyByRhoE = jacobianMatrix_.dSByS(1,1);

    // Calculate diagonal contribution of the geometric boundaries

    forAll(mesh.boundaryMesh(), patchi)
    {

        const fvPatch& patch = mesh.boundary()[patchi];

        if(!patch.coupled() && !isA<emptyFvPatch>(patch))
        {

            tmp<scalarField> dContFluxdp;
            tmp<vectorField> dContFluxdU;
            tmp<scalarField> dContFluxdT;

            tmp<vectorField> dMomFluxdp;
            tmp<tensorField> dMomFluxdU;
            tmp<vectorField> dMomFluxdT;

            tmp<scalarField> dEnergyFluxdp;
            tmp<vectorField> dEnergyFluxdU;
            tmp<scalarField> dEnergyFluxdT;

            jacobianInviscid_->boundaryJacobian
            (
                patchi,
                dContFluxdp, dContFluxdU, dContFluxdT,
                dMomFluxdp, dMomFluxdU, dMomFluxdT,
                dEnergyFluxdp, dEnergyFluxdU, dEnergyFluxdT
            );
            if (jacobianViscous_.valid())
            {
                tmp<scalarField> dContFluxdpVisc;
                tmp<vectorField> dContFluxdUVisc;
                tmp<scalarField> dContFluxdTVisc;

                tmp<vectorField> dMomFluxdpVisc;
                tmp<tensorField> dMomFluxdUVisc;
                tmp<vectorField> dMomFluxdTVisc;

                tmp<scalarField> dEnergyFluxdpVisc;
                tmp<vectorField> dEnergyFluxdUVisc;
                tmp<scalarField> dEnergyFluxdTVisc;

                jacobianViscous_->boundaryJacobian
                (
                    patchi,
                    dContFluxdpVisc, dContFluxdUVisc, dContFluxdTVisc,
                    dMomFluxdpVisc, dMomFluxdUVisc, dMomFluxdTVisc,
                    dEnergyFluxdpVisc, dEnergyFluxdUVisc, dEnergyFluxdTVisc
                );

                dContFluxdp.ref() += dContFluxdpVisc;
                dContFluxdU.ref() += dContFluxdUVisc;
                dContFluxdT.ref() += dContFluxdTVisc;

                dMomFluxdp.ref() += dMomFluxdpVisc;
                dMomFluxdU.ref() += dMomFluxdUVisc;
                dMomFluxdT.ref() += dMomFluxdTVisc;

                dEnergyFluxdp.ref() += dEnergyFluxdpVisc;
                dEnergyFluxdU.ref() += dEnergyFluxdUVisc;
                dEnergyFluxdT.ref() += dEnergyFluxdTVisc;
            }

            tmp<scalarField> rhoI = rho.boundaryField()[patchi].patchInternalField();
            tmp<vectorField> UI = U.boundaryField()[patchi].patchInternalField();
            tmp<scalarField> rhoEI = rhoE.boundaryField()[patchi].patchInternalField();
            tmp<scalarField> gammaI = gamma.boundaryField()[patchi].patchInternalField();
            tmp<scalarField> cvI = cv.boundaryField()[patchi].patchInternalField();

            tmp<scalarField> dPdRho = 0.5*(gammaI()-1)*magSqr(UI());
            tmp<vectorField> dUdRho = -UI()/rhoI();
            tmp<scalarField> dTdRho = -1.0/(cvI()*rhoI())*(rhoEI()/rhoI()-magSqr(UI()));

            tmp<vectorField> dPdRhoU = -(gammaI()-1)*UI();
            tmp<sphericalTensorField> dUdRhoU = 1.0/rhoI()*sphericalTensor::I;
            tmp<vectorField> dTdRhoU = -UI()/(cvI()*rhoI());

            tmp<scalarField> dPdRhoE = gammaI()-1;
            tmp<vectorField> dUdRhoE(new vectorField(mesh_.boundary()[patchi].size(), vector::zero));
            tmp<scalarField> dTdRhoE = 1.0/(cvI()*rhoI());

            scalarField& diagContByRho = dContByRho.diag();
            vectorField& diagContByRhoU = dContByRhoU.diag();
            scalarField& diagContByRhoE = dContByRhoE.diag();

            vectorField& diagMomByRho = dMomByRho.diag();
            tensorField& diagMomByRhoU = dMomByRhoU.diag();
            vectorField& diagMomByRhoE = dMomByRhoE.diag();

            scalarField& diagEnergyByRho = dEnergyByRho.diag();
            vectorField& diagEnergyByRhoU = dEnergyByRhoU.diag();
            scalarField& diagEnergyByRhoE = dEnergyByRhoE.diag();

            forAll(patch, bfacei)
            {

                label iIntCell = patch.faceCells()[bfacei];

                // Continuity coupled to rho
                diagContByRho[iIntCell] += dContFluxdp()[bfacei] * dPdRho()[bfacei];
                diagContByRho[iIntCell] += dContFluxdU()[bfacei] & dUdRho()[bfacei];
                diagContByRho[iIntCell] += dContFluxdT()[bfacei] * dTdRho()[bfacei];

                // Continuity coupled to rhoU
                diagContByRhoU[iIntCell] += dContFluxdp()[bfacei] * dPdRhoU()[bfacei];
                diagContByRhoU[iIntCell] += dContFluxdU()[bfacei] & dUdRhoU()[bfacei];
                diagContByRhoU[iIntCell] += dContFluxdT()[bfacei] * dTdRhoU()[bfacei];

                // Continuity coupled to rhoE
                diagContByRhoE[iIntCell] += dContFluxdp()[bfacei] * dPdRhoE()[bfacei];
                diagContByRhoE[iIntCell] += dContFluxdU()[bfacei] & dUdRhoE()[bfacei];
                diagContByRhoE[iIntCell] += dContFluxdT()[bfacei] * dTdRhoE()[bfacei];


                // Momentum coupled to rho
                diagMomByRho[iIntCell] += dMomFluxdp()[bfacei] * dPdRho()[bfacei];
                diagMomByRho[iIntCell] += dMomFluxdU()[bfacei] & dUdRho()[bfacei];
                diagMomByRho[iIntCell] += dMomFluxdT()[bfacei] * dTdRho()[bfacei];

                // Momentum coupled to rhoU
                diagMomByRhoU[iIntCell] += dMomFluxdp()[bfacei] * dPdRhoU()[bfacei];
                diagMomByRhoU[iIntCell] += dMomFluxdU()[bfacei] & dUdRhoU()[bfacei];
                diagMomByRhoU[iIntCell] += dMomFluxdT()[bfacei] * dTdRhoU()[bfacei];

                // Momentum coupled to rhoE
                diagMomByRhoE[iIntCell] += dMomFluxdp()[bfacei] * dPdRhoE()[bfacei];
                diagMomByRhoE[iIntCell] += dMomFluxdU()[bfacei] & dUdRhoE()[bfacei];
                diagMomByRhoE[iIntCell] += dMomFluxdT()[bfacei] * dTdRhoE()[bfacei];


                // Energy coupled to rho
                diagEnergyByRho[iIntCell] += dEnergyFluxdp()[bfacei] * dPdRho()[bfacei];
                diagEnergyByRho[iIntCell] += dEnergyFluxdU()[bfacei] & dUdRho()[bfacei];
                diagEnergyByRho[iIntCell] += dEnergyFluxdT()[bfacei] * dTdRho()[bfacei];

                // Energy coupled to rhoU
                diagEnergyByRhoU[iIntCell] += dEnergyFluxdp()[bfacei] * dPdRhoU()[bfacei];
                diagEnergyByRhoU[iIntCell] += dEnergyFluxdU()[bfacei] & dUdRhoU()[bfacei];
                diagEnergyByRhoU[iIntCell] += dEnergyFluxdT()[bfacei] * dTdRhoU()[bfacei];

                // Energy coupled to rhoE
                diagEnergyByRhoE[iIntCell] += dEnergyFluxdp()[bfacei] * dPdRhoE()[bfacei];
                diagEnergyByRhoE[iIntCell] += dEnergyFluxdU()[bfacei] & dUdRhoE()[bfacei];
                diagEnergyByRhoE[iIntCell] += dEnergyFluxdT()[bfacei] * dTdRhoE()[bfacei];


            }

        }

    }

}


const compressibleJacobianMatrix& jacobian::constructJacobian()
{
    jacobianInviscid_->addInviscidJacobian(jacobianMatrix_);
    if (jacobianViscous_.valid())
    {
        jacobianViscous_->addViscousJacobian(jacobianMatrix_);
    }
    addBoundaryTerms();
    return jacobianMatrix_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

jacobian::jacobian
(
    const dictionary& dict,
    const fvMesh& mesh,
    const volScalarField& rho,
    const volVectorField& rhoU,
    const volScalarField& rhoE,
    const psiThermo& thermo,
    const scalarField& ddtCoeff,
    const bool& inviscid,
    const compressible::turbulenceModel& turbulence
)
:
    mesh_(mesh),
    rho_(rho),
    rhoU_(rhoU),
    rhoE_(rhoE),
    thermo_(thermo),
    jacobianMatrix_(mesh)
{
    jacobianInviscid_ =
        jacobianInviscid::New
        (
            dict,
            mesh,
            rho,
            rhoU,
            rhoE,
            thermo,
            ddtCoeff
        );
    if (!inviscid)
    {
        jacobianViscous_ =
            jacobianViscous::New
            (
                dict,
                mesh,
                rho,
                rhoU,
                rhoE,
                thermo,
                turbulence
            );
    }
    constructJacobian();
}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

jacobian::~jacobian()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
