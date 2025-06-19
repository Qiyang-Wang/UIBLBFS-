/*---------------------------------------------------------------------------*\

    HiSA: High Speed Aerodynamic solver

    Copyright (C) 2015-2017 Johan Heyns - CSIR, South Africa
    Copyright (C) 2015-2017 Oliver Oxtoby - CSIR, South Africa

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

Class
    laxFriedrich

Description
    Lax-Friedrich approximation of Inviscid LHS Jacobian
    as per Luo, Baum, Lohner, JCP 146 664-690 (1998).

Authors
    Johan Heyns
    Oliver Oxtoby
        Council for Scientific and Industrial Research, South Africa

\*---------------------------------------------------------------------------*/

#ifndef laxFriedrich_H
#define laxFriedrich_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "jacobianInviscid.H"
#include "jacobianMatrix.H"
#include "fvCFD.H"
#include "turbulenceModel.H"
#include "psiThermo.H"
#include "fvjMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class laxFriedrich Declaration
\*---------------------------------------------------------------------------*/

class laxFriedrich
:
    public jacobianInviscid
{

    protected:

        const fvMesh& mesh_;
        const volScalarField& rho_;
        const volVectorField& rhoU_;
        const volScalarField& rhoE_;
        const psiThermo& thermo_;
        tmp< surfaceScalarField > lambdaConv_;
        tmp< volScalarField > gamma_;
        const scalarField& ddtCoeff_;  // 1/dt or similar
        autoPtr< tmp<surfaceScalarField> > meshPhi_;


        //- Private member functions

            void addFluxTerms(compressibleJacobianMatrix& jacobian);
            void addTemporalTerms(compressibleJacobianMatrix& jacobian);


    public:
        //- Runtime type information
        TypeName("LaxFriedrich");

        laxFriedrich
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const volScalarField& rho,
            const volVectorField& rhoU,
            const volScalarField& rhoE,
            const psiThermo& thermo,
            const scalarField& ddtCoeff
        );

        virtual void addInviscidJacobian(compressibleJacobianMatrix& jacobian);

        virtual void boundaryJacobian
        (
            label patchi,
            tmp<scalarField>& dContFluxdp, tmp<vectorField>& dContFluxdU, tmp<scalarField>& dContFluxdT,
            tmp<vectorField>& dMomFluxdp, tmp<tensorField>& dMomFluxdU, tmp<vectorField>& dMomFluxdT,
            tmp<scalarField>& dEnergyFluxdp, tmp<vectorField>& dEnergyFluxdU, tmp<scalarField>& dEnergyFluxdT
        );

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
