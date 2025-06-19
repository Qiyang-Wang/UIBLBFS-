/*---------------------------------------------------------------------------*\

    HiSA: High Speed Aerodynamic solver

    Copyright (C) 2014-2018 Johan Heyns - CSIR, South Africa
    Copyright (C) 2014-2018 Oliver Oxtoby - CSIR, South Africa
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

Class
    Foam::jacobianViscous

Description
    Base class for Jacobian schemes

SourceFiles
    jacobianViscous.C

Authors
    Oliver Oxtoby
    Johan Heyns
    Council for Scientific and Industrial Research, South Africa

\*---------------------------------------------------------------------------*/

#ifndef jacobianViscous_H
#define jacobianViscous_H

#include "jacobianMatrix.H"
#include "IOdictionary.H"
#include "volFields.H"
#include "fvCFD.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "turbulentFluidThermoModel.H"
#include "psiThermo.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class jacobianViscous Declaration
\*---------------------------------------------------------------------------*/

class jacobianViscous
{

protected:
    // Protected data

private:
    // Private Member Functions

        //- Disallow copy construct
        jacobianViscous(const jacobianViscous&);

        //- Disallow default bitwise assignment
        void operator=(const jacobianViscous&);


public:
    //- Runtime type information
    TypeName("jacobianViscous");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            jacobianViscous,
            dictionary,
            (
                const dictionary& dict,
                const fvMesh& mesh,
                const volScalarField& rho,
                const volVectorField& rhoU,
                const volScalarField& rhoE,
                const psiThermo& thermo,
                const compressible::turbulenceModel& turbulence
            ),
            (dict, mesh, rho, rhoU, rhoE, thermo, turbulence)
        );


	// Constructors

        //- Construct from components
        jacobianViscous
        (
            const word& type,
            const dictionary& dict
        );


    // Selectors

        //- Return a reference to the selected gas law
        static autoPtr<jacobianViscous> New
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const volScalarField& rho,
            const volVectorField& rhoU,
            const volScalarField& rhoE,
            const psiThermo& thermo,
            const compressible::turbulenceModel& turbulence
        );


    // Destructor

        virtual ~jacobianViscous();


    // Member Functions

        //- Add the viscous part of the Jacobian
        virtual void addViscousJacobian(compressibleJacobianMatrix& result) = 0;

        //- Return Jacobians of the boundary flux with respect to the primitive variables.
        // This is designed to be used with jacobian class which converts to Jacobians with repect to conserved
        // variables and adds to the main Jacobian.
        virtual void boundaryJacobian
        (
            label patchi,
            tmp<scalarField>& dContFluxdp, tmp<vectorField>& dContFluxdU, tmp<scalarField>& dContFluxdT,
            tmp<vectorField>& dMomFluxdp, tmp<tensorField>& dMomFluxdU, tmp<vectorField>& dMomFluxdT,
            tmp<scalarField>& dEnergyFluxdp, tmp<vectorField>& dEnergyFluxdU, tmp<scalarField>& dEnergyFluxdT
        ) = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
