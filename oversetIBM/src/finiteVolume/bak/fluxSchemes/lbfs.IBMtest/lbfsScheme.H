/*---------------------------------------------------------------------------*\

    HiSA: High Speed Aerodynamic solver

    Copyright (C) 2014-2016 Johan Heyns - CSIR, South Africa
    Copyright (C) 2014-2016 Oliver Oxtoby - CSIR, South Africa
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
    Foam::lbfsScheme

Description
    AUSM+Up flux splitting scheme as described by Liou 2006, A sequal to
    AUSM, Part II: AUSM+-up for all speeds

SourceFiles
    lbfsScheme.C

Authors
    Johan Heyns
    Oliver Oxtoby
        Council for Scientific and Industrial Research, South Africa

\*---------------------------------------------------------------------------*/

#ifndef lbfsScheme_H
#define lbfsScheme_H

#include "fluxScheme.H"
#include "psiThermo.H"
//#include "upwind.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class lbfsScheme Declaration
\*---------------------------------------------------------------------------*/

class lbfsScheme
:
    public fluxScheme
{

protected:
    // Protected data
        const fvMesh& mesh_;
        const psiThermo& thermo_;
        const volScalarField& rho_;
        const volVectorField& U_;
        const volVectorField& rhoU_;
        const volScalarField& rhoE_;
        const dictionary& dict_;

public:
    //- Runtime type information
    TypeName("LBFS");


	// Constructors

        //- Construct from components
        lbfsScheme
        (
            const dictionary& dict,
            const psiThermo& thermo,
            const volScalarField& rho,
            const volVectorField& U,
            const volVectorField& rhoU,
            const volScalarField& rhoE,
            const compressibleTurbulenceModel& turbulence
        );


    // Destructor

        ~lbfsScheme();


    // Member Functions

        //- Calculate mass, momentum and energy flux
        virtual void calcFlux(surfaceScalarField& phi, surfaceVectorField& phiUp, surfaceScalarField& phiEp, surfaceVectorField& Up);

private:
        void lbfs_D1Q4
        (
            scalar slx,
            scalar sly,
            scalar slz,
            scalar rho_L,
            scalar rho_R,
            vector U_L,
            vector U_R,
            scalar p_L,
            scalar p_R,
            scalar& F_1,
            vector& F_24,
            scalar& F_5
        );

        void feq_d1q4
        (
           int side,
           scalar rou,
           scalar uu,
           scalar p,
           scalar gamma,
           scalar node[4],
           scalar d1,
           scalar d2,
           scalar &ep 
        );

        template<class Type>
        void recon
        (
           const GeometricField<Type, fvPatchField, volMesh>& vf,
           tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>& tsf_L,
           tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>& tsf_R
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
