/*---------------------------------------------------------------------------*\

    HiSA: High Speed Aerodynamic solver

    Copyright (C) 2014 Oliver Oxtoby - CSIR, South Africa
    Copyright (C) 2014 Johan Heyns - CSIR, South Africa
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
    Foam::fluxScheme

Description
    Base class for flux schemes

SourceFiles
    fluxScheme.C

Authors
    Oliver Oxtoby
    Johan Heyns
        Council for Scientific and Industrial Research, South Africa

\*---------------------------------------------------------------------------*/

#ifndef fluxScheme_incomp_H
#define fluxScheme_incomp_H

#include "IOdictionary.H"
#include "volFields.H"
#include "fvCFD.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "psiThermo.H"
#include "turbulentFluidThermoModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class fluxScheme Declaration
\*---------------------------------------------------------------------------*/

class fluxScheme_incomp

{

protected:
    // Protected data

private:
    // Private Member Functions

        //- Disallow copy construct
        fluxScheme_incomp(const fluxScheme_incomp&);

        //- Disallow default bitwise assignment
        void operator=(const fluxScheme_incomp&);


public:
    //- Runtime type information
    TypeName("fluxScheme_incomp");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            fluxScheme_incomp,
            dictionary,
            (
                const dictionary& dict,
                const psiThermo& thermo,
                const volScalarField& rho,
                const volVectorField& U,
                const volVectorField& rhoU,
                const compressibleTurbulenceModel& turbulence
            ),
            (dict, thermo, rho, U, rhoU, turbulence)
        );


	// Constructors

        //- Construct from components
        fluxScheme_incomp
        (
            const word& type,
            const dictionary& dict
        );


    // Selectors

        //- Return a reference to the selected gas law
        static autoPtr<fluxScheme_incomp> New
        (
            const dictionary& dict,
            const psiThermo& thermo,
            const volScalarField& rho,
            const volVectorField& U,
            const volVectorField& rhoU,
            const compressibleTurbulenceModel& turbulence
        );


    // Destructor

        virtual ~fluxScheme_incomp();


    // Member Functions

        //- Calculate mass, momentum and energy flux
        virtual void calcFlux(surfaceScalarField& phi, surfaceVectorField& phiUp) = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
