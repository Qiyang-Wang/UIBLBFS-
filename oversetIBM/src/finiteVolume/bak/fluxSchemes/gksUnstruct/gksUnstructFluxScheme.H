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
    Foam::gksUnstructFluxScheme

Description
    AUSM+Up flux splitting scheme as described by Liou 2006, A sequal to
    AUSM, Part II: AUSM+-up for all speeds

SourceFiles
    gksUnstructFluxScheme.C

Authors
    Johan Heyns
    Oliver Oxtoby
        Council for Scientific and Industrial Research, South Africa

\*---------------------------------------------------------------------------*/

#ifndef gksUnstructFluxScheme_H
#define gksUnstructFluxScheme_H

#include "fluxScheme.H"
#include "psiThermo.H"
#include "turbulentFluidThermoModel.H"
//#include "upwind.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class gksUnstructFluxScheme Declaration
\*---------------------------------------------------------------------------*/

class gksUnstructFluxScheme
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

        const compressibleTurbulenceModel& turbulence_;
        surfaceScalarField slen_;
        

public:
    //- Runtime type information
    TypeName("GKS3DUnstruct");


	// Constructors

        //- Construct from components
        gksUnstructFluxScheme
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

        ~gksUnstructFluxScheme();


    // Member Functions

        //- Calculate mass, momentum and energy flux
        virtual void calcFlux(surfaceScalarField& phi, surfaceVectorField& phiUp, surfaceScalarField& phiEp, surfaceVectorField& Up);
        
        void viscousFlux(vector n, scalar rho, scalar P, tensor gradU, scalar nu, scalar& Grho, vector& Gm);

	void icoGKSFlux
	(
		    const scalar& rhoLeft,
		    const scalar& rhoRight,
		    const vector& ULeft,
		    const vector& URight,
		    const tensor& gradULeft,
		    const tensor& gradURight,
		    const vector& gradRhoLeft,
		    const vector& gradRhoRight,
		    const vector& leftCenter,
		    const vector& rightCenter,
		    const vector& faceCenter,
		    const vector& normalVector,
		    const vector& t1,
		    const vector& t2,
		    const scalar& nuEff,
		    const scalar& streamingLength,
		    const scalar& Cspeed,
		    scalar& rhoFlux,
		    vector& rhoUFlux,
		    scalar& rhoEFlux,
                    int ii,
		    word name
	);
 
        void calcStreamingDistant(surfaceScalarField& slen);
        scalar triangleHeight(vector c, vector f1, vector f2);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
