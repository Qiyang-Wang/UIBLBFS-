/*---------------------------------------------------------------------------*\

    HiSA: High Speed Aerodynamic solver

    Copyright (C) 2014-2018 Oliver Oxtoby - CSIR, South Africa
    Copyright (C) 2014-2018 Johan Heyns - CSIR, South Africa
    Copyright (C) 2011-2012 OpenFOAM Foundation

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
    Foam::characteristicWallPressureFvPatchScalarField

Group
    grpWallBoundaryConditions

Description
    Characteristic pressure at walls, as per Liou, JCP 129, 364 (1996)

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            characteristicWallPressure;
    }
    \endverbatim

SeeAlso
    Foam::fixedGradientFvPatchField

SourceFiles
    characteristicWallPressureFvPatchScalarField.C

Authors
    Oliver Oxtoby
    Johan Heyns
        Council for Scientific and Industrial Research, South Africa

\*---------------------------------------------------------------------------*/

#ifndef characteristicWallPressureFvPatchScalarFields_H
#define characteristicWallPressureFvPatchScalarFields_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class characteristicWallPressureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class characteristicWallPressureFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{
    // Private data

        //- Current time index (used for updating)
        label curTimeIndex_;


public:

    //- Runtime type information
    TypeName("characteristicWallPressure");


    // Constructors

        //- Construct from patch and internal field
        characteristicWallPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        characteristicWallPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given characteristicWallPressureFvPatchScalarField onto
        //  a new patch
        characteristicWallPressureFvPatchScalarField
        (
            const characteristicWallPressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        characteristicWallPressureFvPatchScalarField
        (
            const characteristicWallPressureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new characteristicWallPressureFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        characteristicWallPressureFvPatchScalarField
        (
            const characteristicWallPressureFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new characteristicWallPressureFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        //- Update the patch pressure gradient field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
