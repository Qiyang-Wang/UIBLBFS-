/*---------------------------------------------------------------------------*\

    HiSA: High Speed Aerodynamic solver

    Copyright (C) 2014-2018 Johan Heyns - CSIR, South Africa
    Copyright (C) 2014-2018 Oliver Oxtoby - CSIR, South Africa
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
    Foam::characteristicFarfieldPressureFvPatchScalarField

Group
    grpInletBoundaryConditions grpOutletBoundaryConditions

Description
    This boundary condition provides a characteristic farfield condition. For
    more detail please refer to Blazek (J. Blazek, Computational Fluid Dynamics:
    Principles and Applications, Elsevier, 2001) and Lohner (R. Löhner,
    An Introduction to Applied Computational Fluid Dynamics Techniques, John
    Wiley and sons, 2001).

    \heading Patch usage

    \table
        Property   | Description                 | Required  | Default value
        TName      | Temperature field name      | no        | T
        UName      | Velocity field name         | no        | U
        phiName    | Flux field name             | no        | phi
        psiName    | Compressibility field name  | no        | thermo:psi
        UInf       | free-stream velocity        | yes       |
        pInf       | free-stream pressure        | yes       |
        TInf       | free-stream temperature     | yes       |
        gamma      | heat capacity ratio (cp/Cv) | yes       |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            characteristicFarfieldPressure;
        UInf            500;
        pInf            1e4;
        TInf            265;
        gamma           1.4;
    }
    \endverbatim

Note
    This boundary condition is intended for a farfield BC.

SourceFiles
    characteristicFarfieldPressureFvPatchScalarField.C

Authors
    Johan Heyns
    Oliver Oxtoby
        Council for Scientific and Industrial Research, South Africa

\*---------------------------------------------------------------------------*/

#ifndef characteristicFarfieldPressureFvPatchScalarFields_H
#define characteristicFarfieldPressureFvPatchScalarFields_H

#include "fvPatchFields.H"
#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class characteristicFarfieldPressureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class characteristicFarfieldPressureFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private data

        //- Name of temperature field, default = "T"
        word TName_;

        //- Name of velocity field, default = "U"
        word UName_;

        //- Name of flux field, default = "phi"
        word phiName_;

        //- Name of compressibility field field, default = "thermo:psi"
        word psiName_;

        //- Velocity of the free stream
        vector UInf_;

        //- Pressure of the free stream
        scalar pInf_;

        //- Temperature of the free stream
        scalar TInf_;

        //- Heat capacity ratio
        scalar gamma_;


public:

    //- Runtime type information
    TypeName("characteristicFarfieldPressure");


    // Constructors

        //- Construct from patch and internal field
        characteristicFarfieldPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        characteristicFarfieldPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given characteristicFarfieldPressureFvPatchScalarField
        //  onto a new patch
        characteristicFarfieldPressureFvPatchScalarField
        (
            const characteristicFarfieldPressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        characteristicFarfieldPressureFvPatchScalarField
        (
            const characteristicFarfieldPressureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new characteristicFarfieldPressureFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        characteristicFarfieldPressureFvPatchScalarField
        (
            const characteristicFarfieldPressureFvPatchScalarField&,
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
                new characteristicFarfieldPressureFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the velocity at infinity
            const vector& UInf() const
            {
                return UInf_;
            }

            //- Return reference to the velocity at infinity to allow adjustment
            vector& UInf()
            {
                return UInf_;
            }

            //- Return the pressure at infinity
            scalar pInf() const
            {
                return pInf_;
            }

            //- Return reference to the pressure at infinity to allow adjustment
            scalar& pInf()
            {
                return pInf_;
            }

            //- Return the temperature at infinity
            scalar TInf() const
            {
                return TInf_;
            }

            //- Return reference to the temperature at infinity
            //  to allow adjustment
            scalar& TInf()
            {
                return TInf_;
            }


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
