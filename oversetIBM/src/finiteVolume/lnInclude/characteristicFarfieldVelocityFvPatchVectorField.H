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
    Foam::characteristicFarfieldVelocityFvPatchVectorField

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
        pName      | Pressure field name         | no        | p
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
        type            characteristicFarfieldVelocity;
        UInf            500;
        pInf            1e4;
        TInf            265;
        gamma           1.4;
    }
    \endverbatim

Note
    This boundary condition is intended for farfield BC and uses a first order
    approximation on the boundary.

SourceFiles
    characteristicFarfieldVelocityFvPatchVectorField.C

Authors
    Johan Heyns
    Oliver Oxtoby
        Council for Scientific and Industrial Research, South Africa

\*---------------------------------------------------------------------------*/

#ifndef characteristicFarfieldVelocityFvPatchVectorFields_H
#define characteristicFarfieldVelocityFvPatchVectorFields_H

#include "fvPatchFields.H"
#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class characteristicFarfieldVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class characteristicFarfieldVelocityFvPatchVectorField
:
    public mixedFvPatchVectorField
{
    // Private data

        //- Name of temperature field, default = "T"
        word TName_;

        //- Name of pressure field, default = "p"
        word pName_;

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
    TypeName("characteristicFarfieldVelocity");


    // Constructors

        //- Construct from patch and internal field
        characteristicFarfieldVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        characteristicFarfieldVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given characteristicFarfieldVelocityFvPatchVectorField
        //  onto a new patch
        characteristicFarfieldVelocityFvPatchVectorField
        (
            const characteristicFarfieldVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        characteristicFarfieldVelocityFvPatchVectorField
        (
            const characteristicFarfieldVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new characteristicFarfieldVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        characteristicFarfieldVelocityFvPatchVectorField
        (
            const characteristicFarfieldVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new characteristicFarfieldVelocityFvPatchVectorField(*this, iF)
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
