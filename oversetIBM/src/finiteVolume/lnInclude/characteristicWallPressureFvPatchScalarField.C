/*---------------------------------------------------------------------------*\

    HiSA: High Speed Aerodynamic solver

    Copyright (C) 2014-2018 Oliver Oxtoby - CSIR, South Africa
    Copyright (C) 2014-2018 Johan Heyns - CSIR, South Africa
    Copyright (C) 2011-2013 OpenFOAM Foundation

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

#include "characteristicWallPressureFvPatchScalarField.H"
#include "fvCFD.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fluidThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::characteristicWallPressureFvPatchScalarField::characteristicWallPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF)
{}


Foam::characteristicWallPressureFvPatchScalarField::characteristicWallPressureFvPatchScalarField
(
    const characteristicWallPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(p, iF)
{
    patchType() = ptf.patchType();

    // Evaluate the value field from the gradient if the internal field is valid
    if (notNull(iF) && iF.size())
    {
        scalarField::operator=
        (
            //patchInternalField() + gradient()/patch().deltaCoeffs()
            // ***HGW Hack to avoid the construction of mesh.deltaCoeffs
            // which fails for AMI patches for some mapping operations
            patchInternalField() + gradient()*(patch().nf() & patch().delta())
        );
    }
    else
    {
        // Enforce mapping of values so we have a valid starting value. This
        // constructor is used when reconstructing fields
        this->map(ptf, mapper);
    }
}


Foam::characteristicWallPressureFvPatchScalarField::characteristicWallPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF)
{
    if (dict.found("value") && dict.found("gradient"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
        gradient() = scalarField("gradient", dict, p.size());
    }
    else
    {
        fvPatchField<scalar>::operator=(patchInternalField());
        gradient() = 0.0;
    }
}


Foam::characteristicWallPressureFvPatchScalarField::characteristicWallPressureFvPatchScalarField
(
    const characteristicWallPressureFvPatchScalarField& wbppsf
)
:
    fixedGradientFvPatchScalarField(wbppsf)
{}


Foam::characteristicWallPressureFvPatchScalarField::characteristicWallPressureFvPatchScalarField
(
    const characteristicWallPressureFvPatchScalarField& wbppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(wbppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::characteristicWallPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    vectorField Uif =
        patch().lookupPatchField<volVectorField, vector>("U").patchInternalField();
    scalarField psiif =
        patch().lookupPatchField<volScalarField, scalar>("thermo:psi").patchInternalField();

    const fvMesh& mesh = patch().boundaryMesh().mesh();

    const fluidThermo& thermo =
        mesh.lookupObject<fluidThermo>("thermophysicalProperties");
    tmp< volScalarField > gamma = thermo.gamma();
    scalarField gammaif =
        gamma->boundaryField()[patch().index()].patchInternalField();

    scalarField pif =
        patch().lookupPatchField<volScalarField, scalar>("p").patchInternalField();

    // Reference values (Blazek suggests using internal values at cell centres)
    const scalarField cif (sqrt(gammaif/psiif));
    const scalarField rhoif (psiif*pif);

    gradient() = rhoif*cif*(Uif&patch().nf())*patch().deltaCoeffs();

    fixedGradientFvPatchScalarField::updateCoeffs();
}


void Foam::characteristicWallPressureFvPatchScalarField::write(Ostream& os) const
{
    fixedGradientFvPatchScalarField::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        characteristicWallPressureFvPatchScalarField
    );
}


// ************************************************************************* //
