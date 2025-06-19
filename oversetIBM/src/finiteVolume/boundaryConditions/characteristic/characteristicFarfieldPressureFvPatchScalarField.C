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

\*---------------------------------------------------------------------------*/

#include "characteristicFarfieldPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::characteristicFarfieldPressureFvPatchScalarField::
characteristicFarfieldPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    TName_("T"),
    UName_("U"),
    phiName_("phi"),
    psiName_("thermo:psi"),
    UInf_(vector::zero),
    pInf_(0),
    TInf_(0),
    gamma_(0)
{
    refValue() = patchInternalField();
    refGrad() = 0;
    valueFraction() = 1;
}


Foam::characteristicFarfieldPressureFvPatchScalarField::
characteristicFarfieldPressureFvPatchScalarField
(
    const characteristicFarfieldPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    TName_(ptf.TName_),
    UName_(ptf.UName_),
    phiName_(ptf.phiName_),
    psiName_(ptf.psiName_),
    UInf_(ptf.UInf_),
    pInf_(ptf.pInf_),
    TInf_(ptf.TInf_),
    gamma_(ptf.gamma_)
{}


Foam::characteristicFarfieldPressureFvPatchScalarField::
characteristicFarfieldPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    TName_(dict.lookupOrDefault<word>("T", "T")),
    UName_(dict.lookupOrDefault<word>("U", "U")),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    psiName_(dict.lookupOrDefault<word>("psi", "thermo:psi")),
    UInf_(dict.lookup("UInf")),
    pInf_(readScalar(dict.lookup("pInf"))),
    TInf_(readScalar(dict.lookup("TInf"))),
    gamma_(readScalar(dict.lookup("gamma")))
{

    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<scalar>::operator=(patchInternalField());
    }

    refValue() = *this;
    refGrad() = 0;
    valueFraction() = 1;

    if (pInf_ < SMALL)
    {
        FatalIOErrorIn
        (
            "characteristicFarfieldPressureFvPatchScalarField::"
            "characteristicFarfieldPressureFvPatchScalarField"
            "(const fvPatch&, const vectorField&, const dictionary&)",
            dict
        )   << "    unphysical pInf specified (pInf <= 0.0)"
            << "\n    on patch " << this->patch().name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalIOError);
    }
}


Foam::characteristicFarfieldPressureFvPatchScalarField::
characteristicFarfieldPressureFvPatchScalarField
(
    const characteristicFarfieldPressureFvPatchScalarField& sfspvf
)
:
    mixedFvPatchScalarField(sfspvf),
    TName_(sfspvf.TName_),
    UName_(sfspvf.UName_),
    phiName_(sfspvf.phiName_),
    psiName_(sfspvf.psiName_),
    UInf_(sfspvf.UInf_),
    pInf_(sfspvf.pInf_),
    TInf_(sfspvf.TInf_),
    gamma_(sfspvf.gamma_)
{}


Foam::characteristicFarfieldPressureFvPatchScalarField::
characteristicFarfieldPressureFvPatchScalarField
(
    const characteristicFarfieldPressureFvPatchScalarField& sfspvf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(sfspvf, iF),
    TName_(sfspvf.TName_),
    UName_(sfspvf.UName_),
    phiName_(sfspvf.phiName_),
    psiName_(sfspvf.psiName_),
    UInf_(sfspvf.UInf_),
    pInf_(sfspvf.pInf_),
    TInf_(sfspvf.TInf_),
    gamma_(sfspvf.gamma_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::characteristicFarfieldPressureFvPatchScalarField::updateCoeffs()
{
    if (!size() || updated())
    {
        return;
    }

    const fvPatchField<vector>& pU =
        patch().lookupPatchField<volVectorField, vector>(UName_);

    const fvsPatchField<scalar>& pphi =
        patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);

    const fvPatchField<scalar>& ppsi =
        patch().lookupPatchField<volScalarField, scalar>(psiName_);


    scalarField& pp = refValue();
    valueFraction() = 1;

    // get the near patch internal cell values
    const scalarField p(patchInternalField());
    const vectorField U(pU.patchInternalField());

    // Patch inward pointing face unit vector (Same convenction as Lohner)
    const vectorField np(-patch().nf());

    // Velocity pointing inwards (Same convenction as Lohner)
    //const vectorField Unp(cmptMultiply(Up,np));  // Normal
    //const vectorField Utp(pU - Unp); //OO
    //const vectorField Utp(pU - (pU&np)*np);      // tangential

    // Patch normal Mach number
    const scalarField cp(sqrt(gamma_/ppsi));
    const scalarField Mp(pphi/(patch().magSf()*cp));

    // Reference values (Blazek suggests using internal values at cell centres)
    const scalarField cO (sqrt(gamma_/ppsi.patchInternalField()));
    const scalarField rhoO (ppsi.patchInternalField()*p);

    // Set the patch boundary condition based on the Mach number and direction
    // of the flow dictated by the boundary/free-stream pressure difference

    forAll(pp, facei)
    {
        // Info << Mp[facei] << endl;
        if (Mp[facei] <= -1.0)                       // Supersonic inflow
        {
            pp[facei] = pInf_;
        }
        else if (Mp[facei] >= 1.0)                   // Supersonic outflow
        {
            valueFraction()[facei] = 0;
        }
        else if (Mp[facei] <= 0.0)                   // Subsonic inflow
        {
            // NOTE:
            // Lohner uses vertex centred and therefore * refers to the vertex value on
            // the BP, whereas Blazek uses cell-centred and therefore the internal cell
            // value.
            // pp[facei] = 0.5*(pInf_ + p[facei] + (rhoO[facei]*cO[facei]) * ((UInf_ - U[facei]) & np[facei])); // Blazek (cell-centred)
            valueFraction()[facei] = 0.5;
            pp[facei] = pInf_ + (rhoO[facei]*cO[facei]) * ((UInf_ - U[facei]) & np[facei]); // Blazek (cell-centred)
            // pp = 0.5*(pInf_ + p[facei] + (rhoO[facei]*cO[facei]) * ((UInf_ - pU[facei]) & np[facei]));    // Lohner (vertex-centred)
        }
        else                                         // Subsonic outflow
        {
            // // Blazek / Lohner c2 (prescribed p)
            // pp[facei] = pInf_;

            // // Lohner c1 - prescribed p and u_n
            // Provide better reflection properties when compared to c1 (see testBc case)
            valueFraction()[facei] = 0.5;
            pp[facei] = pInf_ + (rhoO[facei]*cO[facei]) * ((UInf_ - U[facei]) & np[facei]); // Blazek (cell-centred)

        }

    }

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::characteristicFarfieldPressureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    #ifdef FOAM_VERSION_1712
        os.writeEntryIfDifferent<word>("T", "T", TName_);
        os.writeEntryIfDifferent<word>("U", "U", UName_);
        os.writeEntryIfDifferent<word>("phi", "phi", phiName_);
        os.writeEntryIfDifferent<word>("psi", "thermo:psi", psiName_);
    #else
        writeEntryIfDifferent<word>(os, "T", "T", TName_);
        writeEntryIfDifferent<word>(os, "U", "U", UName_);
        writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
        writeEntryIfDifferent<word>(os, "psi", "thermo:psi", psiName_);
    #endif
    os.writeKeyword("UInf") << UInf_ << token::END_STATEMENT << nl;
    os.writeKeyword("pInf") << pInf_ << token::END_STATEMENT << nl;
    os.writeKeyword("TInf") << TInf_ << token::END_STATEMENT << nl;
    os.writeKeyword("gamma") << gamma_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        characteristicFarfieldPressureFvPatchScalarField
    );
}

// ************************************************************************* //
