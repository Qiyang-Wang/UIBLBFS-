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

#include "characteristicFarfieldTemperatureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::characteristicFarfieldTemperatureFvPatchScalarField::
characteristicFarfieldTemperatureFvPatchScalarField
(
    const fvPatch& t,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(t, iF),
    pName_("p"),
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


Foam::characteristicFarfieldTemperatureFvPatchScalarField::
characteristicFarfieldTemperatureFvPatchScalarField
(
    const characteristicFarfieldTemperatureFvPatchScalarField& ptf,
    const fvPatch& t,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, t, iF, mapper),
    pName_(ptf.pName_),
    UName_(ptf.UName_),
    phiName_(ptf.phiName_),
    psiName_(ptf.psiName_),
    UInf_(ptf.UInf_),
    pInf_(ptf.pInf_),
    TInf_(ptf.TInf_),
    gamma_(ptf.gamma_)
{}


Foam::characteristicFarfieldTemperatureFvPatchScalarField::
characteristicFarfieldTemperatureFvPatchScalarField
(
    const fvPatch& t,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(t, iF),
    pName_(dict.lookupOrDefault<word>("p", "p")),
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
            scalarField("value", dict, t.size())
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
            "characteristicFarfieldTemperatureFvPatchScalarField::"
            "characteristicFarfieldTemperatureFvPatchScalarField"
            "(const fvPatch&, const vectorField&, const dictionary&)",
            dict
        )   << "    unphysical pInf specified (pInf <= 0.0)"
            << "\n    on patch " << this->patch().name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalIOError);
    }
}


Foam::characteristicFarfieldTemperatureFvPatchScalarField::
characteristicFarfieldTemperatureFvPatchScalarField
(
    const characteristicFarfieldTemperatureFvPatchScalarField& sfspvf
)
:
    mixedFvPatchScalarField(sfspvf),
    pName_(sfspvf.pName_),
    UName_(sfspvf.UName_),
    phiName_(sfspvf.phiName_),
    psiName_(sfspvf.psiName_),
    UInf_(sfspvf.UInf_),
    pInf_(sfspvf.pInf_),
    TInf_(sfspvf.TInf_),
    gamma_(sfspvf.gamma_)
{}


Foam::characteristicFarfieldTemperatureFvPatchScalarField::
characteristicFarfieldTemperatureFvPatchScalarField
(
    const characteristicFarfieldTemperatureFvPatchScalarField& sfspvf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(sfspvf, iF),
    pName_(sfspvf.pName_),
    UName_(sfspvf.UName_),
    phiName_(sfspvf.phiName_),
    psiName_(sfspvf.psiName_),
    UInf_(sfspvf.UInf_),
    pInf_(sfspvf.pInf_),
    TInf_(sfspvf.TInf_),
    gamma_(sfspvf.gamma_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::characteristicFarfieldTemperatureFvPatchScalarField::updateCoeffs()
{
    if (!size() || updated())
    {
        return;
    }

    const fvPatchField<scalar>& pp =
        patch().lookupPatchField<volScalarField, scalar>(pName_);

    const fvPatchField<vector>& pU =
        patch().lookupPatchField<volVectorField, vector>(UName_);

    const fvsPatchField<scalar>& pphi =
        patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);

    const fvPatchField<scalar>& ppsi =
        patch().lookupPatchField<volScalarField, scalar>(psiName_);


    scalarField& pT = refValue();
    valueFraction() = 1;

    // get the near patch internal cell values
    const scalarField T(patchInternalField());
    const scalarField p(pp.patchInternalField());
    const vectorField U(pU.patchInternalField());

    // Patch inward pointing unit vector (Same convention as Lohner)
    const vectorField np(-patch().nf());

    // Velocity pointing inwards (Same convenction as Lohner)
    //const vectorField Unp(cmptMultiply(Up,np));  // Normal
    //const vectorField Utp(pU - Unp); //OO
    //const vectorField Utp(pU - (pU&np)*np);      // tangential

    // Need R of the free-stream flow. Assume R is independent of location
    // along patch so use face 0
    scalar R = 1.0/(ppsi.patchInternalField()()[0]*patchInternalField()()[0]);
    scalar rhoInf = pInf_/(R*TInf_);

    // Patch normal Mach number
    const scalarField cp(sqrt(gamma_/ppsi));
    const scalarField Mp(pphi/(patch().magSf()*cp));

    // Reference values (Blazek suggests using internal values at cell centres)
    scalarField cO (sqrt(gamma_/ppsi.patchInternalField()));
    const scalarField rhoO (ppsi.patchInternalField()*p);

    // Set the patch boundary condition based on the Mach number and direction
    // of the flow dictated by the boundary/free-stream pressure difference

    forAll(pT, facei)
    {
        if (Mp[facei] <= -1.0)                       // Supersonic inflow
        {
            pT[facei] = TInf_;
        }
        else if (Mp[facei] >= 1.0)                   // Supersonic outflow
        {
            valueFraction()[facei] = 0;
        }
        else if (Mp[facei] <= 0.0)                   // Subsonic inflow
        {
            // NOTE: Based on updated pp
            scalar prho = rhoInf + (pp[facei] - pInf_)/(cO[facei]*cO[facei]);
            pT[facei] = pp[facei]/(R*prho);
        }
        else                                         // Subsonic outflow
        {

            // // Blazek / Lohner c2 (prescribed p)
            // pT[facei] = (pInf_*T[facei]*cO[facei]*cO[facei])/
            //             (p[facei]*cO[facei]*cO[facei] + R*T[facei]*(pInf_ - p[facei]));

            // Lohner c1 - prescribed p and u_n
            pT[facei] = pp[facei]*T[facei]*cO[facei]*cO[facei]/(p[facei]*cO[facei]*cO[facei] + R*T[facei]*(pp[facei]-p[facei]));
        }
    }

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::characteristicFarfieldTemperatureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    #ifdef FOAM_VERSION_1712
        os.writeEntryIfDifferent<word>("p", "p", pName_);
        os.writeEntryIfDifferent<word>("U", "U", UName_);
        os.writeEntryIfDifferent<word>("phi", "phi", phiName_);
        os.writeEntryIfDifferent<word>("psi", "thermo:psi", psiName_);
    #else
        writeEntryIfDifferent<word>(os, "p", "p", pName_);
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
        characteristicFarfieldTemperatureFvPatchScalarField
    );
}

// ************************************************************************* //
