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

#include "characteristicFarfieldVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::characteristicFarfieldVelocityFvPatchVectorField::
characteristicFarfieldVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchVectorField(p, iF),
    TName_("T"),
    pName_("p"),
    phiName_("phi"),
    psiName_("thermo:psi"),
    UInf_(vector::zero),
    pInf_(0),
    TInf_(0),
    gamma_(0)
{
    refValue() = patchInternalField();
    refGrad() = vector::zero;
    valueFraction() = 1;
}


Foam::characteristicFarfieldVelocityFvPatchVectorField::
characteristicFarfieldVelocityFvPatchVectorField
(
    const characteristicFarfieldVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchVectorField(ptf, p, iF, mapper),
    TName_(ptf.TName_),
    pName_(ptf.pName_),
    phiName_(ptf.phiName_),
    psiName_(ptf.psiName_),
    UInf_(ptf.UInf_),
    pInf_(ptf.pInf_),
    TInf_(ptf.TInf_),
    gamma_(ptf.gamma_)
{}


Foam::characteristicFarfieldVelocityFvPatchVectorField::
characteristicFarfieldVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchVectorField(p, iF),
    TName_(dict.lookupOrDefault<word>("T", "T")),
    pName_(dict.lookupOrDefault<word>("p", "p")),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    psiName_(dict.lookupOrDefault<word>("psi", "thermo:psi")),
    UInf_(dict.lookup("UInf")),
    pInf_(readScalar(dict.lookup("pInf"))),
    TInf_(readScalar(dict.lookup("TInf"))),
    gamma_(readScalar(dict.lookup("gamma")))
{
    if (dict.found("value"))
    {
        fvPatchField<vector>::operator=
        (
            vectorField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<vector>::operator=(patchInternalField());
    }

    refValue() = *this;
    refGrad() = vector::zero;
    valueFraction() = 1;

    if (pInf_ < SMALL)
    {
        FatalIOErrorIn
        (
            "characteristicFarfieldVelocityFvPatchVectorField::"
            "characteristicFarfieldVelocityFvPatchVectorField"
            "(const fvPatch&, const vectorField&, const dictionary&)",
            dict
        )   << "    unphysical pInf specified (pInf <= 0.0)"
            << "\n    on patch " << this->patch().name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalIOError);
    }
}


Foam::characteristicFarfieldVelocityFvPatchVectorField::
characteristicFarfieldVelocityFvPatchVectorField
(
    const characteristicFarfieldVelocityFvPatchVectorField& sfspvf
)
:
    mixedFvPatchVectorField(sfspvf),
    TName_(sfspvf.TName_),
    pName_(sfspvf.pName_),
    phiName_(sfspvf.phiName_),
    psiName_(sfspvf.psiName_),
    UInf_(sfspvf.UInf_),
    pInf_(sfspvf.pInf_),
    TInf_(sfspvf.TInf_),
    gamma_(sfspvf.gamma_)
{}


Foam::characteristicFarfieldVelocityFvPatchVectorField::
characteristicFarfieldVelocityFvPatchVectorField
(
    const characteristicFarfieldVelocityFvPatchVectorField& sfspvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchVectorField(sfspvf, iF),
    TName_(sfspvf.TName_),
    pName_(sfspvf.pName_),
    phiName_(sfspvf.phiName_),
    psiName_(sfspvf.psiName_),
    UInf_(sfspvf.UInf_),
    pInf_(sfspvf.pInf_),
    TInf_(sfspvf.TInf_),
    gamma_(sfspvf.gamma_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::characteristicFarfieldVelocityFvPatchVectorField::updateCoeffs()
{
    if (!size() || updated())
    {
        return;
    }

    const fvPatchField<scalar>& pp =
        patch().lookupPatchField<volScalarField, scalar>(pName_);

    const fvsPatchField<scalar>& pphi =
        patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);

    const fvPatchField<scalar>& ppsi =
        patch().lookupPatchField<volScalarField, scalar>(psiName_);


    vectorField& Up = refValue();
    valueFraction() = 1;

    // get the near patch internal cell values
    const vectorField U(patchInternalField());
    const scalarField p(pp.patchInternalField());

    // Patch inward pointing unit vector (Same convention as Lohner)
    const vectorField np(-patch().nf());

    // Velocity pointing inwards (Same convenction as Lohner)
    //const vectorField Unp(cmptMultiply(Up,np));  // Normal
    //const vectorField Utp(pU - Unp); //OO
    //const vectorField Utp(pU - (pU&np)*np);      // tangential

    // Patch normal Mach number
    const scalarField cp(sqrt(gamma_/ppsi));
    const scalarField Mp(pphi/(patch().magSf()*cp));//OO

    // Reference values (Blazek suggests using internal values at cell centres)
    scalarField cO (sqrt(gamma_/ppsi.patchInternalField()));
    scalarField rhoO (ppsi.patchInternalField()*p);

    // Set the patch boundary condition based on the Mach number and direction
    // of the flow dictated by the boundary/free-stream pressure difference

    forAll(Up, facei)
    {
        if (Mp[facei] <= -1.0)                       // Supersonic inflow
        {
            Up[facei] = UInf_;
        }
        else if (Mp[facei] >= 1.0)                  // Supersonic outflow
        {
            valueFraction()[facei] = 0;
        }
        else if (Mp[facei] <= 0.0)                   // Subsonic inflow
        {
            // NOTE: Based on updated pp
            Up[facei] = UInf_ + np[facei]*(pInf_ - pp[facei])/(rhoO[facei]*cO[facei]);
        }
        else                                         // Subsonic outflow
        {
            // // Prescribed p
            // Up[facei] = U[facei] + np[facei]*(pInf_ - p[facei])/(rhoO[facei]*cO[facei]); // Blazek
            // //Up[facei] = U[facei] + np[facei]*(pInf_ - pp[facei])/(rhoO[facei]*cO[facei]); // Lohner

            //// Prescribed p and u_n
            //scalar pp = 0.5*(pInf_ + p[facei] + (rhoO[facei]*cO[facei]) * ((UInf_ - U[facei]) & np[facei]));
            scalar pp = p[facei]; //+ (rhoO[facei]*cO[facei]) * ((UInf_ - U[facei]) & np[facei]);
            Up[facei] = (UInf_&np[facei])*np[facei] + (U[facei]-(U[facei]&np[facei])*np[facei]) + np[facei]*(pInf_ - pp)/(rhoO[facei]*cO[facei]); //Lohner c1

        }

    }

    mixedFvPatchVectorField::updateCoeffs();
}


void Foam::characteristicFarfieldVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    #ifdef FOAM_VERSION_1712
        os.writeEntryIfDifferent<word>("T", "T", TName_);
        os.writeEntryIfDifferent<word>("p", "p", pName_);
        os.writeEntryIfDifferent<word>("phi", "phi", phiName_);
        os.writeEntryIfDifferent<word>("psi", "thermo:psi", psiName_);
    #else
        writeEntryIfDifferent<word>(os, "T", "T", TName_);
        writeEntryIfDifferent<word>(os, "p", "p", pName_);
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
        fvPatchVectorField,
        characteristicFarfieldVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
