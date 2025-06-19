/*---------------------------------------------------------------------------*\

    HiSA: High Speed Aerodynamic solver

    Copyright (C) 2014-2018 Johan Heyns - CSIR, South Africa
    Copyright (C) 2014-2018 Oliver Oxtoby - CSIR, South Africa
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

const volScalarField::Boundary& pbf = p.boundaryField();
wordList rhoBoundaryTypes = pbf.types();

forAll(rhoBoundaryTypes, patchi)
{
    if (rhoBoundaryTypes[patchi] == "waveTransmissive")
    {
        rhoBoundaryTypes[patchi] = zeroGradientFvPatchScalarField::typeName;
    }
    else if (rhoBoundaryTypes[patchi] == "characteristicWallPressure")
    {
        rhoBoundaryTypes[patchi] = zeroGradientFvPatchScalarField::typeName;
    }
    // For characteristic-based calculations fixedRho is used which computes rho
    // directly from p & T
    else if (pbf[patchi].fixesValue() || (rhoBoundaryTypes[patchi] == "characteristicFarfieldPressure") )
    {
        rhoBoundaryTypes[patchi] = fixedRhoFvPatchScalarField::typeName; // NOTE: Based on current values (rho = p * psi)
    }
}
