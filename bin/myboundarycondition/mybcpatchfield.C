/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "mybcpatchfield.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "fvPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mybcpatchfield::
mybcpatchfield
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    // refHeight_(0.0),
    refVelocity_(0.0)
{}


Foam::mybcpatchfield::
mybcpatchfield
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict, IOobjectOption::NO_READ),
    // refHeight_(0.0),
    refVelocity_(0.0)
{
    // dict.lookup("refHeight")>> refHeight_;
    dict.lookup("refVelocity")>> refVelocity_;

    updateCoeffs();
}


Foam::mybcpatchfield::
mybcpatchfield
(
    const mybcpatchfield& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(p, iF),
    // refHeight_(ptf.refHeight_),
    refVelocity_(ptf.refVelocity_)
{
    // Note: refValue_ will have default value of 0 for unmapped faces. This
    // can temporarily happen during e.g. redistributePar. We should not
    // access ptf.patch() instead since redistributePar has destroyed this
    // at the time of mapping.

    updateCoeffs();
}


Foam::mybcpatchfield::
mybcpatchfield
(
    const mybcpatchfield& ptf
)
:
    fixedValueFvPatchVectorField(ptf),
    // refHeight_(ptf.refHeight_),
    refVelocity_(ptf.refVelocity_)
{}


Foam::mybcpatchfield::
mybcpatchfield
(
    const mybcpatchfield& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    // refHeight_(ptf.refHeight_),
    refVelocity_(ptf.refVelocity_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mybcpatchfield::autoMap
(
    const fvPatchFieldMapper& mapper
)
{
    fixedValueFvPatchVectorField::autoMap(mapper);
}


void Foam::mybcpatchfield::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);

    // const mybcpatchfield& tiptf =
    //     refCast<const mybcpatchfield>(ptf);

}


void Foam::mybcpatchfield::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    vectorField tvalues= patch().nf();

    forAll(tvalues, faceI)
    {
        tvalues[faceI]=refVelocity_*2.0*tvalues[faceI];
    }


    fvPatchVectorField::operator=(tvalues);
    fvPatchVectorField::updateCoeffs();
}


void Foam::mybcpatchfield::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    // os.writeKeyword("refHeight") <<refHeight_ << token::END_STATEMENT<<nl;
    os.writeKeyword("refVelocity") <<refVelocity_ << token::END_STATEMENT<<nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        mybcpatchfield
    );
}

// ************************************************************************* //
