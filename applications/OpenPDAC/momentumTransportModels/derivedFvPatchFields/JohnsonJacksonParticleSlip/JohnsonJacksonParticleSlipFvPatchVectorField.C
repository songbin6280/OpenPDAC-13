/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2024 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "JohnsonJacksonParticleSlipFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "kineticTheoryModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
makePatchTypeField(fvPatchVectorField,
                   JohnsonJacksonParticleSlipFvPatchVectorField);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::JohnsonJacksonParticleSlipFvPatchVectorField::
    JohnsonJacksonParticleSlipFvPatchVectorField(
        const fvPatch& p,
        const DimensionedField<vector, volMesh>& iF,
        const dictionary& dict)
: partialSlipFvPatchVectorField(p, iF),
  specularityCoefficient_(
      dict.lookup<scalar>("specularityCoefficient", unitFraction))
{
    if (specularityCoefficient_ < 0 || specularityCoefficient_ > 1)
    {
        FatalErrorInFunction
            << "The specularity coefficient has to be between 0 and 1"
            << abort(FatalError);
    }

    fvPatchVectorField::operator=(
        vectorField("value", iF.dimensions(), dict, p.size()));
}


Foam::JohnsonJacksonParticleSlipFvPatchVectorField::
    JohnsonJacksonParticleSlipFvPatchVectorField(
        const JohnsonJacksonParticleSlipFvPatchVectorField& ptf,
        const fvPatch& p,
        const DimensionedField<vector, volMesh>& iF,
        const fieldMapper& mapper)
: partialSlipFvPatchVectorField(ptf, p, iF, mapper),
  specularityCoefficient_(ptf.specularityCoefficient_)
{
}


Foam::JohnsonJacksonParticleSlipFvPatchVectorField::
    JohnsonJacksonParticleSlipFvPatchVectorField(
        const JohnsonJacksonParticleSlipFvPatchVectorField& ptf,
        const DimensionedField<vector, volMesh>& iF)
: partialSlipFvPatchVectorField(ptf, iF),
  specularityCoefficient_(ptf.specularityCoefficient_)
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::JohnsonJacksonParticleSlipFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // lookup the fluid model and the phase
    const phaseSystem& fluid =
        db().lookupObject<phaseSystem>(phaseSystem::propertiesName);

    const phaseModel& phase(fluid.phases()[internalField().group()]);

    // lookup the multi-solid crowding term
    // Sigma_s = sum_j(alpha_j*g0_sj)
    const fvPatchScalarField& sumAlphaGs0(
        patch().lookupPatchField<volScalarField, scalar>(IOobject::groupName(
            Foam::typedName<RASModels::kineticTheoryModel>("sumAlphaGs0"),
            phase.name())));

    // compute the total maximum packing field for the solid mixture
    const tmp<volScalarField> talfasMax(fluid.alfasMax());

    // extract the patch field of the total maximum packing
    const fvPatchScalarField& alfasMax(refCast<const fvPatchScalarField>(
        talfasMax().boundaryField()[patch().index()]));

    // lookup the granular viscosity of the current solid phase
    const scalarField nu(patch().lookupPatchField<volScalarField, scalar>(
        IOobject::groupName("nut", phase.name())));

    // lookup the granular temperature of the current solid phase
    word ThetaName(IOobject::groupName("Theta", phase.name()));

    const fvPatchScalarField& Theta(
        db().foundObject<volScalarField>(ThetaName)
            ? patch().lookupPatchField<volScalarField, scalar>(ThetaName)
            : patch().lookupPatchField<volScalarField, scalar>(
                  phase.volScalarField::name()));

    // calculate the slip value fraction using the multi-solid correction:
    // alpha_s*g0_s -> Sigma_s = sum_j(alpha_j*g0_sj)
    scalarField c(constant::mathematical::pi * sumAlphaGs0
                  * specularityCoefficient_ * sqrt(3 * Theta)
                  / max(6 * nu * alfasMax, small));

    this->valueFraction() = c / (c + patch().deltaCoeffs());

    partialSlipFvPatchVectorField::updateCoeffs();
}


void Foam::JohnsonJacksonParticleSlipFvPatchVectorField::write(
    Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry(os, "specularityCoefficient", specularityCoefficient_);
    writeEntry(os, "value", *this);
}


// ************************************************************************* //
