/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenPDAC.
    This file was derived from the multiphaseEuler solver in OpenFOAM.

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

#include "JohnsonJacksonFrictionalStress.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace frictionalStressModels
{
defineTypeNameAndDebug(JohnsonJackson, 0);

addToRunTimeSelectionTable(frictionalStressModel, JohnsonJackson, dictionary);
}
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::kineticTheoryModels::frictionalStressModels::JohnsonJackson::
    readCoeffs(const dictionary& coeffDict)
{
    Fr_.read(coeffDict);
    eta_.read(coeffDict);
    p_.read(coeffDict);

    phi_.read(coeffDict);
    phi_ *= constant::mathematical::pi / 180.0;

    alphaDeltaMin_.read(coeffDict);

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::frictionalStressModels::JohnsonJackson::
    JohnsonJackson(const dictionary& coeffDict)
: frictionalStressModel(coeffDict),
  Fr_("Fr", dimensionSet(1, -1, -2, 0, 0), coeffDict),
  eta_("eta", dimless, coeffDict), p_("p", dimless, coeffDict),
  phi_("phi", dimless, coeffDict),
  alphaDeltaMin_("alphaDeltaMin", dimless, coeffDict)
{
    phi_ *= constant::mathematical::pi / 180.0;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::frictionalStressModels::JohnsonJackson::
    ~JohnsonJackson()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// -------------------------------------------------------------------------
// Split-invariance correction for Johnson-Jackson frictional pressure
//
// Original behaviour:
//   The frictional pressure law was evaluated from the total solid volume
//   fraction alphaTot and the same pf value was returned for every solid phase.
//
// Problem:
//   If one solid phase is artificially split into two identical phases,
//   each phase receives the same total frictional pressure contribution,
//   so the sum of the per-phase contributions is no longer consistent with
//   the single-phase case.
//
// Correction introduced here:
//   pf_s = (alpha_s / alphaTot) * pfTot(alphaTot)
//
// Rationale:
//   This preserves the original model in the single-solid-phase limit
//   (alpha_s = alphaTot), while distributing the total frictional pressure
//   consistently among identical solid phases in split-phase simulations.
//
// Important consequence:
//   frictionalPressurePrime() must be corrected consistently as the derivative
//   of the modified per-phase pressure, otherwise pPrime() becomes
//   inconsistent.
// -------------------------------------------------------------------------

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::frictionalStressModels::JohnsonJackson::
    frictionalPressure(const phaseModel& phase,
                       const phaseModel& continuousPhase,
                       const dimensionedScalar& alphaMinFriction,
                       const volScalarField& alphasMax) const
{
    // alpha    = volume fraction of the current solid phase "s"
    // alphaTot = total solid volume fraction
    //
    // Original behaviour:
    //   the Johnson-Jackson frictional pressure law was evaluated as a
    //   function of the total solids fraction alphaTot, and the same pf value
    //   was returned for every solid phase.
    //
    // This is not consistent with a phase-splitting test:
    //
    //   case A: one solid phase      -> pf = pfTot(alphaTot)
    //   case B: two identical phases -> pf1 = pf2 = pfTot(alphaTot)
    //
    // In case B, the sum of the per-phase frictional contributions no longer
    // matches the single-phase result.
    //
    // To restore consistency, the total frictional pressure is first computed
    // from alphaTot and then partitioned among the solid phases according to
    // their relative volume fraction:
    //
    //   pf_s = (alpha_s / alphaTot) * pfTot(alphaTot)
    //
    // This guarantees that:
    // - for a single solid phase, alpha_s = alphaTot and the original
    //   formulation is recovered exactly;
    // - for two identical split phases, each phase receives only its proper
    //   share of the total frictional pressure.
    const volScalarField& alpha = phase;
    const volScalarField alphaTot = 1.0 - continuousPhase;

    // Numerical safeguard to avoid division by zero in nearly empty cells.
    const dimensionedScalar alphaTotSmall("alphaTotSmall", dimless, SMALL);
    const volScalarField alphaTotEff = max(alphaTot, alphaTotSmall);

    // Original Johnson-Jackson total frictional pressure law, evaluated
    // using the total solids fraction.
    const volScalarField alphaDiff =
        max(alphaTot - alphaMinFriction, scalar(0));
    const volScalarField gap = max(alphasMax - alphaTot, alphaDeltaMin_);

    const volScalarField pfTot = Fr_ * pow(alphaDiff, eta_) / pow(gap, p_);

    // Per-species partition of the total frictional pressure.
    return (alpha / alphaTotEff) * pfTot;
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::frictionalStressModels::JohnsonJackson::
    frictionalPressurePrime(const phaseModel& phase,
                            const phaseModel& continuousPhase,
                            const dimensionedScalar& alphaMinFriction,
                            const volScalarField& alphasMax) const
{
    // alpha    = volume fraction of the current solid phase "s"
    // alphaTot = total solid volume fraction
    const volScalarField& alpha = phase;
    const volScalarField alphaTot = 1.0 - continuousPhase;

    // Numerical safeguard to avoid division by zero in alpha/alphaTot and
    // in terms involving alphaTot^2 in the denominator.
    const dimensionedScalar alphaTotSmall("alphaTotSmall", dimless, SMALL);
    const volScalarField alphaTotEff = max(alphaTot, alphaTotSmall);

    // Build the original total Johnson-Jackson frictional pressure:
    //
    //   pfTot = pfTot(alphaTot)
    //
    // This is the same law already used in the original implementation,
    // but now interpreted as the total frictional pressure of the solid
    // mixture rather than the pressure assigned independently to each phase.
    const volScalarField alphaDiff =
        max(alphaTot - alphaMinFriction, scalar(0));
    const volScalarField gap = max(alphasMax - alphaTot, alphaDeltaMin_);

    const volScalarField pfTot = Fr_ * pow(alphaDiff, eta_) / pow(gap, p_);

    // Derivative of the total Johnson-Jackson law with respect to alphaTot.
    // This is the same derivative already present in the original code.
    const volScalarField pfTotPrime =
        Fr_
        * (eta_ * pow(alphaDiff, eta_ - 1) * gap + p_ * pow(alphaDiff, eta_))
        / pow(gap, p_ + 1);

    // Main correction:
    //
    // The per-species frictional pressure is now
    //
    //   pf_s = (alpha_s / alphaTot) * pfTot(alphaTot)
    //
    // Therefore its derivative with respect to the same phase fraction alpha_s
    // is NOT obtained by simply multiplying the old derivative by
    // alpha/alphaTot.
    //
    // Using the product rule:
    //
    //   d(pf_s)/d(alpha_s)
    // = d(alpha_s/alphaTot)/d(alpha_s) * pfTot
    //   + (alpha_s/alphaTot) * d(pfTot)/d(alphaTot)
    //
    // which gives:
    //
    //   d(pf_s)/d(alpha_s)
    // = ((alphaTot - alpha_s)/alphaTot^2) * pfTot
    //   + (alpha_s/alphaTot) * pfTotPrime
    //
    // This correction is essential because frictionalPressurePrime() feeds
    // directly into pPrime(), and therefore affects both:
    //   1) the phase-pressure force in the momentum equation,
    //   2) the pressure-momentum coupling.
    //
    // If frictionalPressure() were corrected without correcting
    // frictionalPressurePrime() consistently, the model would become
    // internally inconsistent.
    return ((alphaTotEff - alpha) / sqr(alphaTotEff)) * pfTot
         + (alpha / alphaTotEff) * pfTotPrime;
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::frictionalStressModels::JohnsonJackson::nu(
    const phaseModel& phase,
    const phaseModel& continuousPhase,
    const dimensionedScalar& alphaMinFriction,
    const volScalarField& alphasMax,
    const volScalarField& pf,
    const volScalarField& rho,
    const volSymmTensorField& D) const
{
    return volScalarField::New(
        IOobject::groupName(Foam::typedName<frictionalStressModel>("nu"),
                            phase.group()),
        dimensionedScalar(dimTime, 0.5) * pf / rho * sin(phi_));
}


// ************************************************************************* //
