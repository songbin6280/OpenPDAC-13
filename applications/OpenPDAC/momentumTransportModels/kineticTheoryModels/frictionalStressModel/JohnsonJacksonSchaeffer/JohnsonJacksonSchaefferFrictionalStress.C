/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2024 OpenFOAM Foundation
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

#include "JohnsonJacksonSchaefferFrictionalStress.H"
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
defineTypeNameAndDebug(JohnsonJacksonSchaeffer, 0);

addToRunTimeSelectionTable(frictionalStressModel,
                           JohnsonJacksonSchaeffer,
                           dictionary);
}
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::kineticTheoryModels::frictionalStressModels::
    JohnsonJacksonSchaeffer::readCoeffs(const dictionary& coeffDict)
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

Foam::kineticTheoryModels::frictionalStressModels::JohnsonJacksonSchaeffer::
    JohnsonJacksonSchaeffer(const dictionary& coeffDict)
: frictionalStressModel(coeffDict),
  Fr_("Fr", dimensionSet(1, -1, -2, 0, 0), coeffDict),
  eta_("eta", dimless, coeffDict), p_("p", dimless, coeffDict),
  phi_("phi", dimless, coeffDict),
  alphaDeltaMin_("alphaDeltaMin", dimless, coeffDict)
{
    phi_ *= constant::mathematical::pi / 180.0;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::frictionalStressModels::JohnsonJacksonSchaeffer::
    ~JohnsonJacksonSchaeffer()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// -------------------------------------------------------------------------
// Split-invariance correction for frictional pressure
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
Foam::kineticTheoryModels::frictionalStressModels::JohnsonJacksonSchaeffer::
    frictionalPressure(const phaseModel& phase,
                       const phaseModel& continuousPhase,
                       const dimensionedScalar& alphaMinFriction,
                       const volScalarField& alphasMax) const
{
    // alpha    = volume fraction of the current solid phase "s"
    // alphaTot = total solid volume fraction
    //
    // In the original implementation, the frictional pressure law
    // was evaluated directly as a function of alphaTot, and the same pf value
    // was returned for every solid phase.
    //
    // This creates an inconsistency when a single solid phase is artificially
    // split into two identical phases:
    //
    //   case A: one solid phase      -> pf = pfTot(alphaTot)
    //   case B: two identical phases -> pf1 = pf2 = pfTot(alphaTot)
    //
    // In that case, the sum of the frictional contributions of the two species
    // does not reconstruct the contribution of the single-phase case.
    //
    // To avoid this problem, a per-species partition is introduced:
    //
    //   pf_s = (alpha_s / alphaTot) * pfTot(alphaTot)
    //
    // With this choice:
    // - if only one solid phase exists, alpha_s = alphaTot and the original
    //   formulation is exactly recovered;
    // - if one phase is split into two identical phases, each phase receives
    //   a fraction of the total frictional pressure proportional to its own
    //   volume fraction.
    const volScalarField& alpha = phase;
    const volScalarField alphaTot = 1.0 - continuousPhase;

    // Numerical safeguard to avoid division by zero in nearly empty cells.
    // This regularization does not affect the normal solid-present regime,
    // but prevents singularities in the alpha/alphaTot ratio.
    const dimensionedScalar alphaTotSmall("alphaTotSmall", dimless, SMALL);
    const volScalarField alphaTotEff = max(alphaTot, alphaTotSmall);

    // "Total" Johnson-Jackson/Schaeffer frictional pressure law:
    // this is the original model, evaluated as a function of the total
    // solid volume fraction.
    const volScalarField alphaDiff =
        max(alphaTot - alphaMinFriction, scalar(0));
    const volScalarField gap = max(alphasMax - alphaTot, alphaDeltaMin_);

    const volScalarField pfTot =
        min(Fr_ * pow(alphaDiff, eta_) / pow(gap, p_),
            dimensionedScalar("maxPress", Fr_.dimensions(), 1e9));

    // Per-species correction:
    // only the alpha_s / alphaTot share of the total frictional pressure
    // is assigned to the current solid phase.
    return (alpha / alphaTotEff) * pfTot;
}

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::frictionalStressModels::JohnsonJacksonSchaeffer::
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

    // First build the original "total" frictional pressure:
    //
    //   pfTot = pfTot(alphaTot)
    //
    // i.e. the same law already present in the original code, expressed in
    // terms of the total solid volume fraction.
    const volScalarField alphaDiff =
        max(alphaTot - alphaMinFriction, scalar(0));
    const volScalarField gap = max(alphasMax - alphaTot, alphaDeltaMin_);

    const volScalarField pfTot =
        min(Fr_ * pow(alphaDiff, eta_) / pow(gap, p_),
            dimensionedScalar("maxPress", Fr_.dimensions(), 1e9));

    // Derivative of the total frictional pressure law with respect to alphaTot.
    //
    // This is essentially the same derivative already used in the original
    // code.
    const volScalarField pfTotPrime =
        Fr_
        * (eta_ * pow(alphaDiff, eta_ - 1) * gap + p_ * pow(alphaDiff, eta_))
        / pow(gap, p_ + 1);

    // Main correction:
    //
    // It is NOT sufficient to multiply the old derivative by (alpha/alphaTot).
    // The per-species frictional pressure is now defined as
    //
    //   pf_s = (alpha_s / alphaTot) * pfTot(alphaTot)
    //
    // therefore its derivative with respect to the same phase fraction alpha_s
    // is
    //
    //   d(pf_s)/d(alpha_s)
    // = d(alpha_s/alphaTot)/d(alpha_s) * pfTot
    //   + (alpha_s/alphaTot) * d(pfTot)/d(alphaTot)
    //
    // i.e.
    //
    //   d(pf_s)/d(alpha_s)
    // = ((alphaTot - alpha_s)/alphaTot^2) * pfTot
    //   + (alpha_s/alphaTot) * pfTotPrime
    //
    // This term is essential because frictionalPressurePrime() enters pPrime(),
    // and therefore contributes to both the phase-pressure force and the
    // pressure-momentum coupling of the solver.
    //
    // If frictionalPressure() were corrected without correcting
    // frictionalPressurePrime() consistently, the model would become internally
    // inconsistent: the frictional pressure and its effective gradient would no
    // longer match.
    return ((alphaTotEff - alpha) / sqr(alphaTotEff)) * pfTot
         + (alpha / alphaTotEff) * pfTotPrime;
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::frictionalStressModels::JohnsonJacksonSchaeffer::nu(
    const phaseModel& phase,
    const phaseModel& continuousPhase,
    const dimensionedScalar& alphaMinFriction,
    const volScalarField& alphasMax,
    const volScalarField& pf,
    const volScalarField& rho,
    const volSymmTensorField& D) const
{
    const volScalarField alphas = 1.0 - continuousPhase;

    tmp<volScalarField> tnu(volScalarField::New(
        IOobject::groupName(Foam::typedName<frictionalStressModel>("nu"),
                            phase.group()),
        phase.mesh(),
        dimensionedScalar(dimensionSet(0, 2, -1, 0, 0), 0)));

    volScalarField& nuf = tnu.ref();

    forAll(D, celli)
    {
        if (alphas[celli] > alphaMinFriction.value())
        {
            nuf[celli] =
                0.5 * pf[celli] / rho[celli] * sin(phi_.value())
                / (sqrt((1.0 / 3.0) * sqr(tr(D[celli])) - invariantII(D[celli]))
                   + small);
        }
    }

    const fvPatchList& patches = phase.mesh().boundary();
    volScalarField::Boundary& nufBf = nuf.boundaryFieldRef();

    forAll(patches, patchi)
    {
        if (!patches[patchi].coupled())
        {
            nufBf[patchi] =
                (pf.boundaryField()[patchi] / rho.boundaryField()[patchi]
                 * sin(phi_.value())
                 / (mag(phase.U()().boundaryField()[patchi].snGrad()) + small));
        }
    }

    // Correct coupled BCs
    nuf.correctBoundaryConditions();

    return tnu;
}


// ************************************************************************* //
