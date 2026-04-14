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

#include "SchaefferFrictionalStress.H"
#include "addToRunTimeSelectionTable.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace frictionalStressModels
{
defineTypeNameAndDebug(Schaeffer, 0);

addToRunTimeSelectionTable(frictionalStressModel, Schaeffer, dictionary);
}
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::kineticTheoryModels::frictionalStressModels::Schaeffer::readCoeffs(
    const dictionary& coeffDict)
{
    phi_.read(coeffDict);
    phi_ *= constant::mathematical::pi / 180.0;

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::frictionalStressModels::Schaeffer::Schaeffer(
    const dictionary& coeffDict)
: frictionalStressModel(coeffDict), phi_("phi", dimless, coeffDict)
{
    phi_ *= constant::mathematical::pi / 180.0;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::frictionalStressModels::Schaeffer::~Schaeffer() {}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// -------------------------------------------------------------------------
// Split-invariance correction for Schaeffer frictional pressure
//
// Original behaviour:
//   The Schaeffer frictional pressure law was evaluated from the total solid
//   volume fraction alphaTot and the same pf value was returned for every
//   solid phase.
//
// Problem:
//   If one solid phase is artificially split into two identical phases,
//   each phase receives the same total frictional pressure contribution.
//   As a consequence, the sum of the per-phase frictional contributions is
//   no longer consistent with the single-phase case.
//
// Correction introduced here:
//   pf_s = (alpha_s / alphaTot) * pfTot(alphaTot)
//
// Rationale:
//   This preserves the original Schaeffer model in the single-solid-phase
//   limit (alpha_s = alphaTot), while distributing the total frictional
//   pressure consistently among identical solid phases in split-phase cases.
//
// Important consequence:
//   frictionalPressurePrime() must be corrected consistently as the derivative
//   of the modified per-phase pressure, otherwise pPrime() becomes
//   inconsistent.
// -------------------------------------------------------------------------

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::frictionalStressModels::Schaeffer::
    frictionalPressure(const phaseModel& phase,
                       const phaseModel& continuousPhase,
                       const dimensionedScalar& alphaMinFriction,
                       const volScalarField& alphasMax) const
{
    // alpha    = volume fraction of the current solid phase "s"
    // alphaTot = total solid volume fraction
    //
    // Original behaviour:
    //   the Schaeffer frictional pressure was evaluated as a function of the
    //   total solids volume fraction and the same pf value was returned for
    //   every solid phase.
    //
    // This is not consistent with a phase-splitting test:
    //
    //   case A: one solid phase      -> pf = pfTot(alphaTot)
    //   case B: two identical phases -> pf1 = pf2 = pfTot(alphaTot)
    //
    // In case B, the sum of the per-phase frictional contributions no longer
    // matches the single-phase result.
    //
    // To recover consistency, the total frictional pressure is first computed
    // from alphaTot and then partitioned among the solid phases according to
    // their relative volume fraction:
    //
    //   pf_s = (alpha_s / alphaTot) * pfTot(alphaTot)
    //
    // This guarantees that:
    // - for a single solid phase, alpha_s = alphaTot and the original
    //   Schaeffer formulation is recovered exactly;
    // - for two identical split phases, each phase receives only its proper
    //   share of the total frictional pressure.
    const volScalarField& alpha = phase;
    const volScalarField alphaTot = 1.0 - continuousPhase;

    // Numerical safeguard to avoid division by zero in nearly empty cells.
    const dimensionedScalar alphaTotSmall("alphaTotSmall", dimless, SMALL);
    const volScalarField alphaTotEff = max(alphaTot, alphaTotSmall);

    // Original Schaeffer total frictional pressure law, evaluated using the
    // total solids fraction.
    const volScalarField pfTot =
        dimensionedScalar(dimensionSet(1, -1, -2, 0, 0), 1e24)
        * pow(Foam::max(alphaTot - alphasMax, scalar(0)), 10.0);

    // Per-species partition of the total frictional pressure.
    return (alpha / alphaTotEff) * pfTot;
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::frictionalStressModels::Schaeffer::
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

    // Build the original total Schaeffer frictional pressure:
    //
    //   pfTot = 1e24 * max(alphaTot - alphasMax, 0)^10
    //
    // This is exactly the same total law as in the original implementation,
    // but now it is interpreted as the total frictional pressure of the solid
    // mixture rather than the pressure to be assigned independently to each
    // solid phase.
    const volScalarField pfTot =
        dimensionedScalar(dimensionSet(1, -1, -2, 0, 0), 1e24)
        * pow(Foam::max(alphaTot - alphasMax, scalar(0)), 10.0);

    // Derivative of the total Schaeffer law with respect to alphaTot:
    //
    //   d(pfTot)/d(alphaTot) = 1e25 * max(alphaTot - alphasMax, 0)^9
    //
    // which is the same derivative already present in the original code.
    const volScalarField pfTotPrime =
        dimensionedScalar(dimensionSet(1, -1, -2, 0, 0), 1e25)
        * pow(Foam::max(alphaTot - alphasMax, scalar(0)), 9.0);

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
Foam::kineticTheoryModels::frictionalStressModels::Schaeffer::nu(
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
