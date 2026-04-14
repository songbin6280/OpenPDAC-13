/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2022 OpenFOAM Foundation
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

#include "SrivastavaFrictionalStress.H"
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
defineTypeNameAndDebug(Srivastava, 0);

addToRunTimeSelectionTable(frictionalStressModel, Srivastava, dictionary);
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::frictionalStressModels::Srivastava::Srivastava(
    const dictionary& dict)
: frictionalStressModel(dict),
  coeffDict_(dict.optionalSubDict(typeName + "Coeffs")),
  Fr_("Fr", dimensionSet(1, -1, -2, 0, 0), coeffDict_),
  eta_("eta", dimless, coeffDict_), p_("p", dimless, coeffDict_),
  phi_("phi", dimless, coeffDict_),
  alphaDeltaMin_("alphaDeltaMin", dimless, coeffDict_)
{
    phi_ *= constant::mathematical::pi / 180.0;
    sinPhi_ = sin(phi_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::frictionalStressModels::Srivastava::~Srivastava() {}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// -------------------------------------------------------------------------
// Split-invariance correction for Srivastava frictional pressure
//
// Original behaviour:
//   The Srivastava frictional pressure law was evaluated from the total solid
//   volume fraction alphaTot and the same pf value was returned for every
//   solid phase.
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
Foam::kineticTheoryModels::frictionalStressModels::Srivastava::
    frictionalPressure(const phaseModel& phase,
                       const phaseModel& continuousPhase,
                       const dimensionedScalar& alphaMinFriction,
                       const volScalarField& alphasMax) const
{
    // alpha    = volume fraction of the current solid phase "s"
    // alphaTot = total solid volume fraction
    //
    // Original OpenPDAC behaviour:
    //   the Srivastava frictional pressure law was evaluated as a function
    //   of the total solids fraction alphaTot, and the same pf value was
    //   returned for every solid phase.
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
    //
    // The Srivastava model already contains a piecewise definition:
    //   1) a Johnson-Jackson-type branch,
    //   2) a linearized branch close to alphaMax - alphaDeltaMin.
    //
    // Here the same piecewise law is preserved, but it is interpreted as the
    // total frictional pressure of the solid mixture, which is then distributed
    // to the current species.

    const volScalarField& alpha = phase;
    const volScalarField alphaTot = 1.0 - continuousPhase;

    // Numerical safeguard to avoid division by zero in nearly empty cells.
    const dimensionedScalar alphaTotSmall("alphaTotSmall", dimless, SMALL);
    const volScalarField alphaTotEff = max(alphaTot, alphaTotSmall);

    // Total Johnson-Jackson-like branch evaluated from the total solids
    // fraction.
    const volScalarField JohnsonJacksonP =
        Fr_ * pow(max(alphaTot - alphaMinFriction, scalar(0)), eta_)
        / pow(max(alphasMax - alphaTot, alphaDeltaMin_), p_);

    // Total linearized branch used close to alphaMax - alphaDeltaMin.
    const volScalarField LinearCoeff =
        Fr_
        * (pow(alphasMax - alphaDeltaMin_ - alphaMinFriction, eta_)
               / pow(alphaDeltaMin_, p_)
           + eta_ * pow(alphasMax - alphaDeltaMin_ - alphaMinFriction, eta_ - 1)
                 / pow(alphaDeltaMin_, p_)
           + p_ * pow(alphasMax - alphaDeltaMin_ - alphaMinFriction, eta_)
                 / pow(alphaDeltaMin_, p_ - 1));

    const volScalarField LinearP =
        LinearCoeff * (alphaTot - alphasMax + alphaDeltaMin_);

    // Piecewise total frictional pressure law.
    const volScalarField pfTot =
        neg0(alphaTot - (alphasMax - alphaDeltaMin_)) * JohnsonJacksonP
        + pos(alphaTot - (alphasMax - alphaDeltaMin_)) * LinearP;

    // Per-species partition of the total frictional pressure.
    return (alpha / alphaTotEff) * pfTot;
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::frictionalStressModels::Srivastava::
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

    // Rebuild the same total piecewise frictional pressure law used in
    // frictionalPressure(), but interpreted as the total frictional pressure
    // of the solids mixture.
    const volScalarField JohnsonJacksonP =
        Fr_ * pow(max(alphaTot - alphaMinFriction, scalar(0)), eta_)
        / pow(max(alphasMax - alphaTot, alphaDeltaMin_), p_);

    const volScalarField LinearCoeff =
        Fr_
        * (pow(alphasMax - alphaDeltaMin_ - alphaMinFriction, eta_)
               / pow(alphaDeltaMin_, p_)
           + eta_ * pow(alphasMax - alphaDeltaMin_ - alphaMinFriction, eta_ - 1)
                 / pow(alphaDeltaMin_, p_)
           + p_ * pow(alphasMax - alphaDeltaMin_ - alphaMinFriction, eta_)
                 / pow(alphaDeltaMin_, p_ - 1));

    const volScalarField LinearP =
        LinearCoeff * (alphaTot - alphasMax + alphaDeltaMin_);

    const volScalarField pfTot =
        neg0(alphaTot - (alphasMax - alphaDeltaMin_)) * JohnsonJacksonP
        + pos(alphaTot - (alphasMax - alphaDeltaMin_)) * LinearP;

    // Derivative of the total Johnson-Jackson-like branch with respect to
    // alphaTot. This follows the same regularized expression already present in
    // the original implementation.
    const volScalarField JohnsonJacksonPprime =
        Fr_
        * (eta_ * pow(max(alphaTot - alphaMinFriction, scalar(0)), eta_ - 1)
               * max(alphasMax - alphaTot, alphaDeltaMin_)
           + p_ * pow(max(alphaTot - alphaMinFriction, scalar(0)), eta_)
                 * pos(alphasMax - alphaTot - alphaDeltaMin_))
        / pow(max(alphasMax - alphaTot, alphaDeltaMin_), p_ + 1);

    // Derivative of the linearized branch with respect to alphaTot.
    // Since LinearP = LinearCoeff * (alphaTot - alphasMax + alphaDeltaMin),
    // its derivative is simply LinearCoeff.
    const volScalarField LinearPprime = LinearCoeff;

    // Piecewise derivative of the total frictional pressure law.
    const volScalarField pfTotPrime =
        neg0(alphaTot - (alphasMax - alphaDeltaMin_)) * JohnsonJacksonPprime
        + pos(alphaTot - (alphasMax - alphaDeltaMin_)) * LinearPprime;

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
Foam::kineticTheoryModels::frictionalStressModels::Srivastava::nu(
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


bool Foam::kineticTheoryModels::frictionalStressModels::Srivastava::read()
{
    coeffDict_ <<= dict_.optionalSubDict(typeName + "Coeffs");

    Fr_.read(coeffDict_);
    eta_.read(coeffDict_);
    p_.read(coeffDict_);

    phi_.read(coeffDict_);
    phi_ *= constant::mathematical::pi / 180.0;

    alphaDeltaMin_.read(coeffDict_);

    return true;
}


// ************************************************************************* //
