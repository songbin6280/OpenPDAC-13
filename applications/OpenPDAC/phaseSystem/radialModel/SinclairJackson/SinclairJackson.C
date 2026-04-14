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

#include "SinclairJackson.H"
#include "addToRunTimeSelectionTable.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace radialModels
{
defineTypeNameAndDebug(SinclairJackson, 0);

addToRunTimeSelectionTable(radialModel, SinclairJackson, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radialModels::SinclairJackson::SinclairJackson(
    const dictionary& coeffDict, const phaseSystem& fluid)
: radialModel(coeffDict, fluid)
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radialModels::SinclairJackson::~SinclairJackson() {}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::PtrList<Foam::volScalarField> Foam::radialModels::SinclairJackson::g0(
    const phaseModel& phasei,
    const phaseModel& continuousPhase,
    const dimensionedScalar& alphaMinFriction,
    const volScalarField& alphasMax) const
{
    const volScalarField& alphai = phasei;
    const label& indexi = phasei.index();
    const phaseSystem& fluid = phasei.fluid();

    PtrList<volScalarField> g0_mm(fluid.phases().size());
    PtrList<volScalarField> g0_im(fluid.phases().size());

    volScalarField const_sum = alphai / phasei.d();

    forAll(fluid.phases(), phaseIdx)
    {
        const phaseModel& phase = fluid.phases()[phaseIdx];

        if ((&phase != &continuousPhase) and !(phaseIdx == indexi))
        {
            const volScalarField& alpha = phase;
            const_sum += alpha / phase.d();
        }
    }

    volScalarField alphas = 1.0 - continuousPhase;

    volScalarField g0 =
        1.0 / (1 - cbrt(min(alphas, alphaMinFriction) / alphasMax));

    forAll(fluid.phases(), phaseIdx)
    {
        const phaseModel& phase = fluid.phases()[phaseIdx];

        if (&phase != &continuousPhase)
        {
            g0_mm.set(
                phaseIdx,
                volScalarField("g0_mm" + phasei.name() + "_" + phase.name(),
                               g0 + 0.5 * phase.d() * const_sum));
        }
    }

    forAll(g0_im, iter)
    {
        const phaseModel& phase = fluid.phases()[iter];

        if (&phase != &continuousPhase)
        {
            g0_im.set(iter,
                      volScalarField(
                          "g0_im" + phasei.name() + "_" + phase.name(),
                          (phase.d() * g0_mm[indexi] + phasei.d() * g0_mm[iter])
                              / (phasei.d() + phase.d())));
        }
    }

    return g0_im;
}


Foam::PtrList<Foam::volScalarField>
Foam::radialModels::SinclairJackson::g0prime(
    const phaseModel& phasei,
    const phaseModel& continuousPhase,
    const dimensionedScalar& alphaMinFriction,
    const volScalarField& alphasMax) const
{
    // alpha_i of the phase for which the derivative is taken.
    // In the current OpenPDAC design, this function returns:
    //
    //     g0prime_im[m] = d g0_{i,m} / d alpha_i
    //
    // i.e. the derivative of all pair radial functions involving phase i
    // with respect to the volume fraction of the same phase i.
    const volScalarField& alphai = phasei;
    const label& indexi = phasei.index();
    const phaseSystem& fluid = phasei.fluid();

    PtrList<volScalarField> g0prime_mm(fluid.phases().size());
    PtrList<volScalarField> g0prime_im(fluid.phases().size());

    // ---------------------------------------------------------------------
    // Build:
    //
    //     sum_j (alpha_j / d_j)
    //
    // over all solid phases.
    //
    // This is the polydisperse correction used by the Sinclair-Jackson model.
    // For a fixed phase m, the "diagonal" quantity is:
    //
    //     g0_mm^(m) = g0_base + 0.5 * d_m * sum_j(alpha_j/d_j)
    //
    // Therefore:
    //
    //     d/dalpha_i [g0_mm^(m)]
    //   = d(g0_base)/dalpha_i + 0.5 * d_m * d/dalpha_i[sum_j(alpha_j/d_j)]
    //   = g0prime_base + 0.5 * d_m / d_i
    //
    // because:
    //
    //     d/dalpha_i [sum_j(alpha_j/d_j)] = 1/d_i .
    // ---------------------------------------------------------------------
    volScalarField const_sum = alphai / phasei.d();

    forAll(fluid.phases(), phaseIdx)
    {
        const phaseModel& phase = fluid.phases()[phaseIdx];

        if ((&phase != &continuousPhase) && !(phaseIdx == indexi))
        {
            const volScalarField& alpha = phase;
            const_sum += alpha / phase.d();
        }
    }

    // Total solids volume fraction.
    // The base Sinclair-Jackson singularity depends only on alpha_s,tot.
    volScalarField alphas = 1.0 - continuousPhase;

    // ---------------------------------------------------------------------
    // Base radial function:
    //
    //     g0_base = 1 / (1 - cbrt(min(alpha_s,tot,
    //     alphaMinFriction)/alphasMax))
    //
    // The derivative below is the derivative of the active (unclamped) branch.
    // The clamp is then enforced by the mask "posCoeff", which sets the
    // derivative to zero outside the interval where the base function varies.
    // ---------------------------------------------------------------------
    volScalarField aByaMax(
        cbrt(min(max(alphas, scalar(1e-3)), alphaMinFriction) / alphasMax));

    volScalarField g0primeBase((1.0 / (3.0 * alphasMax))
                               / sqr(aByaMax - sqr(aByaMax)));

    // The base term is constant below 1e-3 and above alphaMinFriction because
    // of the clamp, so its derivative must vanish there.
    volScalarField posCoeff(pos(alphaMinFriction - alphas)
                            * pos(alphas - scalar(1e-3)));

    g0primeBase *= posCoeff;

    // ---------------------------------------------------------------------
    // Derivative of the "diagonal" pair quantity g0_mm^(m):
    //
    //     g0_mm^(m) = g0_base + 0.5 * d_m * sum_j(alpha_j/d_j)
    //
    // with respect to alpha_i:
    //
    //     d g0_mm^(m) / d alpha_i = g0primeBase + 0.5 * d_m / d_i
    //
    // IMPORTANT:
    // - phasei.d() is d_i, the diameter of the phase with respect to which
    //   we differentiate;
    // - phase.d()  is d_m, the diameter of the phase labeling g0_mm^(m).
    // ---------------------------------------------------------------------
    forAll(fluid.phases(), phaseIdx)
    {
        const phaseModel& phase = fluid.phases()[phaseIdx];

        if (&phase != &continuousPhase)
        {
            g0prime_mm.set(
                phaseIdx,
                new volScalarField("g0prime_mm" + phasei.name() + "_"
                                       + phase.name(),
                                   g0primeBase + 0.5 * phase.d() / phasei.d()));
        }
    }

    // ---------------------------------------------------------------------
    // Off-diagonal pair quantity used in g0():
    //
    //     g0_{i,m} = [ d_m * g0_mm^(i) + d_i * g0_mm^(m) ] / (d_i + d_m)
    //
    // Since d_i and d_m are treated as parameters here, the derivative with
    // respect to alpha_i must preserve the same linear weights:
    //
    //     d g0_{i,m} / d alpha_i
    //   = [ d_m * d(g0_mm^(i))/dalpha_i + d_i * d(g0_mm^(m))/dalpha_i ]
    //     / (d_i + d_m)
    //
    // Therefore the weights in g0prime_im must be the same as in g0_im.
    // ---------------------------------------------------------------------
    forAll(g0prime_im, iter)
    {
        const phaseModel& phase = fluid.phases()[iter];

        if (&phase != &continuousPhase)
        {
            g0prime_im.set(iter,
                           new volScalarField("g0prime_im" + phasei.name() + "_"
                                                  + phase.name(),
                                              (phase.d() * g0prime_mm[indexi]
                                               + phasei.d() * g0prime_mm[iter])
                                                  / (phasei.d() + phase.d())));
        }
    }

    return g0prime_im;
}


// ************************************************************************* //
