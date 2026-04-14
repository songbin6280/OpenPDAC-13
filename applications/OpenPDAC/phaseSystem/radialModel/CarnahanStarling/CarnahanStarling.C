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

#include "CarnahanStarling.H"
#include "addToRunTimeSelectionTable.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace radialModels
{
defineTypeNameAndDebug(CarnahanStarling, 0);

addToRunTimeSelectionTable(radialModel, CarnahanStarling, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radialModels::CarnahanStarling::CarnahanStarling(const dictionary& dict,
                                                       const phaseSystem& fluid)
: radialModel(dict, fluid)
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radialModels::CarnahanStarling::~CarnahanStarling() {}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::PtrList<Foam::volScalarField> Foam::radialModels::CarnahanStarling::g0(
    const phaseModel& phasei,
    const phaseModel& continuousPhase,
    const dimensionedScalar& alphaMinFriction,
    const volScalarField& alphasMax) const
{
    const volScalarField& alphai = phasei;
    const label& indexi = phasei.index();
    const phaseSystem& fluid = phasei.fluid();

    PtrList<volScalarField> g0_im(fluid.phases().size());

    volScalarField alphas = alphai;
    volScalarField eta2 = alphai / phasei.d();

    forAll(fluid.phases(), phaseIdx)
    {
        const phaseModel& phase = fluid.phases()[phaseIdx];
        if ((&phase != &continuousPhase) and !(phaseIdx == indexi))
        {
            alphas += phase;
            eta2 += phase / phase.d();
        }
    }

    const volScalarField denominatorTerm =
        max(1.0
                - alphas
                      / max(alphasMax,
                            dimensionedScalar("small", dimless, ROOTVSMALL)),
            dimensionedScalar("small", dimless, ROOTVSMALL));

    forAll(g0_im, iter)
    {
        const phaseModel& phasem = fluid.phases()[iter];

        if (&phasem != &continuousPhase)
        {
            const volScalarField di = phasei.d();
            const volScalarField dm = phasem.d();
            const dimensionedScalar smallD("smallD", dimLength, ROOTVSMALL);
            volScalarField term_d = di * dm / (di + dm + smallD);

            g0_im.set(iter,
                      new volScalarField(
                          "g0_im" + phasei.name() + "_" + phasem.name(),
                          1.0 / denominatorTerm
                              + 3.0 * term_d * eta2 / sqr(denominatorTerm)
                              + 2.0 * sqr(term_d) * sqr(eta2)
                                    / pow3(denominatorTerm)));
        }
    }

    return g0_im;
}

Foam::PtrList<Foam::volScalarField>
Foam::radialModels::CarnahanStarling::g0prime(
    const phaseModel& phasei,
    const phaseModel& continuousPhase,
    const dimensionedScalar& alphaMinFriction,
    const volScalarField& alphasMax) const
{
    const phaseSystem& fluid = phasei.fluid();
    const label indexi = phasei.index();

    PtrList<volScalarField> g0prime_im(fluid.phases().size());

    // Total solids volume fraction: alpha_s = sum_j alpha_j = 1 - alpha_g
    volScalarField alphas = phasei;

    // Mixture moment: eta2 = sum_j(alpha_j/d_j)
    volScalarField eta2 = phasei / phasei.d();

    forAll(fluid.phases(), phaseIdx)
    {
        const phaseModel& phase = fluid.phases()[phaseIdx];

        if ((&phase != &continuousPhase) && (phaseIdx != indexi))
        {
            alphas += phase;
            eta2 += phase / phase.d();
        }
    }

    // D = 1 - alpha_s/alpha_s,max, with the same lower clamp used in g0()
    volScalarField denominatorTerm(
        max(scalar(1) - alphas / max(alphasMax, scalar(ROOTVSMALL)),
            scalar(ROOTVSMALL)));

    // Active branch of the denominator clamp.
    // If the clamp is active, D is frozen and dD/dalpha_i must be zero.
    volScalarField activeDenom(
        pos((scalar(1) - alphas / max(alphasMax, scalar(ROOTVSMALL)))
            - scalar(ROOTVSMALL)));

    // Exact derivative of eta2 with respect to alpha_i:
    // d(eta2)/d(alpha_i) = 1/d_i
    volScalarField dEta2dAlphai = scalar(1) / phasei.d();

    // dD/dalpha_i = -1/alpha_s,max, switched off when the clamp is active
    volScalarField dDdAlphai =
        -activeDenom / max(alphasMax, scalar(ROOTVSMALL));

    forAll(g0prime_im, iter)
    {
        const phaseModel& phasem = fluid.phases()[iter];

        if (&phasem != &continuousPhase)
        {
            const volScalarField di = phasei.d();
            const volScalarField dm = phasem.d();

            // Pair-size factor:
            // t_im = d_i d_m / (d_i + d_m)
            const dimensionedScalar smallD("smallD", dimLength, ROOTVSMALL);
            volScalarField term_d = di * dm / (di + dm + smallD);

            // g0_{i,m} = 1/D + 3 t_im eta2 / D^2 + 2 t_im^2 eta2^2 / D^3
            // and g0prime_im[m] = d g0_{i,m} / d alpha_i

            tmp<volScalarField> dT1(-dDdAlphai / sqr(denominatorTerm));

            tmp<volScalarField> dT2(
                3.0 * term_d
                * (dEta2dAlphai / sqr(denominatorTerm)
                   - 2.0 * eta2 * dDdAlphai / pow3(denominatorTerm)));

            tmp<volScalarField> dT3(
                2.0 * sqr(term_d)
                * (2.0 * eta2 * dEta2dAlphai / pow3(denominatorTerm)
                   - 3.0 * sqr(eta2) * dDdAlphai / pow(denominatorTerm, 4.0)));

            g0prime_im.set(iter,
                           new volScalarField("g0prime_im" + phasei.name() + "_"
                                                  + phasem.name(),
                                              dT1 + dT2 + dT3));
        }
    }

    return g0prime_im;
}

// ************************************************************************* //
