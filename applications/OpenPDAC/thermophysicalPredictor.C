/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2025 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenPDAC.
    This file was derived from the multiphaseEuler solver in OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
    Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program. If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "OpenPDAC.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvcSup.H"
#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmSup.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::solvers::OpenPDAC::compositionPredictor()
{
    forAll(fluid.multicomponentPhases(), multicomponentPhasei)
    {
        phaseModel& phase = fluid_.multicomponentPhases()[multicomponentPhasei];

        UPtrList<volScalarField>& Y = phase.YRef();
        const volScalarField& alpha = phase;
        const volScalarField& rho = phase.rho();

        forAll(Y, i)
        {
            if (phase.solveSpecie(i))
            {
                fvScalarMatrix YiEqn(phase.YiEqn(Y[i])
                                     == fvModels().source(alpha, rho, Y[i]));

                YiEqn.relax();

                fvConstraints().constrain(YiEqn);

                YiEqn.solve("Yi");

                Y[i] = max(0.0, min(Y[i], 1.0));

                fvConstraints().constrain(Y[i]);
            }
            else
            {
                Y[i].correctBoundaryConditions();
            }
        }
    }

    fluid_.correctSpecies();
}


void Foam::solvers::OpenPDAC::energyPredictor()
{
    // Retrieve heat transfer matrices from the heat transfer system
    autoPtr<HashPtrTable<fvScalarMatrix>> heatTransferPtr =
        heatTransferSystem_.heatTransfer();
    HashPtrTable<fvScalarMatrix>& heatTransfer = heatTransferPtr();

    // Read switches for drag energy correction
    const bool dragEnergyCorrection =
        pimple.dict().lookupOrDefault<Switch>("dragEnergyCorrection", false);
    const bool totalEnergy =
        pimple.dict().lookupOrDefault<Switch>("totalEnergy", false);

    // Iterate over all thermal phases to solve their energy equations
    forAll(fluid.thermalPhases(), thermalPhasei)
    {
        phaseModel& phase = fluid_.thermalPhases()[thermalPhasei];

        const volScalarField& alpha = phase;
        const volScalarField& rho = phase.rho();

        // Construct the energy equation
        // Includes heat transfer between phases and fvModels sources
        fvScalarMatrix EEqn(
            phase.heEqn()
            == *heatTransfer[phase.name()]
                   + fvModels().source(alpha, rho, phase.thermo().he()));

        // Apply drag energy correction if enabled
        if (dragEnergyCorrection)
        {
            if (totalEnergy)
            {
                // Include kinetic energy transfer due to drag
                PtrList<volScalarField> dragEnergyTransfers(
                    movingPhases.size());
                momentumTransferSystem_.dragEnergy(dragEnergyTransfers);
                EEqn -= dragEnergyTransfers[thermalPhasei];
            }
            else
            {
                // Include frictional dissipation heating due to drag
                PtrList<volScalarField> dragDissipations(movingPhases.size());
                momentumTransferSystem_.dragDissipation(dragDissipations);
                EEqn -= dragDissipations[thermalPhasei];
            }
        }


        // Relax, constrain and solve the energy equation
        EEqn.relax();
        fvConstraints().constrain(EEqn);
        EEqn.solve();
        fvConstraints().constrain(phase.thermo().he());
    }

    // Update thermodynamic properties based on new enthalpy/temperature
    fluid_.correctThermo();

    // Initialize dmdts with zeros (no mass transfer)
    // This is required for continuity error correction
    PtrList<volScalarField::Internal> dmdts(fluid.phases().size());
    forAll(dmdts, i)
    {
        dmdts.set(i,
                  new volScalarField::Internal(
                      IOobject("dmdt", runTime.name(), mesh),
                      mesh,
                      dimensionedScalar(dimDensity / dimTime, 0)));
    }

    // Correct continuity errors using the initialized mass transfer rates
    fluid_.correctContinuityError(dmdts);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::OpenPDAC::thermophysicalPredictor()
{
    // If there are no thermal phases, there is no energy equation to solve.
    // Skip this entire step.
    if (fluid.thermalPhases().empty())
    {
        return;
    }

    // Bypass logic for "ghost" iterations in the PIMPLE loop
    if (forceFinalPimpleIter_ && !pimple.finalIter())
    {
        return;
    }

    // Dynamic number of correctors: fewer at start, more at the end
    int dynamicMaxCorrectors =
        min(nMaxEnergyCorrectors, nEnergyCorrectors + pimpleIter);

    if (pimple.finalIter())
    {
        dynamicMaxCorrectors = nMaxEnergyCorrectors;
    }

    // Store residual from previous iteration
    autoPtr<volScalarField> heDisplPrev;

    const word& continuousPhaseName_ = fluid.continuousPhaseName();

    for (int Ecorr = 0; Ecorr < dynamicMaxCorrectors; Ecorr++)
    {
        const phaseModel& continuousPhase =
            fluid.phases()[continuousPhaseName_];

        // Store the starting enthalpy of the continuous phase (h_k)
        volScalarField heStart(continuousPhase.thermo().he());

        // Pre-corrections for dispersed phases (e.g., temperature clipping)
        forAll(fluid.thermalPhases(), thermalPhasei)
        {
            const phaseModel& phase = fluid.thermalPhases()[thermalPhasei];
            volScalarField& heNew =
                const_cast<volScalarField&>(phase.thermo().he());

            if (&phase == &continuousPhase)
            {
                // For the continuous phase, clip negative enthalpy to a minimum
                // physical value to prevent crashes from non-physical T.
                volScalarField Tmin(continuousPhase.thermo().T());
                Tmin = dimensionedScalar(dimTemperature, 273.0);
                volScalarField heMin = phase.thermo().he(p_, Tmin);

                heNew =
                    pos(heNew - heMin) * heNew + neg0(heNew - heMin) * heMin;
            }
            else if (correctTdispersed) // This is a dispersed phase
            {
                volScalarField heTcont =
                    phase.thermo().he(p_, continuousPhase.thermo().T());

                heNew = pos(heNew) * heNew + neg0(heNew) * heTcont;

                heNew = pos0(phase - phase.residualAlpha()) * heNew
                      + neg(phase - phase.residualAlpha()) * heTcont;
            }
        }

        // Standard OpenFOAM predictor sequence
        fluid_.correctThermo();
        fluid_.predictThermophysicalTransport();
        compositionPredictor();

        // Solve the linear system. OpenFOAM prints "Linear Residuals" here.
        energyPredictor();

        // ======================================================
        // AITKEN ACCELERATION (Continuous Phase Only)
        // Applied to the Displacement Vector (Delta h)
        // ======================================================

        if (useAitkenRelaxation_)
        {
            // 'he' is the tentative value predicted by Picard (tilde_h_{k+1})
            volScalarField& he =
                const_cast<volScalarField&>(continuousPhase.thermo().he());

            // Calculate current DISPLACEMENT: Delta_h_k = tilde_h_{k+1} - h_k
            volScalarField heDisplCurr = he - heStart;

            // Apply relaxation only if history exists (from 2nd iteration
            // onwards)
            if (Ecorr > 0 && heDisplPrev.valid())
            {
                // Calculate the change in displacement: delta = Delta_h_k -
                // Delta_h_{k-1}
                volScalarField deltaDispl = heDisplCurr - heDisplPrev();

                // Compute scalar products using primitiveField() to avoid
                // dimension checks
                scalar num = gSumProd(heDisplCurr.primitiveField(),
                                      deltaDispl.primitiveField());
                scalar den = gSumSqr(deltaDispl.primitiveField());

                if (den > SMALL)
                {
                    scalar omega = 1.0 - (num / den);

                    // Stable clamp for stiff thermodynamic systems
                    omega = max(0.2, min(omega, 1.2));

                    Info << "  Continuous Aitken omega: " << omega << endl;

                    // Final Update: h_{k+1} = h_k + omega * Delta_h_k
                    he = heStart + omega * heDisplCurr;
                }
            }

            // Update displacement history for the next iteration
            if (!heDisplPrev.valid())
            {
                heDisplPrev.reset(new volScalarField(heDisplCurr));
            }
            else
            {
                heDisplPrev() = heDisplCurr;
            }

            // Ensure boundary conditions are correct for the history field
            heDisplPrev().correctBoundaryConditions();

            // Update thermodynamics with the accelerated field (if modified)
            if (Ecorr > 0)
            {
                fluid_.correctThermo();
            }
        }

        // ======================================================

        // Report Min/Max Temperatures
        forAll(fluid.thermalPhases(), thermalPhasei)
        {
            const phaseModel& phase = fluid.thermalPhases()[thermalPhasei];

            Info << phase.name() << " min/max T "
                 << min(phase.thermo().T()).value() << " - "
                 << max(phase.thermo().T()).value() << endl;
        }

        // Check LINEAR RESIDUALS (r0) for Early Exit
        bool checkResidual(true);
        bool doCheck(false);

        forAll(fluid.thermalPhases(), thermalPhasei)
        {
            const phaseModel& phase = fluid.thermalPhases()[thermalPhasei];

            word name(phase.thermo().he().name());
            const DynamicList<SolverPerformance<scalar>>& sp(
                Residuals<scalar>::field(mesh, name));
            label n = sp.size();

            scalar r0 = cmptMax(sp[n - 1].initialResidual());
            Info << name << " initial residual " << r0 << endl;

            if (energyControlDict.found(name))
            {
                doCheck = true;
                scalar residual(energyControlDict.lookup<scalar>(name));
                checkResidual = checkResidual && (r0 <= residual);
            }
        }

        Info << "Iteration " << Ecorr + 1
             << " Check for initial Energy Residual " << checkResidual << endl;

        // If linear residuals are small enough, break the inner loop
        if (doCheck)
        {
            if (checkResidual)
            {
                convergenceFlag = true;
                break;
            }
            else
            {
                convergenceFlag = false;
            }
            Info << "convergenceFlag = " << convergenceFlag << endl;
        }
    }
}
// ************************************************************************* //
