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
#include "localEulerDdtScheme.H"
#include "surfaceFields.H"
#include "fvcDiv.H"
#include "fvcSurfaceIntegrate.H"
#include "fvcMeshPhi.H"
#include "addToRunTimeSelectionTable.H"
#include "myHydrostaticInitialisation.H"
// #include "polyTopoChangeMap.H"
// #include "polyMeshMap.H"
// #include "collidingCloud.H"
#include "IOobjectList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
defineTypeNameAndDebug(OpenPDAC, 0);
addToRunTimeSelectionTable(solver, OpenPDAC, fvMesh);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::solvers::OpenPDAC::read()
{
    fluidSolver::read();

    predictMomentum =
        pimple.dict().lookupOrDefault<bool>("momentumPredictor", false);

    faceMomentum = pimple.dict().lookupOrDefault<Switch>("faceMomentum", false);

    dragCorrection =
        pimple.dict().lookupOrDefault<Switch>("dragCorrection", false);

    nEnergyCorrectors =
        pimple.dict().lookupOrDefault<int>("nEnergyCorrectors", 1);

    alphaControls.read(mesh.solution().solverDict("alpha"));

    lowPressureTimestepCorrection = pimple.dict().lookupOrDefault<Switch>(
        "lowPressureTimestepCorrection", false);

    tempOscillationControl_ = pimple.dict().lookupOrDefault<Switch>(
        "temperatureOscillationControl", false);

    if (tempOscillationControl_)
    {
        maxTempOscillation_ =
            pimple.dict().lookupOrDefault<scalar>("maxTempOscillation", 0.1);
    }

    correctTdispersed =
        pimple.dict().lookupOrDefault<Switch>("correctTdispersed", false);

    nonOrthogonalResidual =
        pimple.dict().lookupOrDefault<scalar>("nonOrthogonalResidual", 0.0);

    innerResidual = pimple.dict().lookupOrDefault<scalar>("innerResidual", 0.0);

    residualRatio =
        pimple.dict().lookupOrDefault<scalar>("residualRatio", 10.0);

    minOuterCorrectors =
        pimple.dict().lookupOrDefault<label>("minOuterCorrectors", 4);

    nMaxEnergyCorrectors =
        pimple.dict().lookupOrDefault<label>("nMaxEnergyCorrectors", 20);

    useAitkenRelaxation_ =
        pimple.dict().lookupOrDefault<Switch>("useAitkenRelaxation", false);

    if (pimple.dict().found("energyControl"))
    {
        energyControlDict = pimple.dict().subDict("energyControl");
    }


    return true;
}


void Foam::solvers::OpenPDAC::correctCoNum()
{
    scalarField sumPhi(fvc::surfaceSum(mag(phi))().primitiveField());

    forAll(movingPhases, movingPhasei)
    {
        sumPhi = max(sumPhi,
                     fvc::surfaceSum(mag(movingPhases[movingPhasei].phi()))()
                         .primitiveField());
    }

    if (lowPressureTimestepCorrection)
    {
        volScalarField alphasMax = fluid_.alfasMax();
        const word& continuousPhaseName = fluid.continuousPhaseName();
        volScalarField alfaCont = fluid.phases()[continuousPhaseName];

        scalarField alfa_ratio =
            pow(max(0 * alphasMax, alphasMax - alfaCont) / alphasMax, 0.5);

        Info << "p_ratio = " << p_ratio << endl;
        Info << "alfa_ratio: min = " << min(alfa_ratio) << endl;
        sumPhi /= sqrt(p_ratio);
    }

    if (tempOscillationControl_ && tempOscillationFactor_ < 1.0)
    {
        Info << "tempOscillationFactor = " << tempOscillationFactor_ << endl;
        sumPhi /= sqrt(tempOscillationFactor_);
    }


    CoNum_ =
        0.5 * gMax(sumPhi / mesh.V().primitiveField()) * runTime.deltaTValue();

    const scalar meanCoNum = 0.5
                           * (gSum(sumPhi) / gSum(mesh.V().primitiveField()))
                           * runTime.deltaTValue();

    if (lowPressureTimestepCorrection)
    {
        Info << "Courant Number mean: " << meanCoNum * sqrt(p_ratio)
             << " max: " << CoNum * sqrt(p_ratio) << endl;
        Info << "Modified Courant Number mean: " << meanCoNum
             << " max: " << CoNum << endl;
    }
    else
    {
        Info << "Courant Number mean: " << meanCoNum << " max: " << CoNum
             << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::OpenPDAC::OpenPDAC(fvMesh& mesh)
: fluidSolver(mesh),

  predictMomentum(
      pimple.dict().lookupOrDefault<Switch>("momentumPredictor", false)),

  faceMomentum(pimple.dict().lookupOrDefault<Switch>("faceMomentum", false)),

  dragCorrection(
      pimple.dict().lookupOrDefault<Switch>("dragCorrection", false)),

  nEnergyCorrectors(pimple.dict().lookupOrDefault<int>("nEnergyCorrectors", 1)),

  lowPressureTimestepCorrection(pimple.dict().lookupOrDefault<Switch>(
      "lowPressureTimestepCorrection", false)),

  tempOscillationControl_(false),

  maxTempOscillation_(0.1),

  tempOscillationFactor_(1.0),

  buoyancy(mesh),

  fluid_(mesh),

  phases_(fluid_.phases()),

  movingPhases_(fluid_.movingPhases()),

  phi_(fluid_.phi()),

  momentumTransferSystem_(fluid_),

  heatTransferSystem_(fluid_),

  p_(movingPhases_[0].fluidThermo().p()),

  p_rgh_(buoyancy.p_rgh),

  rho(IOobject(
          "rho", runTime.name(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
      fluid_.rho()),

  carrierIdx(0),

  muMix(IOobject("muMix",
                 runTime.name(),
                 mesh,
                 IOobject::NO_READ,
                 IOobject::AUTO_WRITE),
        mesh),

  U(IOobject(
        "U", runTime.name(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh),

  // Initialize cloud
  clouds(parcelClouds::New(mesh, rho, U, muMix, buoyancy.g)),

  pressureReference(p_, p_rgh_, pimple.dict(), fluid_.incompressible()),

  MRF(fluid_.MRF()),

  dmdts_(fluid_.phases().size()),

  fluid(fluid_), phases(phases_), movingPhases(movingPhases_),
  momentumTransfer(momentumTransferSystem_), heatTransfer(heatTransferSystem_),
  p(p_), p_rgh(p_rgh_), phi(phi_)
{
    // Read the controls
    read();

    if (tempOscillationControl_)
    {
        const label nThermal = fluid_.thermalPhases().size();
        if (nThermal > 0)
        {
            minTempMin_.setSize(nThermal);
            maxTempMin_.setSize(nThermal);
            sumTempMin_.setSize(nThermal);
            pimpleTempIterCount_.setSize(nThermal);
        }
        else
        {
            Info << "temperatureOscillationControl is active, but no thermal "
                 << "phases were found." << endl;
            tempOscillationControl_ = false;
        }
    }

    mesh.schemes().setFluxRequired(p_rgh.name());

    // create ph_rgh (p_rgh for hydrostatic pressure)
    volScalarField& ph_rgh = regIOobject::store(new volScalarField(
        IOobject("ph_rgh", "0", mesh, IOobject::MUST_READ, IOobject::NO_WRITE),
        mesh));

    // Initialization of hydrostatic pressure profile
    hydrostaticInitialisation(p_rgh_,
                              ph_rgh,
                              p_,
                              buoyancy.g,
                              buoyancy.hRef,
                              buoyancy.gh,
                              buoyancy.ghf,
                              fluid_,
                              pimple.dict());


    // Correct mixture thermodynamics with new pressure
    fluid_.correctThermo();
    rho = fluid_.rho();

    Info << "hRef " << buoyancy.hRef.value() << endl;

    Info << "min p " << min(p_).value() << " max p " << max(p_).value() << endl;

    p_ratio = min(p_).value() / p_.weightedAverage(mesh_.V()).value();

    Info << "min p_rgh " << min(p_rgh).value() << " max p_rgh "
         << max(p_rgh).value() << endl;
    Info << "min rho " << min(rho).value() << " max rho " << max(rho).value()
         << endl;

    // Carrier phase viscosity
    const word& continuousPhaseName = fluid.continuousPhaseName();

    const volScalarField& muC = phases_[continuousPhaseName].fluidThermo().mu();
    Info << "min muC " << min(muC).value() << " max muC " << max(muC).value()
         << endl;

    volScalarField alphasMax = fluid_.alfasMax();
    Info << "min alphasMax " << min(alphasMax).value() << " max alphasMax "
         << max(alphasMax).value() << endl;

    // Mixture viscosity
    volScalarField alphaS = 1.0 - max(0.0, phases[continuousPhaseName]);
    volScalarField base =
        1.0 - min(alphaS / (alphasMax + ROOTVSMALL), 1.0 - 1.0e-6);

    muMix = muC * pow(base, -1.55);

    Info << "min muMix " << min(muMix).value() << " max muMix "
         << max(muMix).value() << endl;

    // Compute mass-weighted mixture velocity
    U = 0.0 * phases_[0].U();
    forAll(phases_, phasei)
    {
        phaseModel& phase = phases_[phasei];
        U += phase * phase.rho() * phase.U() / rho;
    }

    // Initialize dmdts
    forAll(dmdts_, i)
    {
        dmdts_.set(i,
                   new volScalarField::Internal(
                       IOobject("dmdt",
                                runTime.name(),
                                mesh,
                                IOobject::NO_READ,
                                IOobject::NO_WRITE),
                       mesh,
                       dimensionedScalar(dimDensity / dimTime, 0)));
    }

    forAll(clouds, cI)
    {
        clouds[cI].info();
    }

    if (pimple.dict().lookupOrDefault<bool>("hydrostaticInitialisation", false))
    {

        const Time& runTime = mesh().time();
        scalar startTime_ = runTime.startTime().value();
        scalar deltaT = runTime.deltaT().value();

        // set small value for deltaT to evolve particles
        const_cast<Time&>(runTime).setDeltaT(1.e-5 * deltaT);

        // increase time iterator
        const_cast<Time&>(runTime)++;

        // evolve particle cloud
        clouds.evolve();

        // restore startTime
        const_cast<Time&>(runTime).setTime(startTime_, startTime_);

        // restore deltaT
        const_cast<Time&>(runTime).setDeltaT(deltaT);

        // write everything (including lagrangian)
        const_cast<Time&>(runTime).writeNow();
    }

    correctCoNum();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::OpenPDAC::~OpenPDAC() {}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::OpenPDAC::preSolve()
{
    correctCoNum();

    // Store divU from the previous mesh so that it can be
    // mapped and used in correctPhi to ensure the corrected phi
    // has the same divergence.
    if (correctPhi || mesh.topoChanging())
    {
        // Construct and register divU for mapping
        divU = new volScalarField(
            "divU0", fvc::div(fvc::absolute(phi, movingPhases[0].U())));
    }

    if (mesh.topoChanging() || mesh.distributing())
    {
        clouds.storeGlobalPositions();
    }

    fvModels().preUpdateMesh();

    isMeshChanging_ = mesh_.update();

    if (isMeshChanging_)
    {
        Info << "MESH UPDATED. Disabling PIMPLE loop early exit for this "
             << "time step." << endl;
    }

    pimpleIter = 0;
    forceFinalPimpleIter_ = false;
    ratioFirstCheck = false;

    if (lowPressureTimestepCorrection)
    {
        // Initialize p_ratio to find the minimum value during PIMPLE iterations
        p_ratio = 1.0;
    }

    if (tempOscillationControl_)
    {
        // Initialize temp oscillation stats
        forAll(minTempMin_, i)
        {
            minTempMin_[i] = GREAT;
            maxTempMin_[i] = -GREAT;
            sumTempMin_[i] = 0.0;
            pimpleTempIterCount_[i] = 0;
        }
    }
}


void Foam::solvers::OpenPDAC::prePredictor()
{

    pimpleIter++;

    if (forceFinalPimpleIter_ && !pimple.finalIter())
    {
        Info << "PIMPLE iter: " << pimpleIter
             << " -> Bypassing (waiting for final iteration)." << endl;
        return; // Non fare nulla in questa iterazione
    }

    if (pimple.thermophysics() || pimple.flow())
    {
        alphaControls.correct(CoNum);

        fluid_.solve(alphaControls, rAs, momentumTransferSystem_);

        fluid_.correct();

        // Reset dmdts to zero
        forAll(dmdts_, i)
        {
            dmdts_[i] = dimensionedScalar(dimDensity / dimTime, 0);
        }

        fluid_.correctContinuityError(dmdts_);
    }
}


void Foam::solvers::OpenPDAC::momentumTransportPredictor()
{

    if (forceFinalPimpleIter_ && !pimple.finalIter())
    {
        return; // Non fare nulla in questa iterazione
    }

    fluid_.predictMomentumTransport();
}


void Foam::solvers::OpenPDAC::thermophysicalTransportPredictor()
{
    // Moved inside the nEnergyCorrectors loop in thermophysicalPredictor()
    // fluid_.predictThermophysicalTransport();
}


void Foam::solvers::OpenPDAC::momentumTransportCorrector()
{

    if (forceFinalPimpleIter_ && !pimple.finalIter())
    {
        return; // Non fare nulla in questa iterazione
    }

    fluid_.correctMomentumTransport();
}


void Foam::solvers::OpenPDAC::thermophysicalTransportCorrector()
{

    if (forceFinalPimpleIter_ && !pimple.finalIter())
    {
        return; // Non fare nulla in questa iterazione
    }

    fluid_.correctThermophysicalTransport();
}


void Foam::solvers::OpenPDAC::postSolve()
{
    if (tempOscillationControl_)
    {
        // Reset for next timestep
        tempOscillationFactor_ = 1.0;

        // Calculate oscillation metric for each thermal phase and find the
        // minimum ratio
        forAll(fluid_.thermalPhases(), i)
        {
            const phaseModel& phase = fluid_.thermalPhases()[i];
            if (pimpleTempIterCount_[i] > 1)
            {
                scalar T_min_min = minTempMin_[i];
                scalar T_max_min = maxTempMin_[i];
                scalar T_avg_min = sumTempMin_[i] / pimpleTempIterCount_[i];

                if (T_avg_min > ROOTVSMALL)
                {
                    scalar oscillation = (T_max_min - T_min_min) / T_avg_min;
                    Info << "Phase " << phase.name()
                         << " minTemp oscillation: " << oscillation
                         << " (min: " << T_min_min << ", max: " << T_max_min
                         << ", avg: " << T_avg_min << ")" << endl;

                    if (oscillation > maxTempOscillation_)
                    {
                        scalar ratio =
                            max(0.01, maxTempOscillation_ / oscillation);
                        tempOscillationFactor_ =
                            min(tempOscillationFactor_, ratio);
                    }
                }
            }
        }
    }

    divU.clear();

    if (mesh.topoChanged())
    {
        muMix.primitiveFieldRef() = 0.0; // Pulisce i valori vecchi
    }

    volScalarField alphasMax =
        max(fluid_.alfasMax(), dimensionedScalar(dimless, 0.01));
    const word& continuousPhaseName = fluid.continuousPhaseName();

    volScalarField alphaS =
        1.0 - max(0.0, min(1.0, phases[continuousPhaseName]));
    volScalarField base = max(0.01, 1.0 - alphaS / alphasMax);

    const volScalarField& muC = phases[continuousPhaseName].fluidThermo().mu();
    volScalarField muTemp = muC * pow(base, -1.55);

    muMix = max(muTemp, muC);

    muMix.correctBoundaryConditions();

    rho = fluid_.rho();

    U *= 0.0;
    forAll(phases_, phasei)
    {
        phaseModel& phase = phases_[phasei];
        U += phase * phase.rho() * phase.U() / rho;
    }

    clouds.evolve();
}


// ************************************************************************* //
