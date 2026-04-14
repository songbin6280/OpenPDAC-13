/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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

#include "myHydrostaticInitialisation.H"

#include "phaseSystem.H"
#include "fvmLaplacian.H"
#include "fvcDiv.H"
#include "fvcSnGrad.H"
#include "surfaceInterpolate.H"
#include "constrainPressure.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::hydrostaticInitialisation(volScalarField& p_rgh,
                                     volScalarField& ph_rgh,
                                     volScalarField& p,
                                     const uniformDimensionedVectorField& g,
                                     dimensionedScalar& hRef,
                                     const volScalarField& gh,
                                     const surfaceScalarField& ghf,
                                     phaseSystem& fluid,
                                     const dictionary& dict)
{
    // Reference boundary pressure value, initialized to zero.
    dimensionedScalar pBdry(
        "pBdry", dimensionSet(1, -1, -2, 0, 0), scalar(0.0));

    // The entire procedure is executed only if specified in the case
    // dictionary.
    if (dict.lookupOrDefault<bool>("hydrostaticInitialisation", false))
    {
        const fvMesh& mesh = p_rgh.mesh();

        // Create local copies of the mixture density and velocity fields.
        // These are used for hydrostatic calculations without modifying the
        // main solver fields until the end.
        volScalarField rho("rho", fluid.rho());
        volVectorField U("U", fluid.U());

        // Reference height on the boundary, initialized to a very small number.
        dimensionedScalar hBdry("hBdry", dimLength, -VGREAT);

        // Name of the reference patch. Use the name (word), not the label!
        word local_patchName = ""; // Use the name (word), not the label!

        // Hydrostatic initialization is performed only at startup (time=0),
        // not when restarting from a previous time.
        if (!mesh.time().restart())
        {

            volScalarField ph(p);

            // --- Find Reference Patch ---
            // Iterate over all boundary patches to find a suitable patch
            // to use as a reference for the hydrostatic pressure.
            forAll(ph.boundaryField(), bdryID)
            {
                if ((mag(mesh.Sf().boundaryField()[bdryID][0] ^ g).value()
                     <= 1.0e-8 * mag(mesh.Sf().boundaryField()[bdryID][0])
                            * mag(g).value())
                    && (ph.boundaryField()[bdryID].type() == "fixedValue"))
                {
                    // The patch is suitable if:
                    // 1. It is "horizontal" with respect to gravity (cross
                    // product is zero).
                    // 2. It has a 'fixedValue' boundary condition for pressure.

                    const fvPatchScalarField& ph_p = ph.boundaryField()[bdryID];

                    // Additional check: the pressure value must be uniform
                    // across the patch.
                    if (min(ph_p) == max(ph_p))
                    {
                        local_patchName =
                            ph.boundaryField()[bdryID].patch().name();
                        // Store the pressure value and the reference height
                        // coordinate from the found patch.
                        pBdry.value() = min(ph_p);
                        Sout << "pBdry " << pBdry.value() << endl;

                        if (g.component(0).value() != 0.0)
                        {
                            hBdry.value() = ph_p.patch().Cf()[0].component(0);
                        }
                        else if (g.component(1).value() != 0.0)
                        {
                            hBdry.value() = ph_p.patch().Cf()[0].component(1);
                        }
                        else
                        {
                            hBdry.value() = ph_p.patch().Cf()[0].component(2);
                        }
                    }
                }
            }


            // --- Parallel Synchronization ---
            // Ensure all processors use the same reference patch
            // and the same pressure and height values.
            reduce(local_patchName, maxOp<word>());
            const word& patchName = local_patchName;

            Info << "Reference patch " << patchName << endl;


            if (patchName.empty())
            {
                FatalErrorInFunction << "Could not find a suitable reference "
                                        "patch for hydrostatic initialisation."
                                     << exit(FatalError);
            }

            reduce(pBdry, maxOp<dimensionedScalar>());
            reduce(hBdry, maxOp<dimensionedScalar>());

            // Initialize the 'ph' pressure field with the reference value.
            ph = pBdry;

            // --- Iterative Loop for Pressure-Density Coupling ---
            // This loop calculates a hydrostatic pressure profile, accounting
            // for the variation of density with pressure (for compressible
            // fluids).
            for (label i = 0; i < 10; i++)
            {
                // Update pressure and recalculate thermophysical properties
                // (density).
                p = ph;
                fluid.correctThermo();
                rho = fluid.rho();

                ph = pBdry
                   + ((g & mesh.C()) - (-mag(g) * hBdry)) * 0.5
                         * (min(rho) + max(rho));
                // An average density (min+max)/2 is used to stabilize the
                // calculation.
                Info << "min ph " << min(ph).value() << " max ph "
                     << max(ph).value() << endl;
            }

            // Set hRef to 0, it is no longer needed after calculating ph.
            hRef.value() = 0.0;
            Info << "hRef " << hRef.value() << endl;

            // The new hydrostatic pressure profile is used to update the
            // density field.
            p = ph;
            fluid.correctThermo();
            rho = fluid.rho();

            // We initialize the field ph_rgh (pressure minus the hydrostatic
            // part) with the computed pressure and density. This field should
            // be nearly uniform under hydrostatic conditions.
            ph_rgh = ph - rho * gh;

            Info << "Updating p_rgh boundary condition on patch " << patchName
                 << " to be consistent with the reference pressure." << endl;

            // --- Correction of Boundary Condition for p_rgh ---
            // The BC for p_rgh on the reference patch must be made consistent
            // with the calculated ph, rho, and gh values.
            label patchID = mesh.boundaryMesh().findIndex(patchName);

            // we change the fixed value b.c. of ph_rgh at the top face, in
            // order to be consistent with the values of ph, rho and gh
            if (patchID != -1)
            {
                forAll(ph_rgh.boundaryField()[patchID], faceI)
                {
                    ph_rgh.boundaryFieldRef()[patchID][faceI] =
                        ph.boundaryField()[patchID][faceI]
                        - rho.boundaryField()[patchID][faceI]
                              * gh.boundaryField()[patchID][faceI];
                }
            }

            // Calculate a fictitious flux 'phig' due to density gradients
            // to enforce pressure consistency.
            surfaceScalarField rhof("rhof", fvc::interpolate(rho));
            surfaceScalarField phig(
                "phig", -rhof * ghf * fvc::snGrad(rho) * mesh.magSf());

            // Update the pressure BCs (e.g., 'fixedFluxPressure') to ensure
            // flux consistency.
            constrainPressure(ph_rgh, rho, U, phig, rhof);
            p = ph_rgh + rho * gh;

            fluid.correctThermo();
            rho = fluid.rho();

            const fvPatchScalarField& ph_rgh_top =
                ph_rgh.boundaryField()[patchID];


            scalar min_ph_rgh_top(min(ph_rgh_top));
            reduce(min_ph_rgh_top, minOp<scalar>());

            scalar max_ph_rgh_top(max(ph_rgh_top));
            reduce(max_ph_rgh_top, maxOp<scalar>());

            Info << "min ph_rgh top " << min_ph_rgh_top << " max ph_rgh top "
                 << max_ph_rgh_top << endl;

            p = ph_rgh + rho * gh;
            fluid.correctThermo();
            rho = fluid.rho();

            // --- Final Correction Loop for p_rgh ---
            // Solves a Poisson equation for p_rgh to smooth the field and
            // ensure that pressure forces exactly balance buoyancy forces
            // due to density stratification.
            label nCorr(
                dict.lookupOrDefault<label>("nHydrostaticCorrectors", 5));

            for (label i = 0; i < nCorr; i++)
            {
                volScalarField ph_rgh_old("ph_rgh_old", ph_rgh);

                surfaceScalarField rhof("rhof", fvc::interpolate(rho));

                surfaceScalarField phig(
                    "phig", -rhof * ghf * fvc::snGrad(rho) * mesh.magSf());

                // Update the pressure BCs again to ensure flux consistency.
                constrainPressure(ph_rgh, rho, U, phig, rhof);

                // Solve: ∇ ⋅ (rhof ∇p_rgh) = ∇ ⋅ phig
                fvScalarMatrix ph_rghEqn(fvm::laplacian(rhof, ph_rgh)
                                         == fvc::div(phig));

                ph_rghEqn.solve();

                p = ph_rgh + rho * gh;
                fluid.correctThermo();
                rho = fluid.rho();

                Info << "Hydrostatic pressure convergence (max delta) "
                     << max(mag(ph_rgh - ph_rgh_old)).value() << endl;

                Info << "min p " << min(p).value() << " max p "
                     << max(p).value() << endl;
                Info << "min rho " << min(rho).value() << " max rho "
                     << max(rho).value() << endl;
            }

            // Write the initialized fields to usa for the BC.
            ph_rgh.write();
            p.write();

            // Assign the calculated p_rgh field to the main solver field.
            p_rgh = ph_rgh;
        }
    }
    else
    {
        // If hydrostatic initialization is disabled,
        // simply ensure that p_rgh is consistent with p and rho.
        volScalarField rho("rho", fluid.rho());
        // Force p_rgh to be consistent with p
        p_rgh = p - rho * gh;
    }
}


// ************************************************************************* //
