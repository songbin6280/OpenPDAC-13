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

Class
    Foam::kineticTheoryModels::radialModel

SourceFiles
    radialModel.C

\*---------------------------------------------------------------------------*/

#ifndef radialModel_H
#define radialModel_H

#include "dictionary.H"
#include "volFields.H"
#include "dimensionedTypes.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phaseSystem;
class phaseModel;

/*---------------------------------------------------------------------------*\
                         Class radialModel Declaration
\*---------------------------------------------------------------------------*/

class radialModel
{
    // Private member functions

    //- Read the coefficients from coeffDict
    virtual bool readCoeffs(const dictionary& coeffDict) { return true; }

protected:
    //- Reference to the phase system
    const phaseSystem& fluid_;

public:
    //- Runtime type information
    TypeName("radialModel");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable(autoPtr,
                                 radialModel,
                                 dictionary,
                                 (const dictionary& dict,
                                  const phaseSystem& fluid),
                                 (dict, fluid));


    // Constructors

    //- Construct from the coefficients dictionary
    radialModel(const dictionary& coeffDict, const phaseSystem& fluid);

    //- Disallow default bitwise copy construction
    radialModel(const radialModel&) = delete;


    // Selectors

    static autoPtr<radialModel> New(const dictionary& dict,
                                    const phaseSystem& fluid);


    //- Destructor
    virtual ~radialModel();


    // Member Functions
    bool read(const dictionary& dict);

    virtual PtrList<volScalarField>
    g0(const phaseModel& phase,
       const phaseModel& continuousPhase,
       const dimensionedScalar& alphaMinFriction,
       const volScalarField& alphasMax) const = 0;

    virtual PtrList<volScalarField>
    g0prime(const phaseModel& phase,
            const phaseModel& continuousPhase,
            const dimensionedScalar& alphaMinFriction,
            const volScalarField& alphasMax) const = 0;


    // Member Operators
    void operator=(const radialModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
