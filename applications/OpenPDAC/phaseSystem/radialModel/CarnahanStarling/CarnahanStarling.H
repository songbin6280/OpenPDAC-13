/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    Foam::kineticTheoryModels::radialModels::Lebowitz

Description

SourceFiles
    LebowitzRadial.C

\*---------------------------------------------------------------------------*/

#ifndef CarnahanStarling_H
#define CarnahanStarling_H

#include "radialModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radialModels
{

/*---------------------------------------------------------------------------*\
                       Class Lebowitz Declaration
\*---------------------------------------------------------------------------*/

class CarnahanStarling : public radialModel
{

public:
    //- Runtime type information
    TypeName("CarnahanStarling");


    // Constructors

    //- Construct from the coefficients dictionary
    CarnahanStarling(const dictionary& coeffDict, const phaseSystem& fluid);


    //- Destructor
    virtual ~CarnahanStarling();


    // Member Functions

    PtrList<volScalarField> g0(const phaseModel& phase,
                               const phaseModel& continuousPhase,
                               const dimensionedScalar& alphaMinFriction,
                               const volScalarField& alphasMax) const;

    PtrList<volScalarField> g0prime(const phaseModel& phase,
                                    const phaseModel& continuousPhase,
                                    const dimensionedScalar& alphaMinFriction,
                                    const volScalarField& alphasMax) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radialModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
