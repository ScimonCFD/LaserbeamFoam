/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
    Foam::phaseModel

Description
    Single incompressible phase derived from the phase-fraction.
    Used as part of the multiPhaseMixture for interface-capturing multi-phase
    simulations.

SourceFiles
    phaseModel.C

\*---------------------------------------------------------------------------*/

#ifndef phaseModel_H
#define phaseModel_H

#include "rhoThermo.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class phaseModel Declaration
\*---------------------------------------------------------------------------*/

class phaseModel
:
    public volScalarField
{

    public:

   

        //- Phase state types
        enum class phaseState
        {
            CONDENSED,
            GASEOUS,
            SOLID,
            PLASMA
        };


    // Private data

        word name_;
        const volScalarField& p_;
        const volScalarField& T_;
        autoPtr<rhoThermo> thermo_;
        volScalarField dgdt_;
        phaseState state_;


        // Private Member Functions

        //- Read phase state from dictionary
        phaseState readPhaseState(const dictionary& dict) const;


public:

    // Constructors

        //- Construct from components
        phaseModel
        (
            const word& phaseName,
            const volScalarField& p,
            const volScalarField& T
        );

        //- Return clone
        autoPtr<phaseModel> clone() const;

        //- Return a pointer to a new phaseModel created on freestore
        //  from Istream
        class iNew
        {
            const volScalarField& p_;
            const volScalarField& T_;

        public:

            iNew
            (
                const volScalarField& p,
                const volScalarField& T
            )
            :
                p_(p),
                T_(T)
            {}

            autoPtr<phaseModel> operator()(Istream& is) const
            {
                return autoPtr<phaseModel>::New(word(is), p_, T_);
            }
        };


    // Member Functions

        const word& name() const
        {
            return name_;
        }

        const word& keyword() const
        {
            return name();
        }

        //- Return const-access to phase rhoThermo
        const rhoThermo& thermo() const
        {
            return *thermo_;
        }

        //- Return access to phase rhoThermo
        rhoThermo& thermo()
        {
            return *thermo_;
        }

        //- Return const-access to phase divergence
        const volScalarField& dgdt() const
        {
            return dgdt_;
        }

        //- Return access to phase divergence
        volScalarField& dgdt()
        {
            return dgdt_;
        }


        //- Return the phase state
        phaseState state() const;

        //- Check if phase is condensed
        bool isCondensed() const;

        //- Check if phase is gaseous
        bool isGaseous() const;

        //- Check if phase is solid
        bool isSolid() const;

        //- Check if phase is plasma
        bool isPlasma() const;

        //- Return state as string
        word stateAsWord() const;



        void correct();
};



inline phaseModel::phaseState phaseModel::state() const
{
    return state_;
}

inline bool phaseModel::isCondensed() const
{
    return state_ == phaseState::CONDENSED;
}

inline bool phaseModel::isGaseous() const
{
    return state_ == phaseState::GASEOUS;
}

inline bool phaseModel::isSolid() const
{
    return state_ == phaseState::SOLID;
}

inline bool phaseModel::isPlasma() const
{
    return state_ == phaseState::PLASMA;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
