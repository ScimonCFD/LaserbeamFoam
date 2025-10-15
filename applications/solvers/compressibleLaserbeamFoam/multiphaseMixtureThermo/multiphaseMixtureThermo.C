/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "multiphaseMixtureThermo.H"
#include "alphaContactAngleFvPatchScalarField.H"
#include "Time.H"
#include "subCycle.H"
#include "MULES.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "fvcFlux.H"
#include "fvcMeshPhi.H"
#include "surfaceInterpolate.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multiphaseMixtureThermo, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::multiphaseMixtureThermo::calcAlphas()
{
    scalar level = 0.0;
    alphas_ == 0.0;

    for (const phaseModel& phase : phases_)
    {
        alphas_ += level * phase;
        level += 1.0;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiphaseMixtureThermo::multiphaseMixtureThermo
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    psiThermo(U.mesh(), word::null),
    phases_(lookup("phases"), phaseModel::iNew(p_, T_)),

    mesh_(U.mesh()),
    U_(U),
    phi_(phi),

    rhoPhi_
    (
        IOobject
        (
            "rhoPhi",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/dimTime, Zero)
    ),

    alphas_
    (
        IOobject
        (
            "alphas",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE,
            IOobject::REGISTER
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    ),

    sigmas_(lookup("sigmas")),
    dimSigma_(1, 0, -2, 0, 0),
    dsigmadT_(lookup("dsigmadT")),
    dimdSigmadT_(1, 0, -2, -1, 0),
    boils_(lookup("boils")),
    dimBoil_(0, 0, 0, 1, 0),
    LatentHeatGass_(lookup("LatentHeatGas")),
    dimLatentHeatGas_(0, 2, -2, 0, 0),
    cAlphas_(lookup("interfaceCompression")),
    dAlphas_(lookup("interfaceDiffusion")),
    dimdiff_(0, 2, -1, 0, 0),
    deltaN_
    (
        "deltaN",
        1e-8/cbrt(average(mesh_.V()))
    ),
    phaseChangeSources_(0),
    vDotPhaseChange_
    (
        IOobject
        (
            "vDotPhaseChange",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/dimTime, Zero)
    )
{
    rhoPhi_.setOriented();
    calcAlphas();
    alphas_.write();
    correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::multiphaseMixtureThermo::correct()
{
    for (phaseModel& phase : phases_)
    {
        phase.correct();
    }

    auto phasei = phases_.cbegin();

    psi_ = phasei()*phasei().thermo().psi();
    mu_ = phasei()*phasei().thermo().mu();
    alpha_ = phasei()*phasei().thermo().alpha();

    for (++phasei; phasei != phases_.cend(); ++phasei)
    {
        psi_ += phasei()*phasei().thermo().psi();
        mu_ += phasei()*phasei().thermo().mu();
        alpha_ += phasei()*phasei().thermo().alpha();
    }
}


void Foam::multiphaseMixtureThermo::correctRho(const volScalarField& dp)
{
    for (phaseModel& phase : phases_)
    {
        phase.thermo().rho() += phase.thermo().psi()*dp;
    }
}


Foam::word Foam::multiphaseMixtureThermo::thermoName() const
{
    auto phasei = phases_.cbegin();

    word name = phasei().thermo().thermoName();

    for (++phasei; phasei != phases_.cend(); ++phasei)
    {
        name += ',' + phasei().thermo().thermoName();
    }

    return name;
}


bool Foam::multiphaseMixtureThermo::incompressible() const
{
    for (const phaseModel& phase : phases_)
    {
        if (!phase.thermo().incompressible())
        {
            return false;
        }
    }

    return true;
}


bool Foam::multiphaseMixtureThermo::isochoric() const
{
    for (const phaseModel& phase : phases_)
    {
        if (!phase.thermo().isochoric())
        {
            return false;
        }
    }

    return true;
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseMixtureThermo::he
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    auto phasei = phases_.cbegin();

    tmp<volScalarField> the(phasei()*phasei().thermo().he(p, T));

    for (++phasei; phasei != phases_.cend(); ++phasei)
    {
        the.ref() += phasei()*phasei().thermo().he(p, T);
    }

    return the;
}


Foam::tmp<Foam::scalarField> Foam::multiphaseMixtureThermo::he
(
    const scalarField& p,
    const scalarField& T,
    const labelList& cells
) const
{
    auto phasei = phases_.cbegin();

    tmp<scalarField> the
    (
        scalarField(phasei(), cells)*phasei().thermo().he(p, T, cells)
    );

    for (++phasei; phasei != phases_.cend(); ++phasei)
    {
        the.ref() +=
            scalarField(phasei(), cells)*phasei().thermo().he(p, T, cells);
    }

    return the;
}


Foam::tmp<Foam::scalarField> Foam::multiphaseMixtureThermo::he
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    auto phasei = phases_.cbegin();

    tmp<scalarField> the
    (
        phasei().boundaryField()[patchi]*phasei().thermo().he(p, T, patchi)
    );

    for (++phasei; phasei != phases_.cend(); ++phasei)
    {
        the.ref() +=
            phasei().boundaryField()[patchi]*phasei().thermo().he(p, T, patchi);
    }

    return the;
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseMixtureThermo::hc() const
{
    auto phasei = phases_.cbegin();

    tmp<volScalarField> thc(phasei()*phasei().thermo().hc());

    for (++phasei; phasei != phases_.cend(); ++phasei)
    {
        thc.ref() += phasei()*phasei().thermo().hc();
    }

    return thc;
}


Foam::tmp<Foam::scalarField> Foam::multiphaseMixtureThermo::THE
(
    const scalarField& h,
    const scalarField& p,
    const scalarField& T0,
    const labelList& cells
) const
{
    NotImplemented;
    return T0;
}


Foam::tmp<Foam::scalarField> Foam::multiphaseMixtureThermo::THE
(
    const scalarField& h,
    const scalarField& p,
    const scalarField& T0,
    const label patchi
) const
{
    NotImplemented;
    return T0;
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseMixtureThermo::rho() const
{
    auto phasei = phases_.cbegin();

    tmp<volScalarField> trho(phasei()*phasei().thermo().rho());

    for (++phasei; phasei != phases_.cend(); ++phasei)
    {
        trho.ref() += phasei()*phasei().thermo().rho();
    }

    return trho;
}


Foam::tmp<Foam::scalarField> Foam::multiphaseMixtureThermo::rho
(
    const label patchi
) const
{
    auto phasei = phases_.cbegin();

    tmp<scalarField> trho
    (
        phasei().boundaryField()[patchi]*phasei().thermo().rho(patchi)
    );

    for (++phasei; phasei != phases_.cend(); ++phasei)
    {
        trho.ref() +=
            phasei().boundaryField()[patchi]*phasei().thermo().rho(patchi);
    }

    return trho;
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseMixtureThermo::Cp() const
{
    auto phasei = phases_.cbegin();

    tmp<volScalarField> tCp(phasei()*phasei().thermo().Cp());

    for (++phasei; phasei != phases_.cend(); ++phasei)
    {
        tCp.ref() += phasei()*phasei().thermo().Cp();
    }

    return tCp;
}


Foam::tmp<Foam::scalarField> Foam::multiphaseMixtureThermo::Cp
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    auto phasei = phases_.cbegin();

    tmp<scalarField> tCp
    (
        phasei().boundaryField()[patchi]*phasei().thermo().Cp(p, T, patchi)
    );

    for (++phasei; phasei != phases_.cend(); ++phasei)
    {
        tCp.ref() +=
            phasei().boundaryField()[patchi]*phasei().thermo().Cp(p, T, patchi);
    }

    return tCp;
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseMixtureThermo::Cv() const
{
    auto phasei = phases_.cbegin();

    tmp<volScalarField> tCv(phasei()*phasei().thermo().Cv());

    for (++phasei; phasei != phases_.cend(); ++phasei)
    {
        tCv.ref() += phasei()*phasei().thermo().Cv();
    }

    return tCv;
}


Foam::tmp<Foam::scalarField> Foam::multiphaseMixtureThermo::Cv
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    auto phasei = phases_.cbegin();

    tmp<scalarField> tCv
    (
        phasei().boundaryField()[patchi]*phasei().thermo().Cv(p, T, patchi)
    );

    for (++phasei; phasei != phases_.cend(); ++phasei)
    {
        tCv.ref() +=
            phasei().boundaryField()[patchi]*phasei().thermo().Cv(p, T, patchi);
    }

    return tCv;
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseMixtureThermo::gamma() const
{
    auto phasei = phases_.cbegin();

    tmp<volScalarField> tgamma(phasei()*phasei().thermo().gamma());

    for (++phasei; phasei != phases_.cend(); ++phasei)
    {
        tgamma.ref() += phasei()*phasei().thermo().gamma();
    }

    return tgamma;
}


Foam::tmp<Foam::scalarField> Foam::multiphaseMixtureThermo::gamma
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    auto phasei = phases_.cbegin();

    tmp<scalarField> tgamma
    (
        phasei().boundaryField()[patchi]*phasei().thermo().gamma(p, T, patchi)
    );

    for (++phasei; phasei != phases_.cend(); ++phasei)
    {
        tgamma.ref() +=
            phasei().boundaryField()[patchi]
           *phasei().thermo().gamma(p, T, patchi);
    }

    return tgamma;
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseMixtureThermo::Cpv() const
{
    auto phasei = phases_.cbegin();

    tmp<volScalarField> tCpv(phasei()*phasei().thermo().Cpv());

    for (++phasei; phasei != phases_.cend(); ++phasei)
    {
        tCpv.ref() += phasei()*phasei().thermo().Cpv();
    }

    return tCpv;
}


Foam::tmp<Foam::scalarField> Foam::multiphaseMixtureThermo::Cpv
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    auto phasei = phases_.cbegin();

    tmp<scalarField> tCpv
    (
        phasei().boundaryField()[patchi]*phasei().thermo().Cpv(p, T, patchi)
    );

    for (++phasei; phasei != phases_.cend(); ++phasei)
    {
        tCpv.ref() +=
            phasei().boundaryField()[patchi]
           *phasei().thermo().Cpv(p, T, patchi);
    }

    return tCpv;
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseMixtureThermo::CpByCpv() const
{
    auto phasei = phases_.cbegin();

    tmp<volScalarField> tCpByCpv(phasei()*phasei().thermo().CpByCpv());

    for (++phasei; phasei != phases_.cend(); ++phasei)
    {
        tCpByCpv.ref() += phasei()*phasei().thermo().CpByCpv();
    }

    return tCpByCpv;
}


Foam::tmp<Foam::scalarField> Foam::multiphaseMixtureThermo::CpByCpv
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    auto phasei = phases_.cbegin();

    tmp<scalarField> tCpByCpv
    (
        phasei().boundaryField()[patchi]*phasei().thermo().CpByCpv(p, T, patchi)
    );

    for (++phasei; phasei != phases_.cend(); ++phasei)
    {
        tCpByCpv.ref() +=
            phasei().boundaryField()[patchi]
           *phasei().thermo().CpByCpv(p, T, patchi);
    }

    return tCpByCpv;
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseMixtureThermo::W() const
{
    auto phasei = phases_.cbegin();

    tmp<volScalarField> tW(phasei()*phasei().thermo().W());

    for (++phasei; phasei != phases_.cend(); ++phasei)
    {
        tW.ref() += phasei()*phasei().thermo().W();
    }

    return tW;
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseMixtureThermo::nu() const
{
    return mu()/rho();
}


Foam::tmp<Foam::scalarField> Foam::multiphaseMixtureThermo::nu
(
    const label patchi
) const
{
    return mu(patchi)/rho(patchi);
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseMixtureThermo::kappa() const
{
    auto phasei = phases_.cbegin();

    tmp<volScalarField> tkappa(phasei()*phasei().thermo().kappa());

    for (++phasei; phasei != phases_.cend(); ++phasei)
    {
        tkappa.ref() += phasei()*phasei().thermo().kappa();
    }

    return tkappa;
}


Foam::tmp<Foam::scalarField> Foam::multiphaseMixtureThermo::kappa
(
    const label patchi
) const
{
    auto phasei = phases_.cbegin();

    tmp<scalarField> tkappa
    (
        phasei().boundaryField()[patchi]*phasei().thermo().kappa(patchi)
    );

    for (++phasei; phasei != phases_.cend(); ++phasei)
    {
        tkappa.ref() +=
            phasei().boundaryField()[patchi]*phasei().thermo().kappa(patchi);
    }

    return tkappa;
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseMixtureThermo::alphahe() const
{
    PtrDictionary<phaseModel>::const_iterator phasei = phases_.begin();

    tmp<volScalarField> talphaEff(phasei()*phasei().thermo().alphahe());

    for (++phasei; phasei != phases_.end(); ++phasei)
    {
        talphaEff.ref() += phasei()*phasei().thermo().alphahe();
    }

    return talphaEff;
}


Foam::tmp<Foam::scalarField> Foam::multiphaseMixtureThermo::alphahe
(
    const label patchi
) const
{
    PtrDictionary<phaseModel>::const_iterator phasei = phases_.begin();

    tmp<scalarField> talphaEff
    (
        phasei().boundaryField()[patchi]
       *phasei().thermo().alphahe(patchi)
    );

    for (++phasei; phasei != phases_.end(); ++phasei)
    {
        talphaEff.ref() +=
            phasei().boundaryField()[patchi]
           *phasei().thermo().alphahe(patchi);
    }

    return talphaEff;
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseMixtureThermo::kappaEff
(
    const volScalarField& alphat
) const
{
    auto phasei = phases_.cbegin();

    tmp<volScalarField> tkappaEff(phasei()*phasei().thermo().kappaEff(alphat));

    for (++phasei; phasei != phases_.cend(); ++phasei)
    {
        tkappaEff.ref() += phasei()*phasei().thermo().kappaEff(alphat);
    }

    return tkappaEff;
}


Foam::tmp<Foam::scalarField> Foam::multiphaseMixtureThermo::kappaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    auto phasei = phases_.cbegin();

    tmp<scalarField> tkappaEff
    (
        phasei().boundaryField()[patchi]
       *phasei().thermo().kappaEff(alphat, patchi)
    );

    for (++phasei; phasei != phases_.cend(); ++phasei)
    {
        tkappaEff.ref() +=
            phasei().boundaryField()[patchi]
           *phasei().thermo().kappaEff(alphat, patchi);
    }

    return tkappaEff;
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseMixtureThermo::alphaEff
(
    const volScalarField& alphat
) const
{
    auto phasei = phases_.cbegin();

    tmp<volScalarField> talphaEff(phasei()*phasei().thermo().alphaEff(alphat));

    for (++phasei; phasei != phases_.cend(); ++phasei)
    {
        talphaEff.ref() += phasei()*phasei().thermo().alphaEff(alphat);
    }

    return talphaEff;
}


Foam::tmp<Foam::scalarField> Foam::multiphaseMixtureThermo::alphaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    auto phasei = phases_.cbegin();

    tmp<scalarField> talphaEff
    (
        phasei().boundaryField()[patchi]
       *phasei().thermo().alphaEff(alphat, patchi)
    );

    for (++phasei; phasei != phases_.cend(); ++phasei)
    {
        talphaEff.ref() +=
            phasei().boundaryField()[patchi]
           *phasei().thermo().alphaEff(alphat, patchi);
    }

    return talphaEff;
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseMixtureThermo::rCv() const
{
    auto phasei = phases_.cbegin();

    tmp<volScalarField> trCv(phasei()/phasei().thermo().Cv());

    for (++phasei; phasei != phases_.cend(); ++phasei)
    {
        trCv.ref() += phasei()/phasei().thermo().Cv();
    }

    return trCv;
}


Foam::tmp<Foam::surfaceScalarField>
Foam::multiphaseMixtureThermo::surfaceTensionForce() const
{
    tmp<surfaceScalarField> tstf
    (
        new surfaceScalarField
        (
            IOobject
            (
                "surfaceTensionForce",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar(dimensionSet(1, -2, -2, 0, 0), Zero)
        )
    );

    surfaceScalarField& stf = tstf.ref();
    stf.setOriented();

    forAllConstIters(phases_, phase1)
    {
        const phaseModel& alpha1 = *phase1;

        auto phase2 = phase1;

        for (++phase2; phase2 != phases_.cend(); ++phase2)
        {
            const phaseModel& alpha2 = *phase2;

            auto sigma = sigmas_.cfind(interfacePair(alpha1, alpha2));

            if (!sigma.good())
            {
                FatalErrorInFunction
                    << "Cannot find interface " << interfacePair(alpha1, alpha2)
                    << " in list of sigma values"
                    << exit(FatalError);
            }

            stf += dimensionedScalar("sigma", dimSigma_, *sigma)
               *fvc::interpolate(K(alpha1, alpha2))*
                (
                    fvc::interpolate(alpha2)*fvc::snGrad(alpha1)
                  - fvc::interpolate(alpha1)*fvc::snGrad(alpha2)
                );
        }
    }

    return tstf;
}


Foam::tmp<Foam::volVectorField>
Foam::multiphaseMixtureThermo::MarangoniForce(
const volScalarField& Temperature
) const
{
    tmp<volVectorField> tsMf
    (
        new volVectorField
        (
            IOobject
            (
                "MarangoniForce",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedVector
            (
                "MarangoniForce",
                dimensionSet(1, -2, -2, 0, 0),
                vector::zero
            )
        )
    );

    volVectorField& sMf = tsMf.ref();
    // sMf.setOriented();

    forAllConstIters(phases_, iter1)
    {
        const phaseModel& alpha1 = iter1();

        auto iter2 = iter1;

        for (++iter2; iter2 != phases_.cend(); ++iter2)
        {
            const phaseModel& alpha2 = iter2();

            auto dsigmadT = dsigmadT_.cfind(interfacePair(alpha1, alpha2));

            if (!dsigmadT.found())
            {
                FatalErrorInFunction
                    << "Cannot find interface " << interfacePair(alpha1, alpha2)
                    << " in list of dsigmadT values"
                    << exit(FatalError);
            }

            // dsigmadTTable::const_iterator dsigmadT =
            //     dsigmadT_.find(interfacePair(alpha1, alpha2));

            // if (dsigmadT == dsigmadT_.end())
            // {
            //     FatalErrorInFunction
            //         << "Cannot find interface " << interfacePair(alpha1, alpha2)
            //         << " in list of dsigmadT values"
            //         << exit(FatalError);
            // }

                // Cell gradient of alpha
        const volVectorField gradAlpha
        (
            alpha2*fvc::grad(alpha1) - alpha1*fvc::grad(alpha2)
        );

        const volVectorField nHatM(gradAlpha/(mag(gradAlpha) + deltaN_));

        const volVectorField gradT(fvc::grad(Temperature));



        sMf -=
            dimensionedScalar("dsigmadT", dimdSigmadT_, dsigmadT())
           *(gradT-(nHatM*(nHatM & gradT)))
           *mag(gradAlpha);
        }
    }

    return tsMf;
}



Foam::tmp<Foam::volScalarField> Foam::multiphaseMixtureThermo::solve
(
    volScalarField *massdotterm
)
{
    tmp<volScalarField> tPCR
    (
        new volScalarField
        (
            IOobject
            (
                "PhaseChangeRate",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar
            (
                "PhaseChangeRate",
                dimensionSet(0, 0, -1, 0, 0),
                0.0
            )
        )
    );

    volScalarField& PCR = tPCR.ref();
    *massdotterm *= 0.0;

    const Time& runTime = mesh_.time();
    const dictionary& alphaControls = mesh_.solverDict("alpha");
    label nAlphaSubCycles(alphaControls.get<label>("nAlphaSubCycles"));

    volScalarField& alpha = phases_.first();


    // CALCULATE PHASE CHANGE SOURCES ONCE BEFORE SUB-CYCLING

    dimensionedScalar outerDeltaT = runTime.deltaT();
    
    PtrList<volScalarField> phaseChangeSources(phases_.size());
    volScalarField vDotPhaseChange
    (
        IOobject
        (
            "vDotPhaseChange",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("vDot", dimless/dimTime, 0.0)
    );
    
    calculatePhaseChangeSources
    (
        phaseChangeSources,
        vDotPhaseChange,
        massdotterm,
        outerDeltaT
    );
    
   
    PCR = vDotPhaseChange;
    
    Info << "Phase change volumetric dilation rate:" << endl
         << "  Mean: " << PCR.weightedAverage(mesh_.V()).value() << " 1/s" << endl
         << "  Max:  " << max(mag(PCR)).value() << " 1/s" << endl;
    

    scalar totalMassBefore = 0.0;
    for (const phaseModel& phase : phases_)
    {
        volScalarField phaseMass = phase * phase.thermo().rho();
        totalMassBefore += gSum(phaseMass.primitiveField() * mesh_.V().field());
    }
    Info << "Total mass before alpha solve: " << totalMassBefore << " kg" << endl;
    



        if (phaseChangeSources_.size() != phases_.size())
    {
        phaseChangeSources_.resize(phases_.size());
    }
    
    // Calculate and STORE sources
    calculatePhaseChangeSources
    (
        phaseChangeSources_,
        vDotPhaseChange_,
        massdotterm,
        outerDeltaT
    );



    

    if (nAlphaSubCycles > 1)
    {
        surfaceScalarField rhoPhiSum(0.0*rhoPhi_);

        for
        (
            subCycle<volScalarField> alphaSubCycle(alpha, nAlphaSubCycles);
            !(++alphaSubCycle).end();
        )
        {
            solveAlphasWithSources(phaseChangeSources, vDotPhaseChange);
            rhoPhiSum += (runTime.deltaT()/outerDeltaT)*rhoPhi_;
        }

        rhoPhi_ = rhoPhiSum;
    }
    else
    {
        solveAlphasWithSources(phaseChangeSources, vDotPhaseChange);
    }
    
    // Calculate total mass after sub-cycling for verification
    scalar totalMassAfter = 0.0;
    for (const phaseModel& phase : phases_)
    {
        volScalarField phaseMass = phase * phase.thermo().rho();
        totalMassAfter += gSum(phaseMass.primitiveField() * mesh_.V().field());
    }
    Info << "Total mass after alpha solve:  " << totalMassAfter << " kg" << endl;
    Info << "Mass change: " << (totalMassAfter - totalMassBefore) << " kg" << endl;
    Info << "Relative mass error: " << (totalMassAfter - totalMassBefore)/totalMassBefore * 100.0 << " %" << endl;

    return tPCR;
}

void Foam::multiphaseMixtureThermo::calculatePhaseChangeSources
(
    PtrList<volScalarField>& phaseChangeSources,
    volScalarField& vDotPhaseChange,
    volScalarField* massdotterm,
    const dimensionedScalar& deltaT
)
{
    // Initialize
    vDotPhaseChange *= 0.0;
    *massdotterm *= 0.0;
    
    // Create condensate field to know what is vapour and not......
    volScalarField condensate
    (
        IOobject("condensate", mesh_.time().timeName(), mesh_),
        mesh_,
        dimensionedScalar("condensate", dimless, 1.0)
    );

    for (const phaseModel& phase : phases_)
    {
        if (phase.name().find("vapour") != string::npos || phase.name() == "air")
        {
            condensate -= phase;
        }
    }
    condensate = max(min(condensate, 1.0), 0.0);

    // Read properties
    IOdictionary phasedictionary
    (
        IOobject
        (
            "thermophysicalProperties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ_IF_MODIFIED
        )
    );

    dimensionedScalar P0("P0", dimPressure, phasedictionary);
    dimensionedScalar int_thickness("interface_thickness", dimLength, phasedictionary);
    dimensionedScalar gasconstant("gasconstant", dimensionSet(1, 2, -2, -1, -1), 8.314);
    
    
    dimensionedScalar maxrate("maxrate", dimless/dimTime, 0.05/deltaT.value());

    // Initialize source storage
    label phasei = 0;
    for (phaseModel& alpha : phases_)
    {
        phaseChangeSources.set
        (
            phasei,
            new volScalarField
            (
                IOobject
                (
                    "phaseChangeSource." + alpha.name(),
                    mesh_.time().timeName(),
                    mesh_
                ),
                mesh_,
                dimensionedScalar("zero", dimless/dimTime, 0.0)
            )
        );
        ++phasei;
    }

    // Create phase index map
    HashTable<label> phaseIndex;
    phasei = 0;
    for (const phaseModel& phase : phases_)
    {
        phaseIndex.insert(phase.name(), phasei++);
    }

    HashTable<bool> processedPairs;

    // Process phase pairs for vaporization/condensation
    for (phaseModel& phase1 : phases_)
    {
        for (phaseModel& phase2 : phases_)
        {
            if (&phase2 == &phase1) continue;
            
            bool isLiquidVapor = false;
            phaseModel* liquidPhase = nullptr;
            phaseModel* vaporPhase = nullptr;
            
            if (phase2.name() == (phase1.name() + "vapour"))
            {
                isLiquidVapor = true;
                liquidPhase = &phase1;
                vaporPhase = &phase2;
            }
            else if (phase1.name() == (phase2.name() + "vapour"))
            {
                isLiquidVapor = true;
                liquidPhase = &phase2;
                vaporPhase = &phase1;
            }
            
            if (!isLiquidVapor) continue;
            
            word pairName = liquidPhase->name() + "_" + vaporPhase->name();
            if (processedPairs.found(pairName)) continue;
            processedPairs.insert(pairName, true);
            
            label liqIndex = phaseIndex[liquidPhase->name()];
            label vapIndex = phaseIndex[vaporPhase->name()];
            
            boilTable::const_iterator boilT = boils_.find(interfacePair(*liquidPhase, *vaporPhase));
            if (boilT == boils_.end()) continue;
            
            dimensionedScalar pair_boil_T("pair_boil_T", dimTemperature, boilT());
            
            LatentHeatGasTable::const_iterator LHG = LatentHeatGass_.find(interfacePair(*liquidPhase, *vaporPhase));
            if (LHG == LatentHeatGass_.end()) continue;
            
            dimensionedScalar pair_LHG("pair_LHG", dimEnergy/dimMass, LHG());
            
            Info << "Calculating phase change for pair: " << liquidPhase->name() << " <-> " 
                 << vaporPhase->name() << " Tboil=" << pair_boil_T.value() << "K" << endl;
            
            volScalarField Psat = P0*Foam::exp(
                (liquidPhase->thermo().W()*pair_LHG/(1000.0*pair_boil_T*gasconstant))
                *(1.0 - pair_boil_T/T_)
            );
            
            const volScalarField& rhoLiq = liquidPhase->thermo().rho();
            const volScalarField& rhoVap = vaporPhase->thermo().rho();
            
            volScalarField evapMassRate = 
                max(*liquidPhase, SMALL)*
                Foam::pos(Psat - p_)*
                Foam::sqrt(liquidPhase->thermo().W()/(2000.0*M_PI*gasconstant*T_))*
                (Psat - p_)/int_thickness;
                
            volScalarField condMassRate = 
                max(*vaporPhase, SMALL)*
                max(condensate, 1e-6)*
                Foam::pos(p_ - Psat)*
                Foam::sqrt(vaporPhase->thermo().W()/(2000.0*M_PI*gasconstant*T_))*
                (p_ - Psat)/int_thickness;
            
            evapMassRate = min(evapMassRate, maxrate*rhoLiq*max(*liquidPhase, SMALL));
            condMassRate = min(condMassRate, maxrate*rhoVap*max(*vaporPhase, SMALL));
            
            volScalarField netMassTransfer = evapMassRate - condMassRate;
            
            phaseChangeSources[liqIndex] -= netMassTransfer/rhoLiq;
            phaseChangeSources[vapIndex] += netMassTransfer/rhoVap;
            
            vDotPhaseChange += netMassTransfer*(1.0/rhoVap - 1.0/rhoLiq);
            
            *massdotterm += pair_LHG*netMassTransfer;
        }
    }
}


Foam::tmp<Foam::surfaceVectorField> Foam::multiphaseMixtureThermo::nHatfv
(
    const volScalarField& alpha1,
    const volScalarField& alpha2
) const
{
    /*
    // Cell gradient of alpha
    volVectorField gradAlpha =
        alpha2*fvc::grad(alpha1) - alpha1*fvc::grad(alpha2);

    // Interpolated face-gradient of alpha
    surfaceVectorField gradAlphaf = fvc::interpolate(gradAlpha);
    */

    surfaceVectorField gradAlphaf
    (
        fvc::interpolate(alpha2)*fvc::interpolate(fvc::grad(alpha1))
      - fvc::interpolate(alpha1)*fvc::interpolate(fvc::grad(alpha2))
    );

    // Face unit interface normal
    return gradAlphaf/(mag(gradAlphaf) + deltaN_);
}


Foam::tmp<Foam::surfaceScalarField> Foam::multiphaseMixtureThermo::nHatf
(
    const volScalarField& alpha1,
    const volScalarField& alpha2
) const
{
    // Face unit interface normal flux
    return nHatfv(alpha1, alpha2) & mesh_.Sf();
}


// Correction for the boundary condition on the unit normal nHat on
// walls to produce the correct contact angle.

// The dynamic contact angle is calculated from the component of the
// velocity on the direction of the interface, parallel to the wall.

void Foam::multiphaseMixtureThermo::correctContactAngle
(
    const phaseModel& alpha1,
    const phaseModel& alpha2,
    surfaceVectorField::Boundary& nHatb
) const
{
    const volScalarField::Boundary& gbf
        = alpha1.boundaryField();

    const fvBoundaryMesh& boundary = mesh_.boundary();

    forAll(boundary, patchi)
    {
        if (isA<alphaContactAngleFvPatchScalarField>(gbf[patchi]))
        {
            const alphaContactAngleFvPatchScalarField& acap =
                refCast<const alphaContactAngleFvPatchScalarField>(gbf[patchi]);

            vectorField& nHatPatch = nHatb[patchi];

            vectorField AfHatPatch
            (
                mesh_.Sf().boundaryField()[patchi]
               /mesh_.magSf().boundaryField()[patchi]
            );

            const auto tp =
                acap.thetaProps().cfind(interfacePair(alpha1, alpha2));

            if (!tp.good())
            {
                FatalErrorInFunction
                    << "Cannot find interface " << interfacePair(alpha1, alpha2)
                    << "\n    in table of theta properties for patch "
                    << acap.patch().name()
                    << exit(FatalError);
            }

            const bool matched = (tp.key().first() == alpha1.name());

            const scalar theta0 = degToRad(tp().theta0(matched));
            scalarField theta(boundary[patchi].size(), theta0);

            const scalar uTheta = tp().uTheta();

            // Calculate the dynamic contact angle if required
            if (uTheta > SMALL)
            {
                const scalar thetaA = degToRad(tp().thetaA(matched));
                const scalar thetaR = degToRad(tp().thetaR(matched));

                // Calculated the component of the velocity parallel to the wall
                vectorField Uwall
                (
                    U_.boundaryField()[patchi].patchInternalField()
                  - U_.boundaryField()[patchi]
                );
                Uwall -= (AfHatPatch & Uwall)*AfHatPatch;

                // Find the direction of the interface parallel to the wall
                vectorField nWall
                (
                    nHatPatch - (AfHatPatch & nHatPatch)*AfHatPatch
                );

                // Normalise nWall
                nWall /= (mag(nWall) + SMALL);

                // Calculate Uwall resolved normal to the interface parallel to
                // the interface
                scalarField uwall(nWall & Uwall);

                theta += (thetaA - thetaR)*tanh(uwall/uTheta);
            }


            // Reset nHatPatch to correspond to the contact angle

            scalarField a12(nHatPatch & AfHatPatch);

            scalarField b1(cos(theta));

            scalarField b2(nHatPatch.size());

            forAll(b2, facei)
            {
                b2[facei] = cos(acos(a12[facei]) - theta[facei]);
            }

            scalarField det(1.0 - a12*a12);

            scalarField a((b1 - a12*b2)/det);
            scalarField b((b2 - a12*b1)/det);

            nHatPatch = a*AfHatPatch + b*nHatPatch;

            nHatPatch /= (mag(nHatPatch) + deltaN_.value());
        }
    }
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseMixtureThermo::K
(
    const phaseModel& alpha1,
    const phaseModel& alpha2
) const
{
    tmp<surfaceVectorField> tnHatfv = nHatfv(alpha1, alpha2);

    correctContactAngle(alpha1, alpha2, tnHatfv.ref().boundaryFieldRef());

    // Simple expression for curvature
    return -fvc::div(tnHatfv & mesh_.Sf());
}


Foam::tmp<Foam::volScalarField>
Foam::multiphaseMixtureThermo::nearInterface() const
{
    auto tnearInt = volScalarField::New
    (
        "nearInterface",
        IOobject::NO_REGISTER,
        mesh_,
        dimensionedScalar(dimless, Zero)
    );

    for (const phaseModel& phase : phases_)
    {
        tnearInt.ref() =
            max(tnearInt(), pos0(phase - 0.01)*pos0(0.99 - phase));
    }

    return tnearInt;
}


void Foam::multiphaseMixtureThermo::solveAlphasWithSources
(
    const PtrList<volScalarField>& phaseChangeSources,
    const volScalarField& vDotPhaseChange
)
{
    const word alphaScheme("div(phi,alpha)");
    const word alpharScheme("div(phirb,alpha)");

    // Create local copy of vDot that includes phase change
    volScalarField vDot
    (
        IOobject
        (
            "vDot",
            mesh_.time().timeName(),
            mesh_
        ),
        vDotPhaseChange  
    );

    // Storage 
    PtrList<surfaceScalarField> alphaPhis(phases_.size());
    PtrList<surfaceScalarField> alphaPhiCorrs(phases_.size());
    
    // Initialize phase fluxes as per base compressibleinterfoam solver
    label phasei = 0;
    for (phaseModel& alpha : phases_)
    {
        // Base convective flux
        alphaPhis.set
        (
            phasei,
            new surfaceScalarField
            (
                "alphaPhi." + alpha.name(),
                fvc::flux(phi_, alpha, alphaScheme)
            )
        );
        
        // Corrections (compression/diffusion) - initialize to zero
        alphaPhiCorrs.set
        (
            phasei,
            new surfaceScalarField
            (
                IOobject
                (
                    "alphaPhiCorr." + alpha.name(),
                    mesh_.time().timeName(),
                    mesh_
                ),
                mesh_,
                dimensionedScalar(phi_.dimensions(), Zero)
            )
        );
        
        ++phasei;
    }

    // Create phase index map for speed
    HashTable<label> phaseIndex;
    phasei = 0;
    for (const phaseModel& phase : phases_)
    {
        phaseIndex.insert(phase.name(), phasei++);
    }


    // phasei = 0;
    // for (phaseModel& alpha1 : phases_)
    // {
    //     surfaceScalarField& alphaPhiCorr = alphaPhiCorrs[phasei];
        
    //     for (const phaseModel& alpha2 : phases_)
    //     {
    //         if (&alpha2 == &alpha1) continue;
            
    //         scalarCoeffSymmCTable::const_iterator cAlpha = 
    //             cAlphas_.find(interfacePair(alpha1, alpha2));
            
    //         if (cAlpha != cAlphas_.end() && cAlpha() > 0)
    //         {
    //             surfaceScalarField phic = mag(phi_)/mesh_.magSf();
    //             surfaceScalarField phir = 
    //                 min(cAlpha()*phic, max(phic))*nHatf(alpha1, alpha2);
                
    //             alphaPhiCorr += fvc::flux
    //             (
    //                 -fvc::flux(-phir, alpha2, alpharScheme),
    //                 alpha1,
    //                 alpharScheme
    //             );
    //         }
    //     }
    //     ++phasei;
    // }


    // HashTable<bool> processedDiffusionPairs;

    // for (phaseModel& alpha1 : phases_)
    // {
    //     label idx1 = phaseIndex[alpha1.name()];
        
    //     for (phaseModel& alpha2 : phases_)
    //     {
    //         if (&alpha2 == &alpha1) continue;
            
    //         label idx2 = phaseIndex[alpha2.name()];
            
    //         word pairName;
    //         if (alpha1.name() < alpha2.name())
    //         {
    //             pairName = alpha1.name() + "_" + alpha2.name();
    //         }
    //         else
    //         {
    //             pairName = alpha2.name() + "_" + alpha1.name();
    //         }
            
    //         if (processedDiffusionPairs.found(pairName)) continue;
    //         processedDiffusionPairs.insert(pairName, true);
            
    //         scalarCoeffSymmDTable::const_iterator dAlpha = 
    //             dAlphas_.find(interfacePair(alpha1, alpha2));
            
    //         if (dAlpha != dAlphas_.end() && dAlpha() > 0)
    //         {
    //             dimensionedScalar valdiff("valdiff", dimdiff_, dAlpha());
                
    //             const volScalarField& rho1 = alpha1.thermo().rho();
    //             const volScalarField& rho2 = alpha2.thermo().rho();
                
    //             surfaceScalarField rho1f = fvc::interpolate(rho1);
    //             surfaceScalarField rho2f = fvc::interpolate(rho2);
                
    //             dimensionedScalar rhoSmall("rhoSmall", dimDensity, VSMALL);
                
    //             surfaceScalarField rhoHarmonic = 
    //                 2.0*rho1f*rho2f/max(rho1f + rho2f, rhoSmall);
                
    //             surfaceScalarField massFlux = 
    //                 -valdiff*mesh_.magSf()*rhoHarmonic*fvc::snGrad(alpha1);
                
    //             surfaceScalarField volFlux1 = massFlux/max(rho1f, rhoSmall);
    //             surfaceScalarField volFlux2 = massFlux/max(rho2f, rhoSmall);
                
    //             alphaPhiCorrs[idx1] -= volFlux1;
    //             alphaPhiCorrs[idx2] += volFlux2;
                
    //             // Add diffusion-induced dilation to vDot - check this - does this make sense?!?!?!
    //             volScalarField massFluxVol = fvc::surfaceIntegrate(massFlux);
    //             vDot += massFluxVol*(1.0/rho2 - 1.0/rho1);
    //         }
    //     }
    // }



    

    // MULES limiting on CORRECTIONS - not on total as this is incorrect I think - need to run by Philip!

    volScalarField divU = fvc::div(fvc::absolute(phi_, U_));

phasei = 0;
for (phaseModel& alpha : phases_)
{
    const volScalarField& alphaSource = phaseChangeSources[phasei];
    surfaceScalarField& alphaPhiCorr = alphaPhiCorrs[phasei];
    
    volScalarField Sp
    (
        IOobject("Sp", mesh_.time().timeName(), mesh_),
        mesh_,
        dimensionedScalar("zero", dimless/dimTime, 0.0)
    );
    
    volScalarField Su = alphaSource;  // Explicit source
    
    // CRITICAL: Let MULES see the sources so it can limit them properly
    MULES::limit
    (
        1.0/mesh_.time().deltaT().value(),
        geometricOneField(),
        alpha,
        phi_,
        alphaPhiCorr,
        Sp,
        Su,  // ← MULES will limit this!
        oneField(),  // max bound
        zeroField(), // min bound
        true
    );
    
    ++phasei;
}

    // Ensure sum of limited corrections maintains compatibility
    MULES::limitSum(alphaPhiCorrs);


    // Combine and solve for fucks sake

    rhoPhi_ = dimensionedScalar(rhoPhi_.dimensions(), Zero);

    volScalarField sumAlpha
    (
        IOobject("sumAlpha", mesh_.time().timeName(), mesh_),
        mesh_,
        dimensionedScalar(dimless, Zero)
    );

    phasei = 0;
for (phaseModel& alpha : phases_)
{
        surfaceScalarField& alphaPhi = alphaPhis[phasei];
        surfaceScalarField& alphaPhiCorr = alphaPhiCorrs[phasei];
        
        alphaPhi += alphaPhiCorr;
        
        // NO explicit sources - set to zero!
        volScalarField::Internal Sp
        (
            IOobject("Sp", mesh_.time().timeName(), mesh_),
            mesh_,
            dimensionedScalar(dimless/dimTime, Zero)
        );
        
        volScalarField::Internal Su
        (
            IOobject("Su", mesh_.time().timeName(), mesh_),
            mesh_,
            dimensionedScalar(dimless/dimTime, Zero)  // ZERO sources!
        );

 
        const scalarField& dgdt = alpha.dgdt();
        forAll(dgdt, celli)
        {
            if (dgdt[celli] < 0.0 && alpha[celli] > 0.0)
            {
                Sp[celli] += dgdt[celli]*alpha[celli];
                Su[celli] -= dgdt[celli]*alpha[celli];
            }
            else if (dgdt[celli] > 0.0 && alpha[celli] < 1.0)
            {
                Sp[celli] -= dgdt[celli]*(1.0 - alpha[celli]);
            }
        }
        
        // Cross-phase dgdt coupling
        for (const phaseModel& alpha2 : phases_)
        {
            if (&alpha2 == &alpha) continue;
            
            const scalarField& dgdt2 = alpha2.dgdt();
            forAll(dgdt2, celli)
            {
                if (dgdt2[celli] > 0.0 && alpha2[celli] < 1.0)
                {
                    Sp[celli] -= dgdt2[celli]*(1.0 - alpha2[celli]);
                    Su[celli] += dgdt2[celli]*alpha[celli];
                }
                else if (dgdt2[celli] < 0.0 && alpha2[celli] > 0.0)
                {
                    Sp[celli] += dgdt2[celli]*alpha2[celli];
                }
            }
        }
        
        // Solve with ONLY dgdt, no explicit sources
        MULES::explicitSolve
        (
            geometricOneField(),
            alpha,
            alphaPhi,
            Sp,
            Su
        );
    
    sumAlpha += alpha;
    
    rhoPhi_ += fvc::interpolate(alpha.thermo().rho())*alphaPhi;
    
    Info << alpha.name() << " volume fraction, min, max = "
         << alpha.weightedAverage(mesh_.V()).value()
         << ' ' << min(alpha).value()
         << ' ' << max(alpha).value()
         << endl;
    
    ++phasei;
}
    
    Info << "Phase-sum volume fraction, min, max = "
         << sumAlpha.weightedAverage(mesh_.V()).value()
         << ' ' << min(sumAlpha).value()
         << ' ' << max(sumAlpha).value()
         << endl;
    
    Info << "Volumetric dilation rate, mean, max = "
         << vDot.weightedAverage(mesh_.V()).value()
         << ' ' << max(mag(vDot)).value()
         << endl;

    // Calculate and print total mass
    Info << "Phase pair (liquid + vapor) masses in domain:" << endl;
    
    HashTable<bool> processedPhases;
    
    for (const phaseModel& phase : phases_)
    {
        if (processedPhases.found(phase.name())) continue;
        
        word vaporName = phase.name() + "vapour";
        bool hasVaporPair = false;
        scalar liquidMass = 0.0;
        scalar vaporMass = 0.0;
        
        volScalarField phaseMassLiq = phase * phase.thermo().rho();
        liquidMass = gSum(phaseMassLiq.primitiveField() * mesh_.V().field());
        
        for (const phaseModel& phase2 : phases_)
        {
            if (phase2.name() == vaporName)
            {
                hasVaporPair = true;
                volScalarField phaseMassVap = phase2 * phase2.thermo().rho();
                vaporMass = gSum(phaseMassVap.primitiveField() * mesh_.V().field());
                processedPhases.insert(phase2.name(), true);
                break;
            }
        }
        
        if (hasVaporPair)
        {
            scalar totalPairMass = liquidMass + vaporMass;
            Info << "    " << phase.name() << " + " << vaporName << ": " 
                 << totalPairMass << " kg (liquid: " << liquidMass 
                 << " kg, vapor: " << vaporMass << " kg)" << endl;
            processedPhases.insert(phase.name(), true);
        }
        else if (phase.name() != "air" && phase.name().find("vapour") == string::npos)
        {
            Info << "    " << phase.name() << " (no vapor): " << liquidMass << " kg" << endl;
            processedPhases.insert(phase.name(), true);
        }
        else if (phase.name() == "air")
        {
            Info << "    " << phase.name() << ": " << liquidMass << " kg" << endl;
            processedPhases.insert(phase.name(), true);
        }
    }

    calcAlphas();
}

Foam::tmp<Foam::volScalarField> Foam::multiphaseMixtureThermo::solveAlphas
(
    volScalarField *massdotterm
)
{
    // This is now just a wrapper that calls solve()
    return solve(massdotterm);
}



























// ************************************************************************* //
