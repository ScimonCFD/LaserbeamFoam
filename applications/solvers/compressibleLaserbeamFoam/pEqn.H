{

     if (correctPhi)
    {
        rAU.ref() = 1.0/UEqn.A();
    }
    else
    {
        rAU = 1.0/UEqn.A();
    }


    // volScalarField rAU("rAU", 1.0/UEqn.A());
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU()));
    volVectorField HbyA(constrainHbyA(rAU()*UEqn.H(), U, p_rgh));
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA)
      + fvc::interpolate(rho*rAU())*fvc::ddtCorr(U, phi)
    );

    surfaceScalarField phig
    (
        (
            mixture.surfaceTensionForce()
          - ghf*fvc::snGrad(rho)
        )*rAUf*mesh.magSf()
    );

    phiHbyA += phig;

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p_rgh, U, phiHbyA, rAUf);

    PtrList<fvScalarMatrix> p_rghEqnComps(mixture.phases().size());

    label phasei = 0;
    forAllConstIters(mixture.phases(), phase)
    {
        const rhoThermo& thermo = phase().thermo();
        const tmp<volScalarField> trho(thermo.rho());
        const volScalarField& rho = trho();

        p_rghEqnComps.set
        (
            phasei,
            (
                fvc::ddt(rho) + thermo.psi()*correction(fvm::ddt(p_rgh))
              + fvc::div(phi, rho) - fvc::Sp(fvc::div(phi), rho)
            ).ptr()
        );

        ++phasei;
    }

    // Cache p_rgh prior to solve for density update
    volScalarField p_rgh_0(p_rgh);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqnIncomp
        (
            fvc::div(phiHbyA)
          - fvm::laplacian(rAUf, p_rgh)
          - vDot
        );

        tmp<fvScalarMatrix> p_rghEqnComp;

        phasei = 0;
        forAllConstIters(mixture.phases(), phase)
        {
            tmp<fvScalarMatrix> hmm
            (
                (max(phase(), scalar(0))/phase().thermo().rho())
               *p_rghEqnComps[phasei]
            );

            if (phasei == 0)
            {
                p_rghEqnComp = hmm;
            }
            else
            {
                p_rghEqnComp.ref() += hmm;
            }

            ++phasei;
        }

        solve
        (
            p_rghEqnComp
          + p_rghEqnIncomp,
            p_rgh.select(pimple.finalInnerIter())
        );

        if (pimple.finalNonOrthogonalIter())
        {
    const PtrList<volScalarField>& pcSources = mixture.phaseChangeSources();
    
    phasei = 0;
    for (phaseModel& phase : mixture.phases())
    {
        volScalarField compressionDgdt = 
            pos0(phase) * (p_rghEqnComps[phasei] & p_rgh) / phase.thermo().rho();
        
        // Add phase change if sources are available
        if (pcSources.size() > 0 && phasei < pcSources.size())
        {
            phase.dgdt() = compressionDgdt + pcSources[phasei];
        }
        else
        {
            phase.dgdt() = compressionDgdt;
        }
        
        ++phasei;
    }

            phi = phiHbyA + p_rghEqnIncomp.flux();

            U = HbyA
              + rAU()*fvc::reconstruct((phig + p_rghEqnIncomp.flux())/rAUf);
            U.correctBoundaryConditions();
        }



//         if (pimple.finalNonOrthogonalIter())
// {
//     // First calculate unconstrained dgdt for each phase
//     PtrList<volScalarField> dgdts(mixture.phases().size());
//     volScalarField sumDgdt
//     (
//         IOobject("sumDgdt", mesh.time().timeName(), mesh),
//         mesh,
//         dimensionedScalar(dimless/dimTime, Zero)
//     );
    
//     phasei = 0;
//     for (phaseModel& phase : mixture.phases())
//     {
//         dgdts.set
//         (
//             phasei,
//             new volScalarField
//             (
//                 "dgdt_" + phase.name(),
//                 pos0(phase) * (p_rghEqnComps[phasei] & p_rgh) / phase.thermo().rho()
//             )
//         );
        
//         // Calculate what the sum would be
//         sumDgdt += dgdts[phasei];
        
//         ++phasei;
//     }
    
//     // Now correct each dgdt to ensure Σ(dgdt) = 0 (maintains Σα = 1)
//     phasei = 0;
//     for (phaseModel& phase : mixture.phases())
//     {
//         // Subtract the mean to ensure sum = 0
//         // Weighted by volume fraction so larger phases absorb more correction
//         phase.dgdt() = dgdts[phasei] - sumDgdt;
        
//         ++phasei;
//     }
    
//     // Verify constraint
//     volScalarField checkSum
//     (
//         IOobject("checkSum", mesh.time().timeName(), mesh),
//         mesh,
//         dimensionedScalar(dimless/dimTime, Zero)
//     );
    
//     for (const phaseModel& phase : mixture.phases())
//     {
//         checkSum += phase.dgdt();
//     }
    
//     Info << "Sum of dgdt: mean = " << checkSum.weightedAverage(mesh.V()).value()
//          << ", max = " << max(mag(checkSum)).value() << endl;
    
//     // Rest of pEqn.H continues...
//     phi = phiHbyA + p_rghEqnIncomp.flux();
    
//     U = HbyA + rAU()*fvc::reconstruct((phig + p_rghEqnIncomp.flux())/rAUf);
//     U.correctBoundaryConditions();
// }




    }

    p = max(p_rgh + mixture.rho()*gh, pMin);

    // Update densities from change in p_rgh
    mixture.correctRho(p_rgh - p_rgh_0);
    rho = mixture.rho();

    // Correct p_rgh for consistency with p and the updated densities
    p_rgh = p - rho*gh;
    p_rgh.correctBoundaryConditions();

    K = 0.5*magSqr(U);

    Info<< "max(U) " << max(mag(U)).value() << endl;
    Info<< "min(p_rgh) " << min(p_rgh).value() << endl;
    Info<< "min(T) " << min(T).value() << endl;

    if (!correctPhi)
    {
        rAU.clear();
    }

}
