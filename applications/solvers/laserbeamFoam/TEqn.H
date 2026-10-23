{
    int iter = 0;
    scalarField residual(epsilon1.primitiveField().size(), 1.0);
    scalar maxResidual = 1.0;
    scalar meanResidual = 1.0;
    //scalar Tresidual = 1.0;

    TRHS *= 0.0;

    rhoCp = rho*cp;
    rhoCp.correctBoundaryConditions();

    rhophicp = fvc::interpolate(cp)*rhoPhi;

    em = (epsilon1*emL) + ((1.0 - epsilon1)*emS);

    #include "checkTEqnAssembly.H"

    // The temperature-independent terms, assembled once per time step when
    // their coefficients do not depend on the temperature of the current
    // correction (see checkTEqnAssembly.H)
    tmp<fvScalarMatrix> tTEqnBase;
    if (TEqnAssembleOnce)
    {
        tTEqnBase =
            fvm::ddt(rhoCp,T)
          + fvm::div(rhophicp, T)
          - fvm::Sp(fvc::ddt(rhoCp) + fvc::div(rhophicp), T)
          - fvm::laplacian(kappa, T);
    }

    // const volTensorField gradU(fvc::grad(U));
    // volTensorField tau
    // (
    //     ((mixture.mu()*gradU) + (mixture.mu()*gradU.T()))
    //   - ((2.0/3.0)*mixture.mu()*(fvc::div(phi))*I)
    // );

    // The evaporation and radiation sources and the thermal damper enter the
    // equation multiplied by mag(gradAlpha), which is zero away from the
    // interface: they are evaluated, with the same expressions, only in the
    // cells where mag(gradAlpha) is not zero, and are zero elsewhere
    const scalarField magGradAlpha(mag(gradAlpha.primitiveField()));

    DynamicList<label> interfaceCells(magGradAlpha.size()/8);
    forAll(magGradAlpha, celli)
    {
        if (magGradAlpha[celli] > 0)
        {
            interfaceCells.append(celli);
        }
    }

    Qv = dimensionedScalar(Qv.dimensions(), Zero);

    volScalarField Qrad_expl
    (
        IOobject
        (
            "Qrad_expl",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(SB.dimensions()*pow4(dimTemperature), Zero)
    );

    volScalarField Qrad_inpl
    (
        IOobject
        (
            "Qrad_inpl",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(SB.dimensions()*pow3(dimTemperature), Zero)
    );

    if (damperSwitch)
    {
        thermalDamper = dimensionedScalar(thermalDamper.dimensions(), Zero);
    }

    // Constants of the expressions, combined in the same order as before
    const scalar QvCoeff = (0.82*LatentHeatVap*Mm*p0).value();
    const scalar QvExpCoeff = (LatentHeatVap*Mm).value();
    const scalar QvRootCoeff = (2*M_PI*Mm*R).value();
    const scalar TRef4 = Foam::pow(TRef, 4).value();
    const scalar TvapValue = Tvap.value();
    const scalar RValue = R.value();
    const scalar SBValue = SB.value();

    // apparentCapacity: in the cells of the mushy range the liquid fraction
    // is a function of T, epsilon1 = (T - TSolidus)/(TLiquidus - TSolidus), so
    // its latent heat is taken implicitly, rho*L/(TLiquidus - TSolidus)
    // *(T - T_k)/deltaT with T_k the temperature of the previous correction (a
    // Newton linearisation), and the correction loop converges in a few
    // corrections. The converged solution is the same as with sourceBased
    const bool newtonLatentHeat = (latentHeatMethod == "apparentCapacity");

    volScalarField latentCoeff
    (
        IOobject
        (
            "latentCoeff",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            rho.dimensions()*LatentHeat.dimensions()/dimTemperature/dimTime,
            Zero
        )
    );

    scalarField TPrev;

    do
    {
        epsilon1.storePrevIter();

        if (newtonLatentHeat)
        {
            // The cells in the mushy range, from the liquid fraction and the
            // temperature of the previous correction (a cell that leaves the
            // range is treated as solid or liquid in the next correction)
            scalarField& latentCoeffI = latentCoeff.primitiveFieldRef();
            const scalarField& TI = T;
            const scalarField& epsilon1I = epsilon1;
            const scalarField& TSolidusI = TSolidus;
            const scalarField& TLiquidusI = TLiquidus;
            const scalarField& rhoI = rho;
            const scalarField& LatentHeatI = LatentHeat;
            const scalar rDeltaT = 1.0/runTime.deltaTValue();

            forAll(latentCoeffI, celli)
            {
                const scalar range = TLiquidusI[celli] - TSolidusI[celli];
                const scalar eps = epsilon1I[celli];

                const bool mushy =
                    range > SMALL
                 && LatentHeatI[celli] > SMALL
                 && (
                        (eps > 0 && eps < 1)
                     || (eps <= 0 && TI[celli] > TSolidusI[celli])
                     || (eps >= 1 && TI[celli] < TLiquidusI[celli])
                    );

                latentCoeffI[celli] =
                    mushy ? rhoI[celli]*LatentHeatI[celli]/range*rDeltaT : 0;
            }

            TPrev = T.primitiveField();
        }

        TRHS =
            LatentHeat*((fvc::ddt(rho, epsilon1) + fvc::div(rhoPhi, epsilon1)));

        TRHS.correctBoundaryConditions();

        {
            scalarField& QvI = Qv.primitiveFieldRef();
            scalarField& Qrad_explI = Qrad_expl.primitiveFieldRef();
            scalarField& Qrad_inplI = Qrad_inpl.primitiveFieldRef();
            const scalarField& TI = T;
            const scalarField& emI = em;

            forAll(interfaceCells, i)
            {
                const label celli = interfaceCells[i];
                const scalar Tc = TI[celli];

                QvI[celli] =
                    QvCoeff
                   *Foam::exp(QvExpCoeff*((Tc - TvapValue)/((RValue*Tc)*TvapValue)))
                   /Foam::pow(QvRootCoeff*Tc, 0.5);

                Qrad_explI[celli] =
                    (SBValue*emI[celli])*(TRef4 + 3.0*Foam::pow(Tc, 4.0));

                Qrad_inplI[celli] =
                    -((4*SBValue)*emI[celli]*Foam::pow(Tc, 3.0));
            }
        }

        Qv.correctBoundaryConditions();

        if (damperSwitch)
        {
            scalarField& thermalDamperI = thermalDamper;
            const scalarField& rhoCpI = rhoCp;
            const scalarField& TI = T;
            forAll(interfaceCells, i)
            {
                const label celli = interfaceCells[i];
                thermalDamperI[celli] =
                    2.0*rhoCpI[celli]
                   /(
                       (rho1.value()*polycp_m.value(TI[celli]))
                     + (rho2.value()*polycp_g.value(TI[celli]))
                   );
            }

            thermalDamper.correctBoundaryConditions();
        }

        const volScalarField gradAlphaDamper (mag(gradAlpha)*thermalDamper);
    const volScalarField& Joulesource =
        mthd.valid()
      ? mthd->Joulesource()
      : zeroJoulesource;

        tmp<fvScalarMatrix> tTEqn;
        if (TEqnAssembleOnce)
        {
            tTEqn = tmp<fvScalarMatrix>::New(tTEqnBase());
        }
        else
        {
            tTEqn =
                fvm::ddt(rhoCp,T)
              + fvm::div(rhophicp, T)
              - fvm::Sp(fvc::ddt(rhoCp) + fvc::div(rhophicp), T)
              - fvm::laplacian(kappa, T);
        }
        fvScalarMatrix& TEqn = tTEqn.ref();

        TEqn -= fvm::Sp(Qrad_inpl*gradAlphaDamper, T);

        if (newtonLatentHeat)
        {
            // rho*L/(TLiquidus - TSolidus)*(T - T_k)/deltaT: implicit in T,
            // explicit in the temperature of the previous correction
            TEqn += fvm::Sp(latentCoeff, T);
            TEqn -= latentCoeff*T;
        }

        TEqn -=
        (
            laser.deposition()
          + Joulesource
          + (Qrad_expl - Qv)*gradAlphaDamper
          - TRHS
        );

        TEqn.relax();

        TEqn.solve();

        T.correctBoundaryConditions();

        // T.relax();

        if (newtonLatentHeat)
        {
            // Liquid fraction from the new temperature: the Newton update in
            // the cells of the mushy range, the equilibrium value elsewhere
            scalarField& epsilon1I = epsilon1.primitiveFieldRef();
            const scalarField& latentCoeffI = latentCoeff;
            const scalarField& TI = T;
            const scalarField& TSolidusI = TSolidus;
            const scalarField& TLiquidusI = TLiquidus;

            forAll(epsilon1I, celli)
            {
                const scalar range = TLiquidusI[celli] - TSolidusI[celli];

                if (range <= SMALL)
                {
                    continue;
                }

                const scalar eps =
                    latentCoeffI[celli] > 0
                  ? epsilon1I[celli] + (TI[celli] - TPrev[celli])/range
                  : (TI[celli] - TSolidusI[celli])/range;

                epsilon1I[celli] = max(min(eps, scalar(1)), scalar(0));
            }
        }
        else
        {
            Tcorr = (TLiquidus - TSolidus)*epsilon1 + TSolidus;
            Tcorr.correctBoundaryConditions();

            epsilon1 =
                max
                (
                    min
                    (
                        epsilon1 + (epsilonRel*cp/LatentHeat)*(T - Tcorr),
                        scalar(1)
                    ),
                    scalar(0)
                );
        }

        epsilon1.correctBoundaryConditions();

        residual =
            mag
            (
                epsilon1.primitiveField()
              - epsilon1.prevIter().primitiveField()
            );

        maxResidual = gMax(residual);

        meanResidual = gSum(residual*mesh.V())/gSum(mesh.V().field());

        //ddte1 = fvc::ddt(epsilon1);

        Info<< "Correcting epsilon1, mean residual = " << meanResidual
            << ", max residual = " << maxResidual
            << endl;
    }
    while
    (
        (
            ++iter < minTCorr
         // || (meanResidual > epsilonTol.value())
         || (maxResidual > epsilonTol.value())
        )
     && iter <= maxTCorr
    );

    gradT = fvc::grad(T);
    gradT.correctBoundaryConditions();  // ADDED - CRITICAL!
}
