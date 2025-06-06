MRF.correctBoundaryVelocity(U);

if (laser.powderSim())
{
    DC =
        DarcyConstantlarge*Foam::pow((1.0-epsilon1mask),2)
       /(Foam::pow(epsilon1mask, 3) + DarcyConstantsmall);
}
else
{
    DC =
        DarcyConstantlarge*Foam::pow((1.0 - epsilon1), 2)
       /(Foam::pow(epsilon1, 3) + DarcyConstantsmall);
}

gradAlpha = fvc::grad(alpha1);

nHatM = gradAlpha/(mag(gradAlpha) + deltaN);

Marangoni =
    epsilon1*Marangoni_Constant*(gradT - nHatM*(nHatM & gradT))*mag(gradAlpha);

const scalar Tliquidus1Value = Tliquidus1.value();
const scalarField& TI = T;
vectorField& MarangoniI = Marangoni;
forAll(MarangoniI, celli)
{
    if(TI[celli] < Tliquidus1Value)
    {
        MarangoniI[celli] = vector::zero;
    }
}

Marangoni.correctBoundaryConditions();

if (damperSwitch)
{
    // To account for high density difference
    damper = 2.0*rho/(rho1 + rho2);
}

pVap = 0.54*p0*Foam::exp(LatentHeatVap*Mm*((T - Tvap)/(R*T*Tvap)));

fvVectorMatrix UEqn
(
    fvm::ddt(rho, U)
  + fvm::div(rhoPhi, U)
  + MRF.DDt(rho, U)
  + turbulence.divDevTau(rho, U)
  + fvm::Sp(DC, U)
  - Marangoni*damper
 ==
    phaseChange.SU(rho, rhoPhi, U)
  + fvModels.source(rho, U)
);

UEqn.relax();

fvConstraints.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve
    (
        UEqn
     ==
        fvc::reconstruct
        (
            (
                mixture.surfaceTensionForce()*fvc::interpolate(damper)
              + fvc::interpolate(pVap)
               *fvc::snGrad(alpha1)
               *fvc::interpolate(damper)
              - ghf*fvc::snGrad(rho*rhok)
              - fvc::snGrad(p_rgh)
            )*mesh.magSf()
        )
    );

    fvConstraints.constrain(U);
}
