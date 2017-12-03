fvScalarMatrix alphaEqn
(
  fvm::ddt(alpharain[phase_no]) 
+ fvm::div(phirain[phase_no],alpharain[phase_no],"div(phirain,alpharain)")
+ fvm::Sp(fvc::div(phirain[phase_no]),alpharain[phase_no]) //added for bounded div scheme
);

alphaEqn.relax();

alphaEqn.solve();

alpharain[phase_no].correctBoundaryConditions();
