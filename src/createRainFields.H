PtrList<volVectorField> Urain(phases.size());
PtrList<surfaceScalarField> phirain(phases.size());
PtrList<volScalarField> alpharain(phases.size());
PtrList<volScalarField> scr(phases.size());

forAll(phases, i) {
  {
    word name = ethz::phase_field_name("U", i + 1);
    Urain.set(i, new volVectorField(
                     IOobject(name, runTime.timeName(), mesh,
                              IOobject::MUST_READ, IOobject::AUTO_WRITE),
                     mesh));
  }

  {
    word name = ethz::phase_field_name("phi", i + 1);
    phirain.set(i, new surfaceScalarField(
                       IOobject(name, runTime.timeName(), mesh),
                       linearInterpolate(Urain[i]) & mesh.Sf()));
  }

  {
    word name = ethz::phase_field_name("alpha", i + 1);
    alpharain.set(i, new volScalarField(
                         IOobject(name, runTime.timeName(), mesh,
                                  IOobject::MUST_READ, IOobject::AUTO_WRITE),
                         mesh));
  }
}
