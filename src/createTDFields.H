PtrList<volScalarField> Ctrain;
Ctrain.setSize(phases.size());

forAll(phases, i) {
  const word name = ethz::phase_field_name("Ct", i + 1);
  Ctrain.set(i, new volScalarField(
                    IOobject(name, runTime.timeName(), mesh,
                             IOobject::READ_IF_PRESENT),
                    mesh, dimensionedScalar(name, dimless, 0.0)));
}

volScalarField k(IOobject("k", runTime.timeName(), mesh, IOobject::MUST_READ),
                 mesh);

volScalarField epsilon(IOobject("epsilon", runTime.timeName(), mesh,
                                IOobject::MUST_READ),
                       mesh);

volScalarField nut(IOobject("nut", runTime.timeName(), mesh,
                            IOobject::MUST_READ),
                   mesh);

// Perform scaling only at t = 0
if (runTime.timeName() == "0") {
  k = k * pow(scalingFactor, 2);
  k.correctBoundaryConditions();
  epsilon = epsilon * pow(scalingFactor, 3);
  epsilon.correctBoundaryConditions();
  nut = nut * scalingFactor;
  nut.correctBoundaryConditions();
}

volScalarField nutrain(IOobject("nutrain", runTime.timeName(), mesh), mesh,
                       dimensionedScalar("nutrain", dimViscosity, 1));

volScalarField tp(IOobject("tp", runTime.timeName(), mesh), mesh,
                  dimensionedScalar("tp", dimTime, 0.0));

volScalarField tfl(IOobject("tfl", runTime.timeName(), mesh), mesh,
                   dimensionedScalar("tfl", dimTime, 0));
