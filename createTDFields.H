PtrList<volScalarField> Ctrain;
Ctrain.setSize(phases.size());

for (int phase_no = 0; phase_no < phases.size(); phase_no++) {
  string phase_Ctname = "Ct";
  string phase_no_str;
  OStringStream convert;
  convert << phase_no + 1;
  phase_no_str = convert.str();
  phase_Ctname.append(phase_no_str);

  volScalarField Cttemp(
    IOobject(phase_Ctname, runTime.timeName(), mesh, IOobject::READ_IF_PRESENT),
    mesh,
    dimensionedScalar(phase_Ctname, dimless, 0.0));
  Ctrain.set(phase_no, Cttemp);
}

volScalarField k(IOobject("k", runTime.timeName(), mesh, IOobject::MUST_READ),
                 mesh);

volScalarField epsilon(
  IOobject("epsilon", runTime.timeName(), mesh, IOobject::MUST_READ),
  mesh);

volScalarField nut(
  IOobject("nut", runTime.timeName(), mesh, IOobject::MUST_READ),
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

volScalarField nutrain(IOobject("nutrain", runTime.timeName(), mesh),
                       mesh,
                       dimensionedScalar("nutrain", dimViscosity, 1));

volScalarField tp(IOobject("tp", runTime.timeName(), mesh),
                  mesh,
                  dimensionedScalar("tp", dimTime, 0.0));

volScalarField tfl(IOobject("tfl", runTime.timeName(), mesh),
                   mesh,
                   dimensionedScalar("tfl", dimTime, 0));
