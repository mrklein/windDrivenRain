Info << nl << "Calculating specific catch ratio fields\n" << endl;

for (int phase_no = 0; phase_no < phases.size(); phase_no++) {
  word name = ethz::phase_field_name("scr", phase_no + 1);

  surfaceScalarField normalvel(
      mag((mesh.Sf() / mesh.magSf()) & fvc::interpolate(Urain[phase_no])));
  surfaceScalarField surfaceScr(
      (normalvel * fvc::interpolate(alpharain[phase_no])) *
      ((3600 * 1E3) / (Rh * phases[phase_no][1])));
  const surfaceScalarField::Boundary& patchSurfaceScr =
      surfaceScr.boundaryField();

  volScalarField scrtemp(IOobject(name, runTime.timeName(), mesh), mesh,
                         dimensionedScalar(name, dimless, 0.0));

  volScalarField::Boundary& scrtempb = scrtemp.boundaryFieldRef();

  forAll(scrtempb, patchi) { scrtempb[patchi] = patchSurfaceScr[patchi]; }

  scr.set(phase_no, scrtemp);
  scr[phase_no].write();
}

Info << nl << "Calculating catch ratio field\n" << endl;

volScalarField gcr(  // initialize global catch ratio field
    IOobject("gcr", runTime.timeName(), mesh), mesh,
    dimensionedScalar("gcr", dimless, 0.0));

for (int phase_no = 1; phase_no < phases.size(); phase_no++) {
  gcr += scr[phase_no] * phases[phase_no][1];
}
gcr.write();
