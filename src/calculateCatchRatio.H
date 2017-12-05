Info << nl << "Calculating specific catch ratio fields\n" << endl;

forAll(phases, i) {
  word name = ethz::phase_field_name("scr", i + 1);

  surfaceScalarField normalvel(
      mag((mesh.Sf() / mesh.magSf()) & fvc::interpolate(Urain[i])));
  surfaceScalarField surfaceScr((normalvel * fvc::interpolate(alpharain[i])) *
                                ((3600 * 1E3) / (Rh * phases[i][1])));

  const surfaceScalarField::Boundary& patchSurfaceScr =
      surfaceScr.boundaryField();

  volScalarField scrtemp(IOobject(name, runTime.timeName(), mesh), mesh,
                         dimensionedScalar(name, dimless, 0.0));

  volScalarField::Boundary& scrtempb = scrtemp.boundaryFieldRef();

  forAll(scrtempb, patchi) { scrtempb[patchi] = patchSurfaceScr[patchi]; }

  scr.set(i, scrtemp);
  scr[i].write();
}

Info << nl << "Calculating catch ratio field\n" << endl;

volScalarField gcr(  // initialize global catch ratio field
    IOobject("gcr", runTime.timeName(), mesh), mesh,
    dimensionedScalar("gcr", dimless, 0.0));

forAll(phases, i) { gcr += scr[i] * phases[i][1]; }
gcr.write();
