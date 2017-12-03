Info << nl << "Calculating specific catch ratio fields\n" << endl;

for (int phase_no = 0; phase_no < phases.size(); phase_no++) {
  string phase_scrname = "scr";
  string phase_no_str;
  OStringStream convert;
  convert << phase_no + 1;
  phase_no_str = convert.str();
  phase_scrname.append(phase_no_str);

  surfaceScalarField normalvel(
      mag((mesh.Sf() / mesh.magSf()) & fvc::interpolate(Urain[phase_no])));
  surfaceScalarField surfaceScr(
      (normalvel * fvc::interpolate(alpharain[phase_no])) *
      ((3600 * 1E3) / (Rh * phases[phase_no][1])));
  const surfaceScalarField::Boundary& patchSurfaceScr =
      surfaceScr.boundaryField();

  volScalarField scrtemp(IOobject(phase_scrname, runTime.timeName(), mesh),
                         mesh, dimensionedScalar(phase_scrname, dimless, 0.0));

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
