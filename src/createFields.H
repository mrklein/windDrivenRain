Info << "Reading transportProperties" << nl << endl;

IOdictionary transport(IOobject("transportProperties", runTime.constant(), mesh,
                                IOobject::MUST_READ_IF_MODIFIED));

dimensionedScalar Rh(transport.lookup("Rh"));
dimensionedScalar temp(transport.lookup("temp"));
dimensionedScalar rhoa("rhoa", dimDensity, 0.0);
dimensionedScalar rhop("rhop", dimDensity, 0.0);
dimensionedScalar mua("mua", dimMass / dimLength / dimTime, 0.0);

scalar scalingFactor = transport.lookupOrDefault("scalingFactor", 1.0);

const bool solveTD = transport.lookupOrDefault("solveTD", false);

List<List<scalar>> phases(transport.lookup("phases"));

Info << "Reading field U\n" << endl;
volVectorField U(IOobject("U", runTime.timeName(), mesh, IOobject::MUST_READ),
                 mesh);

// Perform scaling only at t=0
if (runTime.timeName() == "0") {
  U = U * scalingFactor;
  U.correctBoundaryConditions();
}

surfaceScalarField phiwind(IOobject("phiwind", runTime.timeName(), mesh,
                                    IOobject::READ_IF_PRESENT),
                           linearInterpolate(U) & mesh.Sf());

volScalarField Re(IOobject("Re", runTime.timeName(), mesh,
                           IOobject::READ_IF_PRESENT),
                  mesh, dimensionedScalar("Re", dimless, 500));

volScalarField CdRe(IOobject("CdRe", runTime.timeName(), mesh,
                             IOobject::READ_IF_PRESENT),
                    mesh, dimensionedScalar("CdRe", dimless, 250));
