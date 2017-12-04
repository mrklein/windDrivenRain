#include "updateValues.H"

#include "volFields.H"

// Linear interpolation
inline Foam::scalar interpolate(Foam::label i, Foam::scalar v,
                                const Foam::scalarList& f,
                                const Foam::scalarList& x) {
  const auto df = f[i + 1] - f[i];
  const auto dx = x[i + 1] - x[i];
  const auto dv = v - x[i];
  return f[i] + df / dx * dv;
}

// Check if we are in between i and i + 1 in x
inline bool in(const Foam::label i, const Foam::scalar v,
               const Foam::scalarList& x) {
  return (x[i] <= v && v < x[i + 1]);
}

// Interpolate rhoa, mua, and rhop from temperature
void ethz::getParameters(Foam::dimensionedScalar temp_,
                         Foam::dimensionedScalar& rhoa_,
                         Foam::dimensionedScalar& mua_,
                         Foam::dimensionedScalar& rhop_) {
  using namespace Foam;

  const scalar TC = temp_.value() - 273.15;

  scalarList T_1_M({-18.0, -6.7, 0.0, 4.4, 15.6, 26.7, 37.8});
  scalarList RHO_A_M({1.38, 1.32, 1.293, 1.27, 1.22, 1.18, 1.13});
  scalarList MU_A_M({0.0157E-3, 0.0168E-3, 0.0171E-3, 0.0173E-3, 0.0179E-3,
                     0.0184E-3, 0.0190E-3});
  scalarList T_2_M({0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0,
                    22.0, 24.0, 26.0, 28.0, 30.0, 35.0, 40.0});
  scalarList RHO_P_M({999.87, 999.97, 1000.00, 999.97, 999.88, 999.73, 999.52,
                      999.27, 998.97, 998.62, 998.23, 997.80, 997.33, 996.81,
                      996.26, 995.68, 994.00, 992.00});

  if (TC < Foam::min(T_1_M)) {
    Info << "Air temperature is too low!";
  } else if (TC >= max(T_1_M)) {
    Info << "Air temperature is too high!";
  } else {
    forAll(T_1_M, i) {
      if (not in(i, TC, T_1_M)) continue;
      mua_.value() = interpolate(i, TC, MU_A_M, T_1_M);
      rhoa_.value() = interpolate(i, TC, RHO_A_M, T_1_M);
      break;
    }
  }

  if (TC <= min(T_2_M)) {
    Info << "Air temperature is too low!";
  } else if (TC >= max(T_2_M)) {
    Info << "Air temperature is too high!";
  } else {
    forAll(T_2_M, i) {
      if (not in(i, TC, T_2_M)) continue;
      rhop_.value() = interpolate(i, TC, RHO_P_M, T_2_M);
      break;
    }
  }
}

// Interpolate CdRe volume field
Foam::tmp<Foam::volScalarField> ethz::getCdRe(
    const Foam::volScalarField& Reynolds) {
  using namespace Foam;

  volScalarField* CdRe(new volScalarField(Reynolds));

  //-------------------------------------------------------------------------
  // Matrices containing Reynolds numbers and the corresponding drag
  // coefficient values according to Gunn & Kinzer
  scalarList Re_M({1.80,   9.61,   23.4,   43.2,   68.7,   98.9,   134.0,
                   175.0,  220.0,  269.0,  372.0,  483.0,  603.0,  731.0,
                   866.0,  1013.0, 1164.0, 1313.0, 1461.0, 1613.0, 1764.0,
                   1915.0, 2066.0, 2211.0, 2357.0, 2500.0, 2636.0, 2772.0,
                   2905.0, 3033.0, 3164.0, 3293.0, 3423.0, 3549.0});
  scalarList Cd_M({15.0,  4.2,   2.4,   1.66,  1.28,  1.07,  0.926,
                   0.815, 0.729, 0.671, 0.607, 0.570, 0.545, 0.528,
                   0.517, 0.504, 0.495, 0.494, 0.498, 0.503, 0.511,
                   0.520, 0.529, 0.544, 0.559, 0.575, 0.594, 0.615,
                   0.635, 0.660, 0.681, 0.700, 0.727, 0.751});
  //-------------------------------------------------------------------------

  forAll(Reynolds, celli) {
    const scalar Rei = Reynolds[celli];
    if (Rei < min(Re_M)) {
      (*CdRe)[celli] = interpolate(0, Rei, Cd_M, Re_M)*Rei;
    } else if (Rei >= max(Re_M)) {
      (*CdRe)[celli] = interpolate(32, Rei, Cd_M, Re_M)*Rei;
    } else {
      forAll(Re_M, i) {
        if (not in(i, Rei, Re_M)) continue;
        (*CdRe)[celli] = interpolate(i, Rei, Cd_M, Re_M)*Rei;
        break;
      }
    }
  }

  return tmp<volScalarField>(CdRe);
}
