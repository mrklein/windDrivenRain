#ifndef UPDATEVALUES_H_
#define UPDATEVALUES_H_

#include "dimensionedScalar.H"
#include "scalar.H"
#include "volFieldsFwd.H"

namespace ethz {

extern void
getParameters(Foam::dimensionedScalar,
              Foam::dimensionedScalar&,
              Foam::dimensionedScalar&,
              Foam::dimensionedScalar&);

extern Foam::tmp<Foam::volScalarField>
getCdRe(const Foam::volScalarField&);

} // namespace ethz

#endif /* UPDATEVALUES_H_ */
