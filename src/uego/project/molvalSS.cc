#include "time.h"
#include "uego/uego.h"
//#include "internal.h"
#include "functions/Tanimoto.h"
#include <chrono>
#include <iomanip>
#include <stdlib.h>
#include <unistd.h>
////////////////////////////////////////////////////////////
// $Id: testval.cc,v 2.6 1998/03/29 10:40:32 jelasity Exp $
// testval.cc
// definition of NDimRealElement::Value()
// objective function library for n dimensional real spaces
// this version was originaly created for experimental tests
////////////////////////////////////////////////////////////
// modification history:
////////////////////////////////////////////////////////////

//-------------------------------------------------------------------------
#include "functions/MoveAndRotate.h"
#include "functions/Tanimoto.h"
#include "functions/VolumeOverlap.h"
#include <iostream>

// TODO Possible a memory leak
thread_local double *tmp_atoms =
    new double[INI.getMolVariable()->getNumAtoms() * 3];

inline double ValueInner(const double *x) {
  // Jero: My current fix consist in allocate a single memory block
  // per thread (lines 24-26) to handle the temp. atoms.

  //   double *newAtomsRotated =
  //       &(INI.getMolVariable()
  //             ->getAtomsXYZ()[INI.getMolVariable()->atoms.size() * 3]);

  double *newAtomsRotated = tmp_atoms;

  MoveAndRotate::PreciseRotateMolAccording1Axis(
      INI.getMolVariable()->getAtomsXYZ(), INI.getMolVariable()->atoms.size(),
      x[0], x[1], x[2], x[3], x[4], x[5], x[6], newAtomsRotated);

  MoveAndRotate::PreciseMolToNewPosition(
      newAtomsRotated, INI.getMolVariable()->atoms.size(), x[7], x[8], x[9]);

  /*
   *
   * ESTE METODO USA LOS PESOS DE LAS DOS MOLECULAS
   *
   *
   */
  double VAB = VolumeOverlap::preciseOverlapWEGA(
      INI.getMolQuery()->getAtomsXYZ(), INI.getMolQuery()->atoms.size(),
      INI.getMolQuery()->getWeightAtoms(), INI.getMolQuery()->getRadiusAtoms(),
      newAtomsRotated, INI.getMolVariable()->atoms.size(),
      INI.getMolVariable()->getWeightAtoms(),
      INI.getMolVariable()->getRadiusAtoms(), INI.getSameVanDerWaalsRadius());

  return Tanimoto::calculateTanimotoGeneric(
      INI.getMolQuery()->tanimoto, INI.getMolVariable()->tanimoto, VAB);
}

double NDimMolElement::Value() {
  if (x[1] == x[4] && x[2] == x[5] && x[3] == x[6]) {
    return 0;
  }

#ifdef OP_ENABLE_CUDA
  double result = 0.0;
  if (Acc::get_self().is_using_CUDA()) {
    // Big molecules. Use CUDA version
    result = Acc::get_self().ndim_value(this->x);
  } else {
    result = ValueInner(this->x);
  }

  return result;
#else
  return ValueInner(this->x);
#endif
};

double NDimMolElement::PreciseValue() {
  // time_t t1,t2;

  // t1=time(NULL);
  if (x[1] == x[4] && x[2] == x[5] && x[3] == x[6]) {
    return 0;
  }
  // NUEVO
  // double * newAtomsRotated =
  //		new double[INI.getMolVariable()->atoms.size() * 3];

  /*MoveAndRotate::RotateMolAccording1Axis(INI.getMolVariable()->getAtomsXYZ(),
   INI.getMolVariable()->atoms.size(), x[0], x[1], x[2], x[3], 0, 0,
   0, newAtomsRotated);*/
  double *newAtomsRotated =
      &(INI.getMolVariable()
            ->getAtomsXYZ()[INI.getMolVariable()->atoms.size() * 3]);
  MoveAndRotate::PreciseRotateMolAccording1Axis(
      INI.getMolVariable()->getAtomsXYZ(), INI.getMolVariable()->atoms.size(),
      x[0], x[1], x[2], x[3], x[4], x[5], x[6], newAtomsRotated);

  MoveAndRotate::PreciseMolToNewPosition(
      newAtomsRotated, INI.getMolVariable()->atoms.size(), x[7], x[8], x[9]);

  /*
   *
   * ESTE METODO USA LOS PESOS DE LAS DOS MOLECULAS
   *
   *
   */

  double VAB = VolumeOverlap::preciseOverlapWEGA(
      INI.getMolQuery()->getAtomsXYZ(), INI.getMolQuery()->atoms.size(),
      INI.getMolQuery()->getWeightAtoms(), INI.getMolQuery()->getRadiusAtoms(),
      newAtomsRotated, INI.getMolVariable()->atoms.size(),
      INI.getMolVariable()->getWeightAtoms(),
      INI.getMolVariable()->getRadiusAtoms(), INI.getSameVanDerWaalsRadius());

  double result = Tanimoto::calculateTanimotoGeneric(
      INI.getMolQuery()->tanimoto, INI.getMolVariable()->tanimoto, VAB);
  // cout << INI.getMolQuery()->tanimoto <<", "<<
  // INI.getMolVariable()->tanimoto<<", " << VAB<<", "<<result <<"\n"; delete[]
  // newAtomsRotated;

  return result;
};
