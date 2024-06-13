/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   MoveAndRotate.h
 * Author: savins
 *
 * Created on 9 de enero de 2016, 20:47
 */

#ifndef MOVEANDROTATE_H
#define MOVEANDROTATE_H
#include "model/Molecule.h"
#include "model/Point3DDouble.h"

class MoveAndRotate {
public:
  MoveAndRotate();
  void static MolToNewPosition(Molecule *mol, double deltaX, double deltaY,
                               double deltaZ);

  void static MolToNewPosition(double *atomsMolXYZ, unsigned int size,
                               double deltaX, double deltaY, double deltaZ);

  /*This is a implementation without any optimization*/
  void static PreciseMolToNewPosition(double *atomsMolXYZ, unsigned int size,
                                      double deltaX, double deltaY,
                                      double deltaZ);
  /**
   * Dada una molecula, dos puntos en el espacio y un angulo, se obtiene la
   * molecula rotada en un sentido.
   */
  void static RotateMolAccording1Axis(const Molecule *mol, double theta,
                                      Point3DDouble axisA, Point3DDouble axisB,
                                      Molecule *molRotated);

  void static RotateMolAccording1Axis(const double *molAtomsXYZ,
                                      unsigned int nAtoms, double theta,
                                      double axisAx, double axisAy,
                                      double axisAz, double axisBx,
                                      double axisBy, double axisBz,
                                      double *molAtomsRotated);

  /*This is a implementation without any optimization*/
  void static PreciseRotateMolAccording1Axis(const double *molAtomsXYZ,
                                             unsigned int nAtoms, double theta,
                                             double axisAx, double axisAy,
                                             double axisAz, double axisBx,
                                             double axisBy, double axisBz,
                                             double *molAtomsRotated);

  // void static RotateAtomAccordingAngleXYZ(Atom* atomRotate, double theta,
  // double x, double y, double z); void static RotateMolAccordingAngleXYZ(const
  // Molecule* mol, double theta, double x, double y,double z, Molecule*
  // molRotate);
  /**
   * Dada una molecula, se obtiene la molecula rotada para los dos �ngulos.
   */
  void static RotateMolAccording1Axis(const Molecule *mol, double theta,
                                      Point3DDouble axisA, Point3DDouble axisB,
                                      Molecule *molClockwise,
                                      Molecule *molAnticlockwise);
  void static RotateAtomAccording1Axis(Atom *atomClockwise,
                                       Atom *atomAnticlockwise, double theta,
                                       Point3DDouble axisA,
                                       Point3DDouble axisB);

  virtual ~MoveAndRotate();

private:
};

#endif /* MOVEANDROTATE_H */
