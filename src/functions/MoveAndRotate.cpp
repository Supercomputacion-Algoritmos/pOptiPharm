/*
 * File:   MoveAndRotate.cpp
 * Author: savins
 *
 * Created on 9 de enero de 2016, 20:47
 */

#include "functions/MoveAndRotate.h"
#include "model/Point3DDouble.h"
#include "model/Quaternion.h"
#include <iostream>
#include <math.h>
#include <stdio.h>
MoveAndRotate::MoveAndRotate() {}
/**
 * Dada una molecula, se translada en espacio un delta X,Y,Z.
 */
void MoveAndRotate::MolToNewPosition(Molecule *mol, double deltaX,
                                     double deltaY, double deltaZ) {

  for (unsigned i = 0; i < mol->atoms.size(); i++) {

    mol->atoms[i].x += deltaX;
    mol->atoms[i].y += deltaY;
    mol->atoms[i].z += deltaZ;
  }
}

void MoveAndRotate::MolToNewPosition(double *atomsMolXYZ, unsigned int size,
                                     double deltaX, double deltaY,
                                     double deltaZ) {
  int iterador = 0;
  for (unsigned i = 0; i < size; i++) {
    iterador = i * 3;
    atomsMolXYZ[iterador] += deltaX;
    atomsMolXYZ[iterador + 1] += deltaY;
    atomsMolXYZ[iterador + 2] += deltaZ;
  }
}

void MoveAndRotate::PreciseMolToNewPosition(double *atomsMolXYZ,
                                            unsigned int size, double deltaX,
                                            double deltaY, double deltaZ) {
  int iterador = 0;
  for (unsigned i = 0; i < size; i++) {
    iterador = i * 3;
    atomsMolXYZ[iterador] += deltaX;
    atomsMolXYZ[iterador + 1] += deltaY;
    atomsMolXYZ[iterador + 2] += deltaZ;
  }
}

void MoveAndRotate::RotateMolAccording1Axis(const Molecule *mol, double theta,
                                            Point3DDouble axisA,
                                            Point3DDouble axisB,
                                            Molecule *molRotated) {

  // ACORDARME DE ESTO
  *molRotated = *mol;

  Point3DDouble vector = axisB - axisA;
  Quaternion q1;
  Quaternion q1conjugated;
  Point3DDouble atomPositionMinusA;

  double sinTheta = sin(theta / 2);

  double cosTheta = cos(theta / 2);

  double moduleVector =
      sqrt(vector.x * vector.x + vector.y * vector.y + vector.z * vector.z);

  q1.w = cosTheta;
  q1.x = sinTheta * vector.x / moduleVector;
  q1.y = sinTheta * vector.y / moduleVector;
  q1.z = sinTheta * vector.z / moduleVector;
  q1conjugated = Quaternion::conjugate(q1);

  Quaternion a = Quaternion::point3DToQuaternion(axisA);
  Quaternion pMinuxA;
  Quaternion p;
  Quaternion part2;
  Quaternion part3;

  for (unsigned i = 0; i < molRotated->atoms.size(); i++) {
    atomPositionMinusA.x = molRotated->atoms[i].x - axisA.x;
    atomPositionMinusA.y = molRotated->atoms[i].y - axisA.y;
    atomPositionMinusA.z = molRotated->atoms[i].z - axisA.z;
    pMinuxA = Quaternion::point3DToQuaternion(atomPositionMinusA);
    // CLOCKWISE
    part2 = q1 * pMinuxA;
    part3 = part2 * q1conjugated;
    p = a + part3;

    molRotated->atoms[i].x = p.x;
    molRotated->atoms[i].y = p.y;
    molRotated->atoms[i].z = p.z;
  }
}

inline float fastsqrt(float fval) {
  // Fast square root
  int ival;
  ival = (*(int *)&fval);
  ival &= 0x7FFFFFFF; // fixes sqrt(-0) bug
  ival -= 0x3f800000; // subtract 127
  ival >>= 1;         // requires signed shift to preserve sign (exponent/2)
  ival += 0x3f800000; // rebias new exponent
  fval = (*(float *)&ival);
  return fval;
}

inline double fastCos(double x, int len = 4) { // USE BETWEEN -pi, pi
  double sq_x = x * x, output = 1.0, qu_x = sq_x * sq_x, qu0 = qu_x;
  // double fac[] = {2.0, 24.0, 720.0, 40320.0, 3628800.0, 479001600.0};//2!,
  // 4!, 6!, 8!, 10!, 12!...
  double fac[] = {
      0.5,           0.041666667, 1.388888889e-3, 2.48015873e-5, 2.755731922e-7,
      2.087675699e-9}; // 2!, 4!, 6!, 8!, 10!, 12!...
  for (int i = 0; i < len; i = i + 2) {
    output -= (sq_x * fac[i]);
    sq_x = sq_x * qu_x;
  }
  for (int i = 1; i < len; i = i + 2) {
    output += (qu0 * fac[i]);
    qu0 = qu0 * qu_x;
  }
  return output;
}

inline double fastSin(double x, int len = 4) { // USE BETWEEN -pi, pi
  double tr_x = x * x * x, output = x, qu_x = tr_x * x;
  // double fac[] = { 6.0, 120.0, 5040.0, 362880.0, 39916800.0, 6227020800.0 };
  // //3!, 5!, 7!, 9!, 11!, 13!...
  double fac[] = {
      0.166666667,    8.333333333e-3, 1.984126984e-4, 2.755731922e-6,
      2.505210839e-8, 1.605904384e-10}; // 3!, 5!, 7!, 9!, 11!, 13!...

  for (int i = 0; i < len; i = i + 2) {
    output -= (tr_x * fac[i]);
    tr_x = tr_x * qu_x;
  }
  tr_x = qu_x * x;
  for (int i = 1; i < len; i = i + 2) {
    output += (tr_x * fac[i]);
    tr_x = tr_x * qu_x;
  }
  if (output < 0)
    return 1e-5;
  return output;
}

void MoveAndRotate::RotateMolAccording1Axis(const double *molAtomsXYZ,
                                            unsigned int nAtoms, double theta,
                                            double axisAx, double axisAy,
                                            double axisAz, double axisBx,
                                            double axisBy, double axisBz,
                                            double *molAtomsRotated) {

  double vectorX = axisBx - axisAx;
  double vectorY = axisBy - axisAy;
  double vectorZ = axisBz - axisAz;

  double q1w, q1x, q1y, q1z, q1wConjugated, q1xConjugated, q1yConjugated,
      q1zConjugated, atomPositionMinusAW, atomPositionMinusAX,
      atomPositionMinusAY, atomPositionMinusAZ;

  double angle = theta * 0.5;
  double sinTheta = sin(angle);

  // double cosTheta2 = cos(angle);

  double cosTheta = cos(angle);

  // double sinTheta2 = sin(angle);

  /*	double moduleVector = sqrt(
   vectorX * vectorX + vectorY * vectorY + vectorZ * vectorZ);
   */
  // FASTER IMPLEMENTATION
  double moduleVector =
      fastsqrt(vectorX * vectorX + vectorY * vectorY + vectorZ * vectorZ);

  q1w = cosTheta;
  q1x = sinTheta * vectorX / moduleVector;
  q1y = sinTheta * vectorY / moduleVector;
  q1z = sinTheta * vectorZ / moduleVector;

  q1wConjugated = q1w;
  q1xConjugated = -1 * q1x;
  q1yConjugated = -1 * q1y;
  q1zConjugated = -1 * q1z;

  double part2w, part2x, part2y, part2z;
  double part3w, part3x, part3y, part3z;
  int iterador = 0;

  for (unsigned i = 0; i < nAtoms; i++) {
    iterador = 3 * i;
    atomPositionMinusAW = 0;
    atomPositionMinusAX = molAtomsXYZ[iterador] - axisAx;
    atomPositionMinusAY = molAtomsXYZ[iterador + 1] - axisAy;
    atomPositionMinusAZ = molAtomsXYZ[iterador + 2] - axisAz;

    // CLOCKWISE
    // part2 = q1 * pMinuxA;
    part2w = q1w * atomPositionMinusAW - q1x * atomPositionMinusAX -
             q1y * atomPositionMinusAY - q1z * atomPositionMinusAZ;
    part2x = q1w * atomPositionMinusAX + q1x * atomPositionMinusAW +
             q1y * atomPositionMinusAZ - q1z * atomPositionMinusAY;
    part2y = q1w * atomPositionMinusAY - q1x * atomPositionMinusAZ +
             q1y * atomPositionMinusAW + q1z * atomPositionMinusAX;
    part2z = q1w * atomPositionMinusAZ + q1x * atomPositionMinusAY -
             q1y * atomPositionMinusAX + q1z * atomPositionMinusAW;

    // part3 = part2 * q1conjugated;
    part3w = part2w * q1wConjugated - part2x * q1xConjugated -
             part2y * q1yConjugated - part2z * q1zConjugated;
    part3x = part2w * q1xConjugated + part2x * q1wConjugated +
             part2y * q1zConjugated - part2z * q1yConjugated;
    part3y = part2w * q1yConjugated - part2x * q1zConjugated +
             part2y * q1wConjugated + part2z * q1xConjugated;
    part3z = part2w * q1zConjugated + part2x * q1yConjugated -
             part2y * q1xConjugated + part2z * q1wConjugated;

    // p = a + part3;

    molAtomsRotated[iterador] = axisAx + part3x;
    molAtomsRotated[iterador + 1] = axisAy + part3y;
    molAtomsRotated[iterador + 2] = axisAz + part3z;
  }
}

void MoveAndRotate::PreciseRotateMolAccording1Axis(
    const double *molAtomsXYZ, unsigned int nAtoms, double theta, double axisAx,
    double axisAy, double axisAz, double axisBx, double axisBy, double axisBz,
    double *molAtomsRotated) {

  double vectorX = axisBx - axisAx;
  double vectorY = axisBy - axisAy;
  double vectorZ = axisBz - axisAz;

  double q1w, q1x, q1y, q1z, q1wConjugated, q1xConjugated, q1yConjugated,
      q1zConjugated, atomPositionMinusAW, atomPositionMinusAX,
      atomPositionMinusAY, atomPositionMinusAZ;

  double angle = theta * 0.5;

  double sinTheta = sin(angle);
  double cosTheta = cos(angle);

  // FASTER IMPLEMENTATION
  double moduleVector =
      sqrt(vectorX * vectorX + vectorY * vectorY + vectorZ * vectorZ);

  q1w = cosTheta;
  q1x = sinTheta * vectorX / moduleVector;
  q1y = sinTheta * vectorY / moduleVector;
  q1z = sinTheta * vectorZ / moduleVector;

  q1wConjugated = q1w;
  q1xConjugated = -1 * q1x;
  q1yConjugated = -1 * q1y;
  q1zConjugated = -1 * q1z;

  double part2w, part2x, part2y, part2z;
  double part3w, part3x, part3y, part3z;
  int iterador = 0;

  for (unsigned i = 0; i < nAtoms; i++) {
    iterador = 3 * i;
    atomPositionMinusAW = 0;
    atomPositionMinusAX = molAtomsXYZ[iterador] - axisAx;
    atomPositionMinusAY = molAtomsXYZ[iterador + 1] - axisAy;
    atomPositionMinusAZ = molAtomsXYZ[iterador + 2] - axisAz;

    // CLOCKWISE
    // part2 = q1 * pMinuxA;
    part2w = q1w * atomPositionMinusAW - q1x * atomPositionMinusAX -
             q1y * atomPositionMinusAY - q1z * atomPositionMinusAZ;
    part2x = q1w * atomPositionMinusAX + q1x * atomPositionMinusAW +
             q1y * atomPositionMinusAZ - q1z * atomPositionMinusAY;
    part2y = q1w * atomPositionMinusAY - q1x * atomPositionMinusAZ +
             q1y * atomPositionMinusAW + q1z * atomPositionMinusAX;
    part2z = q1w * atomPositionMinusAZ + q1x * atomPositionMinusAY -
             q1y * atomPositionMinusAX + q1z * atomPositionMinusAW;

    // part3 = part2 * q1conjugated;
    part3w = part2w * q1wConjugated - part2x * q1xConjugated -
             part2y * q1yConjugated - part2z * q1zConjugated;
    part3x = part2w * q1xConjugated + part2x * q1wConjugated +
             part2y * q1zConjugated - part2z * q1yConjugated;
    part3y = part2w * q1yConjugated - part2x * q1zConjugated +
             part2y * q1wConjugated + part2z * q1xConjugated;
    part3z = part2w * q1zConjugated + part2x * q1yConjugated -
             part2y * q1xConjugated + part2z * q1wConjugated;

    // p = a + part3;
    molAtomsRotated[iterador] = axisAx + part3x;
    molAtomsRotated[iterador + 1] = axisAy + part3y;
    molAtomsRotated[iterador + 2] = axisAz + part3z;
  }
}

void MoveAndRotate::RotateMolAccording1Axis(const Molecule *mol, double theta,
                                            Point3DDouble axisA,
                                            Point3DDouble axisB,
                                            Molecule *molClockwise,
                                            Molecule *molAnticlockwise) {
  *molClockwise = *mol;
  *molAnticlockwise = *mol;

  for (unsigned i = 0; i < mol->atoms.size(); i++) {

    RotateAtomAccording1Axis(&molClockwise->atoms[i],
                             &molAnticlockwise->atoms[i], theta, axisA, axisB);
  }
}
/*
 void MoveAndRotate::RotateMolAccording1Axis(const Molecule* mol, double theta,
 Point3DDouble axisA, Point3DDouble axisB, Molecule* molRotated) {
 *molRotated = *mol;

 for (unsigned i = 0; i < mol->atoms.size(); i++) {

 RotateAtomAccording1Axis(&molRotated->atoms[i], theta, axisA, axisB);
 }

 }*/

/*void MoveAndRotate::RotateAtomAccordingAngleXYZ(Atom* atomRotate, double
 theta, double x, double y, double z) {

 Quaternion q1;
 Quaternion q1conjugated;
 Point3DDouble atomPositionMinusA;

 double sinTheta = sin(theta / 2);

 double cosTheta = cos(theta / 2);

 double moduleVector = sqrt(x * x + y *y + z * z);

 atomPositionMinusA.x = atomRotate->x - x;
 atomPositionMinusA.y = atomRotate->y - y;
 atomPositionMinusA.z = atomRotate->z - z;
 q1.w = cosTheta;
 q1.x = sinTheta * x / moduleVector;
 q1.y = sinTheta * y / moduleVector;
 q1.z = sinTheta * z / moduleVector;

 q1conjugated = Quaternion::conjugate(q1);

 Quaternion a = Quaternion::XYZToQuaternion(x,y,z);
 Quaternion pMinuxA = Quaternion::point3DToQuaternion(atomPositionMinusA);

 //CLOCKWISE
 Quaternion part2 = q1*pMinuxA;
 Quaternion part3 = part2*q1conjugated;
 Quaternion p = a + part3;

 atomRotate->x = p.x;
 atomRotate->y = p.y;
 atomRotate->z = p.z;
 }*/
/*void MoveAndRotate::RotateMolAccordingAngleXYZ(const Molecule* mol, double
 theta, double x, double y,double z, Molecule* molRotate) { *molRotate = *mol;

 for (int i = 0; i < mol->atoms.size(); i++) {

 RotateAtomAccordingAngleXYZ(&molRotate->atoms[i], theta, x,y,z);
 }

 }*/

/**
 * Dados dos puntos en el espacio y un angulo se devuelve un atomo rotado en una
 * direccion y otro en la otra
 */
void MoveAndRotate::RotateAtomAccording1Axis(Atom *atomClockwise,
                                             Atom *atomAnticlockwise,
                                             double theta, Point3DDouble axisA,
                                             Point3DDouble axisB) {
  Point3DDouble vector = axisB - axisA;
  Quaternion q1;
  Quaternion q1conjugated;
  Point3DDouble atomPositionMinusA;

  double sinTheta = sin(theta / 2);

  double cosTheta = cos(theta / 2);

  double moduleVector =
      sqrt(vector.x * vector.x + vector.y * vector.y + vector.z * vector.z);

  atomPositionMinusA.x = atomClockwise->x - axisA.x;
  atomPositionMinusA.y = atomClockwise->y - axisA.y;
  atomPositionMinusA.z = atomClockwise->z - axisA.z;
  q1.w = cosTheta;
  q1.x = sinTheta * vector.x / moduleVector;
  q1.y = sinTheta * vector.y / moduleVector;
  q1.z = sinTheta * vector.z / moduleVector;

  q1conjugated = Quaternion::conjugate(q1);

  Quaternion a = Quaternion::point3DToQuaternion(axisA);
  Quaternion pMinuxA = Quaternion::point3DToQuaternion(atomPositionMinusA);

  // CLOCKWISE
  Quaternion part2 = q1 * pMinuxA;
  Quaternion part3 = part2 * q1conjugated;
  Quaternion p = a + part3;

  atomClockwise->x = p.x;
  atomClockwise->y = p.y;
  atomClockwise->z = p.z;

  // ANTICLOCKWISE
  part2 = q1conjugated * pMinuxA;
  part3 = part2 * q1;
  p = a + part3;

  atomAnticlockwise->x = p.x;
  atomAnticlockwise->y = p.y;
  atomAnticlockwise->z = p.z;
}

MoveAndRotate::~MoveAndRotate() {}
