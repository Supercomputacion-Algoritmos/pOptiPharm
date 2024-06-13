/*
 * File:   VolumeOverlap.cpp
 * Author: savins
 *
 * Created on 2 de enero de 2016, 10:45
 */

#include "functions/VolumeOverlap.h"
#include <iostream>
#include <math.h>
#include <stdio.h>

double distsqr(const double x, const double y, const double z,
               const Atom &atom) {
  double tmpx = atom.x - x;
  double tmpy = atom.y - y;
  double tmpz = atom.z - z;
  return tmpx * tmpx + tmpy * tmpy + tmpz * tmpz;
}

double distsqr(Atom atom1, Atom atom2) {
  double tmpx = atom1.x - atom2.x;
  double tmpy = atom1.y - atom2.y;
  double tmpz = atom1.z - atom2.z;
  return tmpx * tmpx + tmpy * tmpy + tmpz * tmpz;
}

double VolumeOverlap::calculateVolumeOverlapHardSphere(Molecule *fixedMol,
                                                       Molecule *variableMol) {
  // Calculate bounding box
  double xmin, xmax, ymin, ymax, zmin, zmax;
  xmin = fixedMol->atoms[0].x - fixedMol->atoms[0].radiusVanDerWaals;
  ymin = fixedMol->atoms[0].y - fixedMol->atoms[0].radiusVanDerWaals;
  zmin = fixedMol->atoms[0].z - fixedMol->atoms[0].radiusVanDerWaals;

  xmax = fixedMol->atoms[0].x + fixedMol->atoms[0].radiusVanDerWaals;
  ymax = fixedMol->atoms[0].y + fixedMol->atoms[0].radiusVanDerWaals;
  zmax = fixedMol->atoms[0].z + fixedMol->atoms[0].radiusVanDerWaals;

  for (unsigned i = 0; i < fixedMol->atoms.size(); i++) {
    xmin =
        min(xmin, fixedMol->atoms[i].x - fixedMol->atoms[i].radiusVanDerWaals);
    ymin =
        min(ymin, fixedMol->atoms[i].y - fixedMol->atoms[i].radiusVanDerWaals);
    zmin =
        min(zmin, fixedMol->atoms[i].z - fixedMol->atoms[i].radiusVanDerWaals);

    xmax =
        max(xmax, fixedMol->atoms[i].x + fixedMol->atoms[i].radiusVanDerWaals);
    ymax =
        max(ymax, fixedMol->atoms[i].y + fixedMol->atoms[i].radiusVanDerWaals);
    zmax =
        max(zmax, fixedMol->atoms[i].z + fixedMol->atoms[i].radiusVanDerWaals);
  }

  for (unsigned i = 0; i < variableMol->atoms.size(); i++) {
    xmin = min(xmin, variableMol->atoms[i].x -
                         variableMol->atoms[i].radiusVanDerWaals);
    ymin = min(ymin, variableMol->atoms[i].y -
                         variableMol->atoms[i].radiusVanDerWaals);
    zmin = min(zmin, variableMol->atoms[i].z -
                         variableMol->atoms[i].radiusVanDerWaals);

    xmax = max(xmax, variableMol->atoms[i].x +
                         variableMol->atoms[i].radiusVanDerWaals);
    ymax = max(ymax, variableMol->atoms[i].y +
                         variableMol->atoms[i].radiusVanDerWaals);
    zmax = max(zmax, variableMol->atoms[i].z +
                         variableMol->atoms[i].radiusVanDerWaals);
  }

  // const double res = 0.33;
  const double res = 0.01;
  int count = 0;

  for (double xcor = xmin; xcor < xmax; xcor += res) {
    for (double ycor = ymin; ycor < ymax; ycor += res) {
      for (double zcor = zmin; zcor < zmax; zcor += res) {
        bool refAtomFound = false;
        for (unsigned i = 0; i < fixedMol->atoms.size() && !refAtomFound; i++) {
          refAtomFound = (fixedMol->atoms[i].radiusVanDerWaals *
                              fixedMol->atoms[i].radiusVanDerWaals >=
                          distsqr(xcor, ycor, zcor, fixedMol->atoms[i]));
        }
        if (!refAtomFound)
          continue;
        bool fitAtomFound = false;
        for (unsigned i = 0; i < variableMol->atoms.size() && !fitAtomFound;
             i++) {
          fitAtomFound = (variableMol->atoms[i].radiusVanDerWaals *
                              variableMol->atoms[i].radiusVanDerWaals >=
                          distsqr(xcor, ycor, zcor, variableMol->atoms[i]));
        }
        if (fitAtomFound)
          count++;
      }
    }
  }
  return count * res * res * res;
}

double VolumeOverlap::calculateOverlapVolumeGaussAnalytic(Molecule *fixed,
                                                          Molecule *variable) {
  const double partialalpha = 2.41798793102;
  const double pi = 3.14159265358;
  double overlap = 0.0;

  for (unsigned i = 0; i < fixed->atoms.size(); i++) {
    Atom atomFixed = fixed->atoms[i];
    double alphai = partialalpha /
                    (atomFixed.radiusVanDerWaals * atomFixed.radiusVanDerWaals);
    for (unsigned j = 0; j < variable->atoms.size(); j++) {
      Atom atomVariable = variable->atoms[j];
      double Rij2 = distsqr(atomFixed, atomVariable);

      // if (Rij2 > 16) continue;

      double alphaj = partialalpha / (atomVariable.radiusVanDerWaals *
                                      atomVariable.radiusVanDerWaals);
      double Kij = exp(-(alphai * alphaj * Rij2) / (alphai + alphaj));
      double Vij = 8 * Kij * pow((pi / (alphai + alphaj)), 1.5);
      // p = 2.7 instead of 2sqrt(2)
      //  double Vij = 7.29 * Kij * pow((pi / (alphai + alphaj)), 1.5);

      // printf("Vij from ref %d (w=%f), fit %d (w=%f) =
      // %f\n",i,atomFixed.radiusVanDerWaals,j,atomVariable.radiusVanDerWaals,Vij);
      overlap += Vij;
    }
  }
  return overlap;
}

double
VolumeOverlap::calculateOverlapVolumeGaussAnalyticWEGA(Molecule *fixed,
                                                       Molecule *variable) {

  const double pi = 3.14159265358;
  double overlap = 0.0;

  for (unsigned i = 0; i < fixed->atoms.size(); i++) {
    Atom atomFixed = fixed->atoms[i];
    double alphai = pow(
        3 * 2 * sqrt(2) * sqrt(pi) /
            (4 * (atomFixed.radiusVanDerWaals * atomFixed.radiusVanDerWaals *
                  atomFixed.radiusVanDerWaals)),
        0.67);
    for (unsigned j = 0; j < variable->atoms.size(); j++) {
      Atom atomVariable = variable->atoms[j];
      double Rij2 = distsqr(atomFixed, atomVariable);

      // if (Rij2 > 16) continue;

      double alphaj = pow(3 * 2 * sqrt(2) * sqrt(pi) /
                              (4 * (atomVariable.radiusVanDerWaals *
                                    atomVariable.radiusVanDerWaals *
                                    atomVariable.radiusVanDerWaals)),
                          0.67);
      double Kij = exp(-(alphai * alphaj * Rij2) / (alphai + alphaj));
      double Vij = 8 * Kij * pow((pi / (alphai + alphaj)), 1.5);
      // p = 2.7 instead of 2sqrt(2)
      //  double Vij = 7.29 * Kij * pow((pi / (alphai + alphaj)), 1.5);

      // printf("Vij from ref %d (w=%f), fit %d (w=%f) =
      // %f\n",i,atomFixed.radiusVanDerWaals,j,atomVariable.radiusVanDerWaals,Vij);
      overlap += Vij;
    }
  }
  return overlap;
}

double
VolumeOverlap::calculateOverlapVolumeGaussAnalyticWeight(Molecule *fixed,
                                                         Molecule *variable) {
  const double partialalpha = 2.41798793102;
  const double pi = 3.14159265358;
  double overlap = 0.0;

  for (unsigned i = 0; i < fixed->atoms.size(); i++) {
    Atom atomFixed = fixed->atoms[i];
    double alphai = partialalpha /
                    (atomFixed.radiusVanDerWaals * atomFixed.radiusVanDerWaals);
    for (unsigned j = 0; j < variable->atoms.size(); j++) {
      Atom atomVariable = variable->atoms[j];
      double Rij2 = distsqr(atomFixed, atomVariable);

      // if (Rij2 > 16) continue;

      double alphaj = partialalpha / (atomVariable.radiusVanDerWaals *
                                      atomVariable.radiusVanDerWaals);
      double Kij = exp(-(alphai * alphaj * Rij2) / (alphai + alphaj));
      double Vij = 8 * Kij * pow((pi / (alphai + alphaj)), 1.5);
      // p = 2.7 instead of 2sqrt(2)
      //  double Vij = 7.29 * Kij * pow((pi / (alphai + alphaj)), 1.5);

      // Vij with the weight
      double WeightVij = calculateWeight(&atomFixed, &fixed->atoms, i) *
                         calculateWeight(&atomVariable, &variable->atoms, j) *
                         Vij;
      // printf("Vij from ref %d (w=%f), fit %d (w=%f) =
      // %f\n",i,atomFixed.radiusVanDerWaals,j,atomVariable.radiusVanDerWaals,Vij);
      overlap += WeightVij;
    }
  }
  return overlap;
}
/** ORIGINAL
 double VolumeOverlap::calculateOverlapVolumeGaussAnalyticWeightWEGA(
 Molecule* fixed, Molecule* variable) {

 const double pi = 3.14159265358;
 double overlap = 0.0;
 double p = 2 * sqrt(2);

 for (unsigned i = 0; i < fixed->atoms.size(); i++) {
 Atom atomFixed = fixed->atoms[i];
 double alphai = pow(
 (3 * p * sqrt(pi)
 / (4
 * (atomFixed.radiusVanDerWaals
 * atomFixed.radiusVanDerWaals
 * atomFixed.radiusVanDerWaals))), 0.67);
 for (unsigned j = 0; j < variable->atoms.size(); j++) {
 Atom atomVariable = variable->atoms[j];
 double Rij2 = distsqr(atomFixed, atomVariable);

 //if (Rij2 > 16) continue;

 double alphaj = pow(
 3 * p * sqrt(pi)
 / (4
 * (atomVariable.radiusVanDerWaals
 * atomVariable.radiusVanDerWaals
 * atomVariable.radiusVanDerWaals)),
 0.67);

 double Kij = exp(-(alphai * alphaj * Rij2) / (alphai + alphaj));

 double Vij = p * p * Kij * pow((pi / (alphai + alphaj)), 1.5);

 //p = 2.7 instead of 2sqrt(2)
 // double Vij = 7.29 * Kij * pow((pi / (alphai + alphaj)), 1.5);
 double WeightVij = calculateWeight(&atomFixed, &fixed->atoms, i)
 * calculateWeight(&atomVariable, &variable->atoms, j) * Vij;
 //printf("Vij from ref %d (w=%f), fit %d (w=%f) =
 %f\n",i,atomFixed.radiusVanDerWaals,j,atomVariable.radiusVanDerWaals,Vij);

 overlap += WeightVij;

 }
 }
 return overlap;
 }*/

double VolumeOverlap::calculateOverlapVolumeGaussAnalyticWeightWEGA(
    Molecule *fixed, Molecule *variable) {

  double overlap = 0.0;

  for (unsigned i = 0; i < fixed->atoms.size(); i++) {
    for (unsigned j = 0; j < variable->atoms.size(); j++) {
      // Calculo Vij con el metodo porque he visto que se calculaba de la misma
      // forma que en el caso del peso.
      double Vij = calculateOverlapVolumeAtomsGaussAnalytic(
          &(fixed->atoms[i]), &(variable->atoms[j]));

      double WeightVij =
          calculateWeight(&(fixed->atoms[i]), &fixed->atoms, i) *
          calculateWeight(&(variable->atoms[j]), &variable->atoms, j) * Vij;
      overlap += WeightVij;
      // overlap += Vij;
    }
  }
  return overlap;
}

double VolumeOverlap::overlapROCS(Molecule *fixed, Molecule *variable) {

  const double partialalpha = 2.41798793102;
  const double pi = 3.14159265358;
  double overlap = 0.0;

  for (unsigned i = 0; i < fixed->atoms.size(); i++) {
    Atom atomFixed = fixed->atoms[i];
    double alphai = partialalpha /
                    (atomFixed.radiusVanDerWaals * atomFixed.radiusVanDerWaals);
    for (unsigned j = 0; j < variable->atoms.size(); j++) {
      Atom atomVariable = variable->atoms[j];
      double Rij2 = distsqr(atomFixed, atomVariable);

      // if (Rij2 > 16) continue;

      double alphaj = partialalpha / (atomVariable.radiusVanDerWaals *
                                      atomVariable.radiusVanDerWaals);
      double Kij = exp(-(alphai * alphaj * Rij2) / (alphai + alphaj));
      double Vij = 8 * Kij * pow((pi / (alphai + alphaj)), 1.5);
      overlap += Vij;
    }
  }
  return overlap;
}

double VolumeOverlap::overlapWEGA(Molecule *fixed, Molecule *variable,
                                  bool sameVanDerWaalsRadius) {
  double overlap = 0.0;

  for (unsigned i = 0; i < fixed->atoms.size(); i++) {
    for (unsigned j = 0; j < variable->atoms.size(); j++) {
      // Calculo Vij con el metodo porque he visto que se calculaba de la misma
      // forma que en el caso del peso.
      double Vij = calculateOverlapVolumeAtomsWEGA(
          &(fixed->atoms[i]), &(variable->atoms[j]), sameVanDerWaalsRadius);

      double WeightVij =
          calculateWeightWEGA(&(fixed->atoms[i]), &fixed->atoms, i,
                              sameVanDerWaalsRadius) *
          calculateWeightWEGA(&(variable->atoms[j]), &variable->atoms, j,
                              sameVanDerWaalsRadius) *
          Vij;

      overlap += WeightVij;
    }
  }

  return overlap;
}

#include <unistd.h>
double VolumeOverlap::overlapWEGA(double *fixed, unsigned int fixedSize,
                                  double *weightFixed, double *radiusFixed,
                                  double *variable, unsigned int variableSize,
                                  double *weightVariable,
                                  double *radiusVariable,
                                  bool sameVanDerWaalsRadius) {

  double overlap = 0.0;

  //#pragma omp parallel for num_threads(2) reduction(+:overlap)
  for (unsigned int i = 0; i < fixedSize; i++) {
    int iteradori = i * 3;
    for (unsigned int j = 0; j < variableSize; j++) {
      int iteradorj = j * 3;
      // Calculo Vij con el metodo porque he visto que se calculaba de la misma
      // forma que en el caso del peso.
      double Vij = calculateOverlapVolumeAtomsWEGA(
          fixed[iteradori], fixed[iteradori + 1], fixed[iteradori + 2],
          variable[iteradorj], variable[iteradorj + 1], variable[iteradorj + 2],
          sameVanDerWaalsRadius, radiusFixed[i], radiusVariable[j]);
      double WeightVij = weightFixed[i] * weightVariable[j] * Vij;
      // double WeightVij = 1 * Vij;
      overlap += WeightVij;
      /*	if (overlap > 1000){
       cout << "############################\n";
       cout << "Vij: "<<Vij << ", " << Vij2<<endl;
       cout << "Overlap: "<<overlap<<endl;
       cout << fixed[iteradori]<<", "<<fixed[iteradori + 1]<<", "<<
       fixed[iteradori + 2]<<", "<<variable[iteradorj]<<", "<<
       variable[iteradorj + 1]
       <<", "<<variable[iteradorj + 2]<<endl;
       cout <<"----------------------------\n";

       usleep(100000);

       }*/
    }
  }

  // cout << overlap << ", " << overlap2<<endl;

  return overlap;
}

double VolumeOverlap::preciseOverlapWEGA(double *fixed, unsigned int fixedSize,
                                         double *weightFixed,
                                         double *radiusFixed, double *variable,
                                         unsigned int variableSize,
                                         double *weightVariable,
                                         double *radiusVariable,
                                         bool sameVanDerWaalsRadius) {

  double overlap = 0.0;

  for (unsigned int i = 0; i < fixedSize; i++) {
    int iteradori = i * 3;
    for (unsigned int j = 0; j < variableSize; j++) {
      int iteradorj = j * 3;
      // Calculo Vij con el metodo porque he visto que se calculaba de la misma
      // forma que en el caso del peso.
      double Vij = preciseCalculateOverlapVolumeAtomsWEGA(
          fixed[iteradori], fixed[iteradori + 1], fixed[iteradori + 2],
          variable[iteradorj], variable[iteradorj + 1], variable[iteradorj + 2],
          sameVanDerWaalsRadius, radiusFixed[i], radiusVariable[j]);

      overlap += weightFixed[i] * weightVariable[j] * Vij;
    }
  }

  return overlap;
}

/*NON ANALYTICS OVVERLAP VOLUME GAUSSIAN*/
inline float distsqr(Point3DDouble &a, Atom &b) {
  float tmpx = b.x - a.x;
  float tmpy = b.y - a.y;
  float tmpz = b.z - a.z;
  return tmpx * tmpx + tmpy * tmpy + tmpz * tmpz;
}

double rho(Atom atom, Point3DDouble gc) {
  double rt22 = 2.82842712475f;
  double partialalpha = -2.41798793102f;
  double alpha =
      partialalpha / (atom.radiusVanDerWaals * atom.radiusVanDerWaals);
  float r2 = distsqr(gc, atom);
  return rt22 * exp(alpha * r2);
}

Point3DDouble gridSpaceToRealSpace(double res, Point3DDouble &lb, int x, int y,
                                   int z) {
  Point3DDouble rs;
  rs.x = lb.x + (x + 0.5) * res;
  rs.y = lb.y + (y + 0.5) * res;
  rs.z = lb.z + (z + 0.5) * res;
  return rs;
}

double VolumeOverlap::getOverlapVolumeGaussianNonAnalytic(Molecule &ref,
                                                          Molecule &fit,
                                                          Box &grid) {
  double volume = 0.0;

  for (int x = 0; x < grid.extent.x; x++) {
    for (int y = 0; y < grid.extent.y; y++) {
      for (int z = 0; z < grid.extent.z; z++) {
        double refgridval = 1.0;
        double fitgridval = 1.0;
        Point3DDouble gc =
            gridSpaceToRealSpace(grid.res, grid.lowestPoint, x, y, z);
        for (unsigned i = 0; i < ref.atoms.size(); i++) {
          refgridval *= (1 - rho(ref.atoms[i], gc));
        }
        refgridval = 1 - refgridval;
        for (unsigned i = 0; i < fit.atoms.size(); i++) {
          fitgridval *= (1 - rho(fit.atoms[i], gc));
        }
        fitgridval = 1 - fitgridval;
        volume += refgridval * fitgridval;
      }
    }
  }
  return volume * grid.res * grid.res * grid.res;
}

/*HERE THERE ARE THE IMPLEMENTATIONS ABOUT GRID*/
void VolumeOverlap::boundingBox(Molecule &mol, double margin, Point3DDouble &lb,
                                Point3DDouble &ub) {
  Atom &atom1 = mol.atoms[0];
  float rad = atom1.radiusVanDerWaals;
  lb.x = atom1.x - rad;
  lb.y = atom1.y - rad;
  lb.z = atom1.z - rad;
  ub.x = atom1.x + rad;
  ub.y = atom1.y + rad;
  ub.z = atom1.z + rad;
  for (unsigned int i = 1; i < mol.atoms.size(); i++) {
    Atom &atomi = mol.atoms[i];
    rad = atomi.radiusVanDerWaals;
    lb.x = min(lb.x, atomi.x - rad);
    lb.y = min(lb.y, atomi.y - rad);
    lb.z = min(lb.z, atomi.z - rad);
    ub.x = max(ub.x, atomi.x + rad);
    ub.y = max(ub.y, atomi.y + rad);
    ub.z = max(ub.z, atomi.z + rad);
  }
  lb.x -= margin;
  lb.y -= margin;
  lb.z -= margin;
  ub.x += margin;
  ub.y += margin;
  ub.z += margin;
  return;
}

Box getHostGridFromBox(Point3DDouble lb, Point3DDouble ub, double res) {
  Box retval;
  retval.lowestPoint = lb;
  Point3DInteger extent;
  Point3DDouble dims = ub - lb;
  // printf("Box dimensions: %fx%fx%f\n",dims.x,dims.y,dims.z);
  extent.x = (int)(ceil(dims.x / res));
  extent.y = (int)(ceil(dims.y / res));
  extent.z = (int)(ceil(dims.z / res));
  retval.res = res;
  retval.extent = extent;
  // printf("Box extent:
  // %ux%ux%u\n",retval.extent.x,retval.extent.y,retval.extent.z);
  retval.points = new float[extent.x * extent.y * extent.z];
  return retval;
}

Box VolumeOverlap::getHostGrid(Molecule &mol, double res, double margin) {
  // Calculate bounding box
  Point3DDouble ub, lb;
  boundingBox(mol, margin, lb, ub);
  return getHostGridFromBox(lb, ub, res);
  // printf("Got bounding box: [%f,%f,%f] -
  // [%f,%f,%f]\n",retval.lb.x,retval.lb.y,retval.lb.z,ub.x,ub.y,ub.z);
}

double VolumeOverlap::tanimoto(Molecule &fixed, Molecule &variable, Box &box) {
  double oAA = getOverlapVolumeGaussianNonAnalytic(fixed, fixed, box);
  double oBB = getOverlapVolumeGaussianNonAnalytic(variable, variable, box);
  double oAB = getOverlapVolumeGaussianNonAnalytic(fixed, variable, box);
  double tanimoto = oAB / (oAA + oBB - oAB);
  return tanimoto;
}

// Esto esta copiado de PAPERS y no creo que sirva. Utiliza otra forma de
// calcular el solapamiento entre dos atomos
/*
 double VolumeOverlap::calculateOverlapVolumeAtomsGaussAnalytic(Atom* atom1,
 Atom* atom2) {
 const double partialalpha = 2.41798793102;
 const double pi = 3.14159265358;
 double overlap = 0.0;
 double p = 2.7;			//2 * sqrt(2);

 double alphai = partialalpha
 / (atom1->radiusVanDerWaals * atom1->radiusVanDerWaals);

 double Rij2 = distsqr(*atom1, *atom2);

 //if (Rij2 > 16) continue;

 double alphaj = partialalpha
 / (atom2->radiusVanDerWaals * atom2->radiusVanDerWaals);
 double Kij = exp(-(alphai * alphaj * Rij2) / (alphai + alphaj));
 double Vij = p * p * Kij * pow((pi / (alphai + alphaj)), 1.5);

 overlap += Vij;

 return overlap;
 }*/
double VolumeOverlap::calculateOverlapVolumeAtomsGaussAnalytic(Atom *atom1,
                                                               Atom *atom2) {

  const double pi = 3.14159265358;
  // double overlap = 0.0;

  double p = 2.828427;
  // double p = 2.7;
  Atom atomFixed = *atom1;
  double alphai =
      pow(3 * p * sqrt(pi) /
              (4 * (atomFixed.radiusVanDerWaals * atomFixed.radiusVanDerWaals *
                    atomFixed.radiusVanDerWaals)),
          2 / 3);

  Atom atomVariable = *atom2;
  double Rij2 = distsqr(atomFixed, atomVariable);

  // if (Rij2 > 16) continue;

  double alphaj = pow(3 * p * sqrt(pi) /
                          (4 * (atomVariable.radiusVanDerWaals *
                                atomVariable.radiusVanDerWaals *
                                atomVariable.radiusVanDerWaals)),
                      2 / 3);

  double Kij = exp(-(alphai * alphaj * Rij2) / (alphai + alphaj));
  // p = 2.7 instead of 2sqrt(2)
  //  double Vij = 7.29 * Kij * pow((pi / (alphai + alphaj)), 1.5);
  double Vij = p * p * Kij * pow((pi / (alphai + alphaj)), 1.5);

  return Vij;
}
double VolumeOverlap::calculateWeight(Atom *atom, vector<Atom> *atoms,
                                      int iAtom) {
  double pi = 3.14159265358;
  double vi = 4 * pi * pow(atom->radiusVanDerWaals, 3) / 3;
  // double p = 2.828427;
  // double p = 2.7;
  // double vi = p * (pow((pi / 0.8665), 1.5));

  double vij = 0;
  double k = 0.8665;
  for (unsigned i = 0; i < atoms->size(); i++) {
    if (i == iAtom)
      continue;
    vij += calculateOverlapVolumeAtomsGaussAnalytic(atom, &atoms->at(i));
  }

  double result = vi / (vi + k * vij);

  return result;
}

double
VolumeOverlap::calculateOverlapVolumeAtomsWEGA(Atom *atom1, Atom *atom2,
                                               bool sameVanDerWaalsRadius) {

  const double pi = 3.14159265358;
  // double overlap = 0.0;
  // double p = 2.8284;
  // double p = 2.7;
  Atom atomFixed = *atom1;
  Atom atomVariable = *atom2;
  double Vij = 0;
  double Rij2 = distsqr(atomFixed, atomVariable);
  if (sameVanDerWaalsRadius) {
    double Kij = exp(Rij2 * -0.3731438999881213);
    Vij = 7.99984656 * Kij * pow(2.104813085306902, 1.5);
  } else {
    double alphai = 2.417972471923026 /
                    (atomFixed.radiusVanDerWaals * atomFixed.radiusVanDerWaals);
    double alphaj = 2.417972471923026 / (atomVariable.radiusVanDerWaals *
                                         atomVariable.radiusVanDerWaals);
    /*double alphai = pow(
     3 * p * sqrt(pi)
     / (4
     * (atomFixed.radiusVanDerWaals
     * atomFixed.radiusVanDerWaals
     * atomFixed.radiusVanDerWaals)), 2 / 3);
     double alphaj = pow(
     3 * p * sqrt(pi)
     / (4
     * (atomVariable.radiusVanDerWaals
     * atomVariable.radiusVanDerWaals
     * atomVariable.radiusVanDerWaals)),
     2 / 3);
     */

    double Kij = exp(-(alphai * alphaj * Rij2) / (alphai + alphaj));
    Vij = 7.99984656 * Kij * pow((pi / (alphai + alphaj)), 1.5);
  }

  return Vij;
}

/**Fast implementation exp
 * https://stackoverflow.com/questions/10552280/fast-exp-calculation-possible-to-improve-accuracy-without-losing-too-much-perfo
 *
 */

inline double fastExp(double x) {
  if (x < -80) {
    x = -80;
  }
  volatile union {
    float f;
    unsigned int i;
  } cvt;

  /* exp(x) = 2^i * 2^f; i = floor (log2(e) * x), 0 <= f <= 1 */
  float t = x * 1.442695041f;
  float fi = floorf(t);
  float f = t - fi;
  int i = (int)fi;
  cvt.f =
      (0.3371894346f * f + 0.657636276f) * f + 1.00172476f; /* compute 2^f */
  cvt.i += (i << 23);                                       /* scale by 2^i */
  return cvt.f;
}

inline double VolumeOverlap::calculateOverlapVolumeAtomsWEGA(
    double atom1X, double atom1Y, double atom1Z, double atom2X, double atom2Y,
    double atom2Z, bool sameVanDerWaalsRadius, double radiusFixed,
    double radiusVariable) {

  // double p = 2.8284;
  // double p = 2.7;

  double Vij = 0;
  double tempX = atom1X - atom2X;
  double tempY = atom1Y - atom2Y;
  double tempZ = atom1Z - atom2Z;
  double Rij2 = tempX * tempX + tempY * tempY + tempZ * tempZ;

  if (sameVanDerWaalsRadius) {
    //		double Kij2 = exp(Rij2 * -0.3731438999881213);
    double Kij = fastExp(Rij2 * -0.3731438999881213);

    // Vij = 7.99984656 * Kij * pow(2.104813085306902, 1.5);
    Vij = 24.428790199 * Kij;
  } else {
    double pi = 3.14159265358;
    double alphai = 2.417972471923026 / (radiusFixed * radiusFixed);
    double alphaj = 2.417972471923026 / (radiusVariable * radiusVariable);

    double Kij = fastExp(-(alphai * alphaj * Rij2) / (alphai + alphaj));
    Vij = 7.99984656 * Kij * pow((pi / (alphai + alphaj)), 1.5);
  }

  return Vij;
}

inline double VolumeOverlap::preciseCalculateOverlapVolumeAtomsWEGA(
    double atom1X, double atom1Y, double atom1Z, double atom2X, double atom2Y,
    double atom2Z, bool sameVanDerWaalsRadius, double radiusFixed,
    double radiusVariable) {

  double Vij = 0;
  double tempX = atom1X - atom2X;
  double tempY = atom1Y - atom2Y;
  double tempZ = atom1Z - atom2Z;
  double Rij2 = tempX * tempX + tempY * tempY + tempZ * tempZ;

  if (sameVanDerWaalsRadius) {
    double Kij = exp(Rij2 * -0.3731438999881213);
    Vij = 24.428790199 * Kij;
  } else {
    double alphai = 2.417972471923026 / (radiusFixed * radiusFixed);
    double alphaj = 2.417972471923026 / (radiusVariable * radiusVariable);

    double Kij = exp(-(alphai * alphaj * Rij2) / (alphai + alphaj));
    Vij = 7.99984656 * Kij * pow((M_PI / (alphai + alphaj)), 1.5);
  }

  return Vij;
}

double VolumeOverlap::calculateWeightWEGA(Atom *atom, vector<Atom> *atoms,
                                          int iAtom,
                                          bool sameVanDerWaalsRadius) {
  double pi = 3.14159265358;
  double vi = 4 * pi * pow(atom->radiusVanDerWaals, 3) / 3;
  // double p = 2.8284;

  double vij = 0;
  double k = 0.8665;
  for (unsigned i = 0; i < atoms->size(); i++) {
    if (i == iAtom)
      continue;
    vij += calculateOverlapVolumeAtomsWEGA(atom, &atoms->at(i),
                                           sameVanDerWaalsRadius);
  }

  double result = vi / (vi + k * vij);

  return result;
}

double VolumeOverlap::calculateWeightWEGA(double atomX, double atomY,
                                          double atomZ, double sizeAtoms,
                                          double *atoms, double *radius,
                                          int iAtom,
                                          bool sameVanDerWaalsRadius) {

  double pi = 3.14159265358;
  double vi = 4 * pi * pow(radius[iAtom], 3) / 3;
  // double p = 2.8284;

  double vij = 0;
  double k = 0.8665;
  int iterador = 0;
  for (unsigned i = 0; i < sizeAtoms; i++) {
    iterador = 3 * i;
    if (i == iAtom)
      continue;

    vij += calculateOverlapVolumeAtomsWEGA(
        atomX, atomY, atomZ, atoms[iterador], atoms[iterador + 1],
        atoms[iterador + 2], sameVanDerWaalsRadius, radius[iAtom], radius[i]);
    // cout << "VIJ de : " <<iAtom<<": "<< vij<<endl;
  }
  double result = vi / (vi + k * vij);

  return result;
}

void VolumeOverlap::longestShortestAxis(Molecule &mol, double *longest,
                                        double *shortest) {
  for (unsigned i = 0; i < mol.atoms.size(); i++) {
    double modulo =
        sqrt(mol.atoms[i].x * mol.atoms[i].x + mol.atoms[i].y * mol.atoms[i].y +
             mol.atoms[i].z * mol.atoms[i].z);
    if (modulo > *longest) {
      *longest = modulo;
    } else if (modulo < *shortest) {
      *shortest = modulo;
    }
  }
}
void VolumeOverlap::boundingBox(Molecule &mol, Point3DDouble *lb,
                                Point3DDouble *ub) {
  Atom &atom1 = mol.atoms[0];
  lb->x = atom1.x;
  lb->y = atom1.y;
  lb->z = atom1.z;
  ub->x = atom1.x;
  ub->y = atom1.y;
  ub->z = atom1.z;
  for (unsigned i = 1; i < mol.atoms.size(); i++) {
    Atom &atomi = mol.atoms[i];
    lb->x = min(lb->x, atomi.x);
    lb->y = min(lb->y, atomi.y);
    lb->z = min(lb->z, atomi.z);
    ub->x = max(ub->x, atomi.x);
    ub->y = max(ub->y, atomi.y);
    ub->z = max(ub->z, atomi.z);
  }
}
