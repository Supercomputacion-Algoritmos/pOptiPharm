/*
 * Parameters.h
 *
 *  Created on: 16 jan. 2017
 *      Author: Savins
 */

#ifndef SRC_FUNCTIONS_PARAMETERS_H_
#define SRC_FUNCTIONS_PARAMETERS_H_

#include "functions/VolumeOverlap.h"
#include "model/Molecule.h"
#include "model/Point3DDouble.h"
#include "uego/uego.h"
#include "uego/usrintf.h"
#include <cmath> // Para abs(double)
#include <string.h>
#include <string>

// Align molecules according the three axis.
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SVD>
using namespace Eigen;

class Parameters {
public:
  struct Configuration {
    bool updateConfig = false;
    bool onlyConfig = true;
    bool align = true;
    bool pointOneZero = false;
    bool saveMolFile = false;
    bool considerHidrogens = false;
    bool sameRadiusVdW = true;
    uint32_t threads_num = 0;
    bool pin_threads = false;
    bool force_CUDA = false;
    uint32_t min_combined_size = 500;
  };
  int static readArguments(int argc, char **argv, int *indexArray,
                           int lenghtIndexArray, Configuration *c);
  string static getNameMol(char *path);
  void static defineDeltaMovement(Ini *ini, double *minX, double *maxX,
                                  double *minY, double *maxY, double *minZ,
                                  double *maxZ);
  void static alignMolecule(Molecule *mol);
};

#endif /* SRC_FUNCTIONS_TANIMOTO_H_ */
