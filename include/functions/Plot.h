/*
 * Plot.h
 *
 *  Created on: 26 sept. 2018
 *      Author: savins
 */

#ifndef SRC_FUNCTIONS_PLOT_H_
#define SRC_FUNCTIONS_PLOT_H_

#include "functions/MoveAndRotate.h"
#include "model/Molecule.h"
#include "model/Point3DDouble.h"
#include "read_write/WriteMolecule.h"

class Plot {
public:
  Plot();
  virtual ~Plot();
  void static MoleculesSolution(Molecule *target, Molecule *query, double angle,
                                double x1, double y1, double z1, double x2,
                                double y2, double z2, double deltaX,
                                double deltaY, double deltaZ, double value,
                                string name);
};

#endif /* SRC_FUNCTIONS_PLOT_H_ */
