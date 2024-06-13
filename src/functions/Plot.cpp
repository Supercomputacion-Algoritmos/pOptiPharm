/*
 * Plot.cpp
 *
 *  Created on: 26 sept. 2018
 *      Author: savins
 */

#include "functions/Plot.h"

Plot::Plot() {
  // TODO Auto-generated constructor stub
}

Plot::~Plot() {
  // TODO Auto-generated destructor stub
}

void Plot::MoleculesSolution(Molecule *query, Molecule *variable, double angle,
                             double x1, double y1, double z1, double x2,
                             double y2, double z2, double deltaX, double deltaY,
                             double deltaZ, double value, string name) {
  cout << "Pintada " << query->mol_name << " and " << variable->mol_name
       << " en " << name << "\n";
  Molecule *moleculeRotated = new Molecule();
  Point3DDouble p1, p2;
  p1.x = x1;
  p1.y = y1;
  p1.z = z1;
  p2.x = x2;
  p2.y = y2;
  p2.z = z2;

  MoveAndRotate::RotateMolAccording1Axis(variable, angle, p1, p2,
                                         moleculeRotated);
  MoveAndRotate::MolToNewPosition(moleculeRotated, deltaX, deltaY, deltaZ);

  WriteMolecule::writeMol(name + ".mol2", moleculeRotated);
  string finalname = query->mol_name + "-" + variable->mol_name + "-" + name;
  string command = "cat " + query->name_file + ".mol2 " + name + ".mol2 > " +
                   finalname + ".mol2";
  cout << command << endl;
  system(command.c_str());
  command = "python scriptPlot.py " + finalname + ".mol2 " + finalname + ".pse";
  system(command.c_str());

  command = "python label.py " + finalname + ".pse.png " + to_string(value);
  system(command.c_str());
}
