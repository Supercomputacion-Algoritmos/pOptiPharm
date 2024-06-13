/*
 * WriteMolecule.h
 *
 *  Created on: 10 de abr. de 2016
 *      Author: savins
 */

#ifndef SRC_READ_WRITE_WRITEMOLECULE_H_
#define SRC_READ_WRITE_WRITEMOLECULE_H_
#include "model/Atom.h"
#include "model/Bond.h"
#include "model/Molecule.h"
#include <fstream>
#include <iostream>
#include <vector>
using namespace std;

class WriteMolecule {
public:
  void static writeMol(string filePath, Molecule *molecule);
  string static MolToString(Molecule *molecule);
};

#endif /* SRC_READ_WRITE_WRITEMOLECULE_H_ */
