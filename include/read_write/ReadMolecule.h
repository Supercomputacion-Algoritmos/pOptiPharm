
/*
 * File:   readMolecule.h
 * Author: savins
 *
 * Created on 2 de enero de 2016, 10:50
 */

#ifndef READMOLECULE_H
#define READMOLECULE_H

#include "model/Atom.h"
#include "model/Bond.h"
#include "model/Molecule.h"
#include "uego/uego.h"
#include <ctype.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h> /* exit, EXIT_FAILURE */
#include <string>
#include <vector>
using namespace std;

class ReadMolecule {
public:
  void static readMol(char *filePath, Molecule *molecule,
                      bool considerHydrogens, bool sameVanDerWaalsRadius);

  void static readDBMol2(char *filePath, vector<Molecule> *conformations,
                         bool considerHydrogens, bool sameVanDerWaalsRadius);
};

#endif /* READMOLECULE_H */
