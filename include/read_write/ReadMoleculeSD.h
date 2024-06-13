/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   ReadMoleculeSD.h
 * Author: savins
 *
 * Created on 25 de enero de 2016, 18:52
 */

#ifndef READMOLECULESD_H
#define READMOLECULESD_H
using namespace std;
#include "model/Molecule.h"
#include <fstream>
#include <iostream>
class ReadMoleculeSD {
public:
public:
  void static readMol(string filePath, vector<Molecule> *molecules);
};

#endif /* READMOLECULESD_H */
