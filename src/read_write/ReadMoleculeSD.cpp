
/*
 * File:   ReadMoleculeSD.cpp
 * Author: savins
 *
 * Created on 25 de enero de 2016, 18:52
 */

#include "read_write/ReadMoleculeSD.h"

void ReadMoleculeSD::readMol(string filePath, vector<Molecule> *molecules) {
  string line = "";
  ifstream infile;

  infile.open(filePath.c_str());

  while (!infile.eof()) {

    Molecule mol;

    getline(infile, mol.mol_name);

    getline(infile, line); // information
    getline(infile, line); // blank
    getline(infile, line);
    istringstream iss(line);
    iss >> mol.num_atoms >> mol.num_bonds;
    for (int i = 0; i < mol.num_atoms; i++) {
      Atom atom;
      getline(infile, line);

      istringstream iss(line);
      string vdv;
      iss >> atom.x >> atom.y >> atom.z >> vdv;
      atom.radiusVanDerWaals = atom.setRadiusVanDerWaals(vdv);
      mol.atoms.push_back(atom);
    }
    for (int i = 0; i < mol.num_bonds; i++) {
      getline(infile, line);
    }
    getline(infile, line);

    getline(infile, line);
    molecules->push_back(mol);

    if (line.compare("$$$$") != 0)
      break;
  }

  infile.close();
  cout << "Read file completed!!" << endl;
}
