/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   readMolecule.cpp
 * Author: savins
 *
 * Created on 2 de enero de 2016, 10:50
 */

#include "read_write/ReadMolecule.h"
#include <chrono>
#include <regex>

std::istream &safeGetline(std::istream &is, std::string &t) {
  t.clear();

  // The characters in the stream are read one-by-one using a std::streambuf.
  // That is faster than reading them one-by-one using the std::istream.
  // Code that uses streambuf this way must be guarded by a sentry object.
  // The sentry object performs various tasks,
  // such as thread synchronization and updating the stream state.

  std::istream::sentry se(is, true);
  std::streambuf *sb = is.rdbuf();

  for (;;) {
    int c = sb->sbumpc();
    switch (c) {
    case '\n':
      return is;
    case '\r':
      if (sb->sgetc() == '\n')
        sb->sbumpc();
      return is;
    case EOF:
      // Also handle the case when the last line has no line ending
      if (t.empty())
        is.setstate(std::ios::eofbit);
      return is;
    default:
      t += (char)c;
    }
  }
}

/*
 void ReadMolecule::readMol(char* filePath, Molecule *molecule,
 bool considerHydrogens, bool sameVanDerWaalsRadius) {
 auto begin = std::chrono::high_resolution_clock::now();

 ifstream infile(filePath);
 if (!infile.is_open()) {
 cerr << "uego: !! The file " << filePath << " can not be open" << endl;
 exit(EXIT_FAILURE);
 }
 string line = "";
 while (safeGetline(infile, line)) {
 if (line.compare("@<TRIPOS>MOLECULE") == 0) {
 //Name
 safeGetline(infile, molecule->mol_name);
 //Number variables
 safeGetline(infile, line);
 istringstream iss(line);
 iss >> molecule->num_atoms >> molecule->num_bonds
 >> molecule->num_subst >> molecule->num_feat
 >> molecule->num_sets;
 //
 safeGetline(infile, molecule->mol_type);
 safeGetline(infile, molecule->charge_type);
 safeGetline(infile, line);
 while (line.compare("@<TRIPOS>ATOM") != 0) {
 molecule->mol_comment += line;
 molecule->mol_comment += "\n";
 safeGetline(infile, line);
 }
 }

 if (line.compare("@<TRIPOS>ATOM") == 0) {
 for (long i = 0; i < molecule->num_atoms; ++i) {
 Atom *atom = new Atom();
 safeGetline(infile, line);
 istringstream iss(line);
 iss >> atom->atom_id >> atom->atom_name >> atom->x >> atom->y
 >> atom->z >> atom->atom_type >> atom->subst_id
 >> atom->subst_name >> atom->charge >> atom->status_bit;

 if (!considerHydrogens && atom->atom_name[0] == 'H') {
 if (atom->atom_name.size() > 1) {
 if (isdigit(atom->atom_name[1])) {
 delete atom;
 continue;
 }
 } else {
 delete atom;
 continue;
 }

 }

 if (sameVanDerWaalsRadius) {
 atom->radiusVanDerWaals = 1.80;

 } else {
 atom->radiusVanDerWaals = atom->setRadiusVanDerWaals(
 atom->atom_type);
 }

 molecule->atoms.push_back(*atom);
 delete atom;
 }
 }

 if (line.compare("@<TRIPOS>BOND") == 0) {
 for (long i = 0; i < molecule->num_bonds; ++i) {
 Bond *bond = new Bond();
 safeGetline(infile, line);
 istringstream iss(line);
 iss >> bond->bond_id >> bond->origin_atom_id
 >> bond->target_atom_id >> bond->bond_type
 >> bond->status_bits;

 if (!considerHydrogens) {
 bool encontrado = false;

 for (unsigned i = 0; i < molecule->atoms.size(); ++i) {
 if (molecule->atoms[i].atom_id == bond->target_atom_id)
 encontrado = true;
 }
 if (!encontrado) {
 delete bond;
 continue;
 }
 }
 molecule->bonds.push_back(*bond);
 delete bond;
 }
 }
 }
 //Update numbers of atoms and bonds
 molecule->num_atoms = molecule->atoms.size();
 molecule->num_bonds = molecule->bonds.size();

 //Reasign id to atoms
 for (int i = 0; i < molecule->num_atoms; i++) {
 molecule->atoms[i].old_atom_id = molecule->atoms[i].atom_id;
 molecule->atoms[i].atom_id = (i + 1);
 }
 //Reasign id to bonds
 for (int i = 0; i < molecule->num_bonds; i++) {
 molecule->bonds[i].bond_id = (i + 1);
 for (int j = 0; j < molecule->num_atoms; j++) {

 if (molecule->bonds[i].origin_atom_id
 == molecule->atoms[j].old_atom_id) {
 molecule->bonds[i].origin_atom_id = molecule->atoms[j].atom_id;
 }
 if (molecule->bonds[i].target_atom_id
 == molecule->atoms[j].old_atom_id) {
 molecule->bonds[i].target_atom_id = molecule->atoms[j].atom_id;
 }
 }
 }

 infile.close();

 auto end = std::chrono::high_resolution_clock::now();
 std::cout << "Tiempo leido: " << std::chrono::duration_cast
 < std::chrono::nanoseconds
 > (end - begin).count() << "ns" << std::endl;
 cerr << "uego: -- File " << filePath << " read successfully" << endl;

 }
 */
void ReadMolecule::readMol(char *filePath, Molecule *molecule,
                           bool considerHydrogens, bool sameVanDerWaalsRadius) {

  ifstream infile(filePath);
  if (!infile.is_open()) {
    cerr << "optipharm: !! The file " << filePath << " can not be open\n";
    exit(EXIT_FAILURE);
  }
  string line = "";
  while (safeGetline(infile, line)) {
    if (line.compare("@<TRIPOS>MOLECULE") == 0) {
      // Name
      safeGetline(infile, molecule->mol_name);
      // Number variables
      safeGetline(infile, line);
      istringstream iss(line);
      iss >> molecule->num_atoms >> molecule->num_bonds >>
          molecule->num_subst >> molecule->num_feat >> molecule->num_sets;
      //
      safeGetline(infile, molecule->mol_type);
      safeGetline(infile, molecule->charge_type);
      safeGetline(infile, line);
      while (line.compare("@<TRIPOS>ATOM") != 0) {
        molecule->mol_comment += line;
        molecule->mol_comment += "\n";
        safeGetline(infile, line);
      }
    }

    if (line.compare("@<TRIPOS>ATOM") == 0) {
      for (long i = 0; i < molecule->num_atoms; ++i) {
        Atom *atom = new Atom();
        safeGetline(infile, line);
        istringstream iss(line);
        iss >> atom->atom_id >> atom->atom_name >> atom->x >> atom->y >>
            atom->z >> atom->atom_type >> atom->subst_id >> atom->subst_name >>
            atom->charge >> atom->status_bit;
        atom->radiusVanDerWaals = 1.80;
        molecule->atoms.push_back(*atom);
        delete atom;
      }
    }

    if (line.compare("@<TRIPOS>BOND") == 0) {
      for (long i = 0; i < molecule->num_bonds; ++i) {
        Bond *bond = new Bond();
        safeGetline(infile, line);
        istringstream iss(line);
        iss >> bond->bond_id >> bond->origin_atom_id >> bond->target_atom_id >>
            bond->bond_type >> bond->status_bits;

        molecule->bonds.push_back(*bond);
        delete bond;
      }
    }
  }
  infile.close();

  if (!considerHydrogens) {
    // arary to save which atoms remove.
    bool *atomsToDelete = new bool[molecule->num_atoms];
    for (unsigned i = 0; i < molecule->num_atoms; i++) {
      atomsToDelete[i] = false;
    }

    // check atoms to delete
    for (int i = molecule->num_atoms - 1; i > -1; i--) {
      if (molecule->atoms[i].atom_name[0] == 'H') {
        if (molecule->atoms[i].atom_name.size() > 1) {
          if (isdigit(molecule->atoms[i].atom_name[1])) {
            atomsToDelete[i] = true;
            molecule->atoms.erase(molecule->atoms.begin() + i);
          }
        } else {
          atomsToDelete[i] = true;
          molecule->atoms.erase(molecule->atoms.begin() + i);
        }
      }
    }

    // New ID of each bonds after the removing.
    unsigned *idminusHydrogen = new unsigned[molecule->num_atoms];
    unsigned contador = 0;
    for (unsigned i = 0; i < molecule->num_atoms; i++) {
      if (atomsToDelete[i]) {
        contador += 1;
      }
      idminusHydrogen[i] = contador;
    }

    // remove bonds where hydrogen are in some side of the bond
    for (int i = (molecule->num_bonds - 1); i > -1; i--) {
      if (atomsToDelete[molecule->bonds[i].origin_atom_id - 1] ||
          atomsToDelete[molecule->bonds[i].target_atom_id - 1]) {
        molecule->bonds.erase(molecule->bonds.begin() + i);
      }
    }
    // updating number of bonds after removing hydrogen.
    molecule->num_bonds = molecule->bonds.size();

    // cout << molecule->toString();
    // Reasign id to bonds
    //  To do so, for each id is compared the number of hydrogens removed until
    //  the id. Thus, that number of hydrogens is remove from the original id.
    // cout << "Renaming"<<endl;
    for (unsigned i = 0; i < molecule->num_bonds; i++) {
      // cout << molecule->bonds[i].bond_id<< " "
      // <<molecule->bonds[i].origin_atom_id <<" "<<
      // molecule->bonds[i].target_atom_id<<endl;
      molecule->bonds[i].bond_id = i + 1;
      molecule->bonds[i].origin_atom_id -=
          idminusHydrogen[molecule->bonds[i].origin_atom_id - 1];
      molecule->bonds[i].target_atom_id -=
          idminusHydrogen[molecule->bonds[i].target_atom_id - 1];
      // cout << molecule->bonds[i].bond_id<< " "
      // <<molecule->bonds[i].origin_atom_id << " "
      // <<molecule->bonds[i].target_atom_id<<endl; cout << "__________"<<endl;
    }

    molecule->num_atoms = molecule->atoms.size();
    for (unsigned i = 0; i < molecule->num_atoms; i++) {
      molecule->atoms[i].atom_id = i + 1;
    }
    delete atomsToDelete;
    delete idminusHydrogen;
  }
  if (!sameVanDerWaalsRadius) {
    for (int i = 0; i < molecule->num_atoms; i++) {
      molecule->atoms[i].radiusVanDerWaals =
          molecule->atoms[i].setRadiusVanDerWaals(molecule->atoms[i].atom_type);
    }
  }

  char msg[250];
  sprintf(msg, "File  %s read successfully", filePath);
  message(msg, MSG_INFORMATION);
}

void ReadMolecule::readDBMol2(char *filePath, vector<Molecule> *conformations,
                              bool considerHydrogens,
                              bool sameVanDerWaalsRadius) {

  ifstream infile(filePath);
  if (!infile.is_open()) {
    cerr << "uego: !! The file " << filePath << " can not be open\n";
    exit(EXIT_FAILURE);
  }

  Molecule *molecule;
  while (!infile.eof()) {

    string line = "";
    while (safeGetline(infile, line)) {
      if (line.compare("@<TRIPOS>MOLECULE") == 0) {
        molecule = new Molecule();
        // Name
        safeGetline(infile, molecule->mol_name);
        // Number variables
        safeGetline(infile, line);
        istringstream iss(line);
        iss >> molecule->num_atoms >> molecule->num_bonds >>
            molecule->num_subst >> molecule->num_feat >> molecule->num_sets;
        //
        safeGetline(infile, molecule->mol_type);
        safeGetline(infile, molecule->charge_type);
        safeGetline(infile, line);
        while (line.compare("@<TRIPOS>ATOM") != 0) {
          molecule->mol_comment += line + "\n";
          safeGetline(infile, line);
        }
      }

      if (line.compare("@<TRIPOS>ATOM") == 0) {

        for (long i = 0; i < molecule->num_atoms; ++i) {
          Atom atom;
          safeGetline(infile, line);
          istringstream iss(line);
          iss >> atom.atom_id >> atom.atom_name >> atom.x >> atom.y >> atom.z >>
              atom.atom_type >> atom.subst_id >> atom.subst_name >>
              atom.charge >> atom.status_bit;

          if (!considerHydrogens && atom.atom_name[0] == 'H') {
            if (atom.atom_name.size() > 1) {
              if (isdigit(atom.atom_name[1])) {
                continue;
              }
            } else {
              continue;
            }
          }

          if (sameVanDerWaalsRadius) {
            atom.radiusVanDerWaals = 1.80;

          } else {
            atom.radiusVanDerWaals = atom.setRadiusVanDerWaals(atom.atom_type);
          }
          molecule->atoms.push_back(atom);
        }
      }

      if (line.compare("@<TRIPOS>BOND") == 0) {
        for (long i = 0; i < molecule->num_bonds; ++i) {
          Bond bond;
          safeGetline(infile, line);
          istringstream iss(line);
          iss >> bond.bond_id >> bond.origin_atom_id >> bond.target_atom_id >>
              bond.bond_type >> bond.status_bits;
          if (!considerHydrogens) {
            bool encontrado = false;

            for (unsigned i = 0; i < molecule->atoms.size(); ++i) {
              if (molecule->atoms[i].atom_id == bond.target_atom_id)
                encontrado = true;
            }
            if (!encontrado) {
              continue;
            }
          }
          molecule->bonds.push_back(bond);
        }
        // Update numbers of atoms and bonds
        molecule->num_atoms = molecule->atoms.size();
        molecule->num_bonds = molecule->bonds.size();

        // Reasign id to bonds
        for (unsigned i = 0; i < molecule->num_atoms; i++) {
          for (unsigned j = 0; j < molecule->num_bonds; j++) {
            // replace origin bond
            if (molecule->bonds[j].origin_atom_id ==
                molecule->atoms[i].atom_id) {
              molecule->bonds[j].origin_atom_id = i + 1;
            }
            // replace target bond
            if (molecule->bonds[j].target_atom_id ==
                molecule->atoms[i].atom_id) {
              molecule->bonds[j].target_atom_id = i + 1;
            }
          }
          molecule->atoms[i].atom_id = i + 1;
        }
        for (unsigned j = 0; j < molecule->num_bonds; j++) {
          molecule->bonds[j].bond_id = j + 1;
        }
        conformations->push_back(*molecule);
        delete molecule;
      }
    }
  }
  infile.close();

  char msg[250];
  sprintf(msg, "DBmolecule  %s read successfully", filePath);
  message(msg, MSG_INFORMATION);
}

/*
 void ReadMolecule::readConfCyndi(char* filePath,
 vector<Molecule> *conformations, bool considerHydrogens, bool
 sameVanDerWaalsRadius) { Molecule molecule; std::string stringPath(filePath);
 int i = 0;
 ifstream infile(filePath);
 if (!infile.is_open()) {
 cerr << "uego: !! The file " << filePath << " can not be open" << endl;
 exit (EXIT_FAILURE);
 }
 while (!infile.eof()) {
 cout <<"contador: "<< i++<<endl;
 string line = "";
 while (safeGetline(infile, line)) {
 if (line.compare("@<TRIPOS>MOLECULE") == 0) {
 //Name
 safeGetline(infile, molecule.mol_name);
 //Number variables
 safeGetline(infile, line);
 istringstream iss(line);
 iss >> molecule.num_atoms >> molecule.num_bonds
 >> molecule.num_subst >> molecule.num_feat
 >> molecule.num_sets;
 //
 safeGetline(infile, molecule.mol_type);
 safeGetline(infile, molecule.charge_type);
 safeGetline(infile, line);
 while (line.compare("@<TRIPOS>ATOM") != 0) {
 molecule.mol_comment += line;
 molecule.mol_comment += "\n";
 safeGetline(infile, line);
 }
 }

 if (line.compare("@<TRIPOS>ATOM") == 0) {

 for (long i = 0; i < molecule.num_atoms; ++i) {
 Atom atom;
 safeGetline(infile, line);
 istringstream iss(line);
 iss >> atom.atom_id >> atom.atom_name >> atom.x >> atom.y
 >> atom.z >> atom.atom_type >> atom.subst_id
 >> atom.subst_name >> atom.charge
 >> atom.status_bit;

 if (!considerHydrogens && atom.atom_name[0] == 'H') {
 if (atom.atom_name.size() > 1) {
 if (isdigit(atom.atom_name[1])) {
 continue;
 }
 } else {
 continue;
 }

 }

 if (sameVanDerWaalsRadius) {
 atom.radiusVanDerWaals = 1.70;

 } else {
 atom.radiusVanDerWaals = atom.setRadiusVanDerWaals(
 atom.atom_type);
 }

 molecule.atoms.push_back(atom);
 }
 }

 if (line.compare("@<TRIPOS>BOND") == 0) {
 for (long i = 0; i < molecule.num_bonds; ++i) {
 Bond bond;
 ;
 safeGetline(infile, line);
 istringstream iss(line);
 iss >> bond.bond_id >> bond.origin_atom_id
 >> bond.target_atom_id >> bond.bond_type
 >> bond.status_bits;
 if (!considerHydrogens) {
 bool encontrado = false;

 for (unsigned i = 0; i < molecule.atoms.size(); ++i) {
 if (molecule.atoms[i].atom_id
 == bond.target_atom_id)
 encontrado = true;
 }
 if (!encontrado) {
 continue;
 }
 }
 molecule.bonds.push_back(bond);
 }
 }
 }
 //Update numbers of atoms and bonds
 molecule.num_atoms = molecule.atoms.size();
 molecule.num_bonds = molecule.bonds.size();

 //Reasign id to atoms
 for (int i = 0; i < molecule.num_atoms; i++) {
 molecule.atoms[i].old_atom_id = molecule.atoms[i].atom_id;
 molecule.atoms[i].atom_id = (i + 1);
 }
 //Reasign id to bonds
 for (int i = 0; i < molecule.num_bonds; i++) {
 molecule.bonds[i].bond_id = (i + 1);
 for (int j = 0; j < molecule.num_atoms; j++) {

 if (molecule.bonds[i].origin_atom_id
 == molecule.atoms[j].old_atom_id) {
 molecule.bonds[i].origin_atom_id =
 molecule.atoms[j].atom_id;
 }
 if (molecule.bonds[i].target_atom_id
 == molecule.atoms[j].old_atom_id) {
 molecule.bonds[i].target_atom_id =
 molecule.atoms[j].atom_id;
 }
 }
 }
 cout << molecule.mol_name;
 conformations->push_back(molecule);
 }
 infile.close();

 cerr << "uego: -- DBmolecule " << filePath << " read successfully" << endl;
 }
 */
