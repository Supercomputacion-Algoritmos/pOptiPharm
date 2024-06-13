/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   Molecule.cpp
 * Author: Savíns
 *
 * Created on 25 de noviembre de 2015, 16:26
 */

#include "model/Molecule.h"
#include "functions/VolumeOverlap.h"

Molecule::Molecule() {
  this->name_file = "";
  this->mol_name = "";
  this->mol_type = "";
  this->charge_type = "";
  this->mol_comment = "";
  this->num_atoms = 0;
  this->num_bonds = 0;
  this->num_subst = 0;
  this->num_feat = 0;
  this->num_sets = 0;
  this->tanimoto = 0;
  this->atomsXYZ = NULL;
  this->radiusAtoms = NULL;
  this->weightAtoms = NULL;
}

Molecule::Molecule(string mol_name, string mol_type, string charge_type,
                   string mol_comment, int num_atoms, int num_bonds,
                   int num_subst, int num_feat, int num_sets) {
  this->mol_name = mol_name;
  this->mol_type = mol_type;
  this->charge_type = charge_type;
  this->mol_comment = mol_comment;
  this->num_atoms = num_atoms;
  this->num_bonds = num_bonds;
  this->num_subst = num_subst;
  this->num_feat = num_feat;
  this->num_sets = num_sets;
  this->atomsXYZ = NULL;
  this->radiusAtoms = NULL;
  this->weightAtoms = NULL;
  this->tanimoto = 0;
}

Molecule::Molecule(const Molecule &mol) {
  this->name_file = mol.name_file;
  this->mol_name = mol.mol_name;
  this->mol_type = mol.mol_type;
  this->mol_comment = mol.mol_comment;
  this->charge_type = mol.charge_type;
  this->num_atoms = mol.num_atoms;
  this->num_bonds = mol.num_bonds;
  this->num_subst = mol.num_subst;
  this->num_feat = mol.num_feat;
  this->num_sets = mol.num_sets;
  this->atoms.clear();
  this->atoms = mol.atoms;
  this->bonds.clear();
  this->bonds = mol.bonds;
  this->tanimoto = 0;
  this->atomsXYZ = NULL;
  this->radiusAtoms = NULL;
  this->weightAtoms = NULL;
}

Molecule::~Molecule() {

  if (!this->atoms.empty())
    this->atoms.clear();

  if (!this->bonds.empty())
    this->bonds.clear();

  if (this->atomsXYZ != NULL)
    delete[] this->atomsXYZ;

  if (this->radiusAtoms != NULL)
    delete[] this->radiusAtoms;

  if (this->weightAtoms != NULL)
    delete[] this->weightAtoms;
}

void Molecule::setTanimoto(double *atomsXYZ, unsigned int size, double *weight,
                           double *radius, bool sameVDWR) {
  this->tanimoto = VolumeOverlap::preciseOverlapWEGA(
      atomsXYZ, size, weight, radius, atomsXYZ, size, weight, radius, sameVDWR);
}

string Molecule::toString() {

  std::stringstream ss;
  ss << this->mol_name << " : " << this->mol_type
     << "\nnumberOfAtoms: " << this->num_atoms
     << "\nnumberOfBonds: " << this->num_bonds << "\n";
  ss << "Atoms:\n";
  for (unsigned i = 0; i < this->atoms.size(); i++) {
    ss << this->atoms[i].toString() << "\n";
  }
  for (unsigned i = 0; i < this->bonds.size(); i++) {
    ss << this->bonds[i].toString() << "\n";
  }
  std::string s = ss.str();
  return s;
}

void Molecule::ObjectToArray(double *atoms, double *radius, double *weights,
                             bool sameVanderWaalsRadius) {

  int iterador = 0;
  int numberOfAtoms = this->atoms.size();
  for (int i = 0; i < numberOfAtoms; i++) {
    iterador = i * 3;
    atoms[iterador] = this->atoms[i].x;
    atoms[iterador + 1] = this->atoms[i].y;
    atoms[iterador + 2] = this->atoms[i].z;
    radius[i] = this->atoms[i].radiusVanDerWaals;
  }
  iterador = 0;
  double val = 0;
  for (int i = 0; i < numberOfAtoms; i++) {
    iterador = i * 3;
    weights[i] = VolumeOverlap::calculateWeightWEGA(
        atoms[iterador], atoms[iterador + 1], atoms[iterador + 2],
        numberOfAtoms, atoms, radius, i, sameVanderWaalsRadius);
    // cout << "Peso de "<<i<<"["<<atoms[iterador]<<", "<<atoms[iterador +
    // 1]<<", "<<atoms[iterador + 2]<<"]: "<<weights[i]<<"\n";
  }
  // cout << val << endl;
}

void Molecule::calculateWeightWEGAAtoms(double *atoms, unsigned int nAtoms,
                                        double *weights, double *radius,
                                        bool sameVanderWaalsRadius) {
  int iterador = 0;
  for (int i = 0; i < nAtoms; i++) {
    iterador = i * 3;
    weights[i] = VolumeOverlap::calculateWeightWEGA(
        atoms[iterador], atoms[iterador + 1], atoms[iterador + 2], nAtoms,
        atoms, radius, i, sameVanderWaalsRadius);
  }
}

bool Molecule::updateAtomCoordinates(double *coordinates) {
  for (int i = 0; i < this->num_atoms; i++) {
    this->atoms[i].x = coordinates[3 * i];
    this->atoms[i].y = coordinates[3 * i + 1];
    this->atoms[i].z = coordinates[3 * i + 2];
  }

  return true;
}
