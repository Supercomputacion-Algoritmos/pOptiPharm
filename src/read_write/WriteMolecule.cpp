/*
 * WriteMolecule.cpp
 *
 *  Created on: 10 de abr. de 2016
 *      Author: savins
 */
#include "read_write/WriteMolecule.h"
#include "uego/uego.h"

#include <chrono>
#include <fstream>
#include <iostream>
void WriteMolecule::writeMol(string filePath, Molecule *molecule) {
  // auto start = std::chrono::high_resolution_clock::now();
  ofstream fs(filePath.c_str());
  fs << "@<TRIPOS>MOLECULE\n"
     << molecule->mol_name << "\n"
     << molecule->num_atoms << " " << molecule->num_bonds << "\n"
     << molecule->mol_type << "\n"
     << molecule->charge_type << "\n\n";

  fs << "@<TRIPOS>ATOM\n";
  for (long i = 0; i < molecule->num_atoms; ++i) {
    fs << "\t" << molecule->atoms.at(i).atom_id << "\t\t"
       << molecule->atoms.at(i).atom_name << "\t\t" << molecule->atoms.at(i).x
       << "\t\t" << molecule->atoms.at(i).y << "\t\t" << molecule->atoms.at(i).z
       << "\t\t" << molecule->atoms.at(i).atom_type << "\t\t"
       << molecule->atoms.at(i).subst_id << "\t\t"
       << molecule->atoms.at(i).subst_name << "\t\t"
       << molecule->atoms.at(i).charge << "\t\t"
       << molecule->atoms.at(i).status_bit << "\n";
  }
  fs << "@<TRIPOS>BOND\n";
  for (long i = 0; i < molecule->num_bonds; ++i) {
    fs << "\t" << molecule->bonds.at(i).bond_id << "\t\t"
       << molecule->bonds.at(i).origin_atom_id << "\t\t"
       << molecule->bonds.at(i).target_atom_id << "\t\t"
       << molecule->bonds.at(i).bond_type << "\t\t"
       << molecule->bonds.at(i).status_bits << "\n";
  }
  fs.close();
  // message("Archivo molecula guardada!!", MSG_INFORMATION);
  // auto end = std::chrono::high_resolution_clock::now();
  // std::cout <<
  // std::chrono::duration_cast<std::chrono::nanoseconds>(end-start).count() <<
  // "ns" << std::endl;
}

string WriteMolecule::MolToString(Molecule *molecule) {
  // auto start = std::chrono::high_resolution_clock::now();
  string archivo = "";
  archivo = "@<TRIPOS>MOLECULE\n" + molecule->mol_name + "\n" +
            to_string(molecule->num_atoms) + " " +
            to_string(molecule->num_bonds) + "\n" + molecule->mol_type + "\n" +
            molecule->charge_type + "\n\n";

  archivo += "@<TRIPOS>ATOM\n";
  for (long i = 0; i < molecule->num_atoms; ++i) {
    archivo += "\t" + to_string(molecule->atoms.at(i).atom_id) + "\t\t" +
               molecule->atoms.at(i).atom_name + "\t\t" +
               to_string(molecule->atoms.at(i).x) + "\t\t" +
               to_string(molecule->atoms.at(i).y) + "\t\t" +
               to_string(molecule->atoms.at(i).z) + "\t\t" +
               molecule->atoms.at(i).atom_type + "\t\t" +
               to_string(molecule->atoms.at(i).subst_id) + "\t\t" +
               molecule->atoms.at(i).subst_name + "\t\t" +
               to_string(molecule->atoms.at(i).charge) + "\t\t" +
               molecule->atoms.at(i).status_bit + "\n";
  }
  archivo += "@<TRIPOS>BOND\n";
  for (long i = 0; i < molecule->num_bonds; ++i) {
    archivo += "\t" + to_string(molecule->bonds.at(i).bond_id) + "\t\t" +
               to_string(molecule->bonds.at(i).origin_atom_id) + "\t\t" +
               to_string(molecule->bonds.at(i).target_atom_id) + "\t\t" +
               molecule->bonds.at(i).bond_type + "\t\t" +
               molecule->bonds.at(i).status_bits + "\n";
  }
  return archivo;
  // message("Archivo molecula guardada!!", MSG_INFORMATION);
  // auto end = std::chrono::high_resolution_clock::now();
  // std::cout <<
  // std::chrono::duration_cast<std::chrono::nanoseconds>(end-start).count() <<
  // "ns" << std::endl;
}
