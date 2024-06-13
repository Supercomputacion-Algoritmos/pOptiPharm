/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   Bond.cpp
 * Author: Savíns
 *
 * Created on 25 de noviembre de 2015, 16:27
 */

#include "model/Bond.h"

Bond::Bond() {
  this->bond_id = 0;
  this->origin_atom_id = 0;
  this->target_atom_id = 0;
  this->bond_type = "";
  this->status_bits = "";
}

Bond::Bond(int bond_id, int origin_atom_id, int target_atom_id,
           string bond_type, string status_bits) {
  this->bond_id = bond_id;
  this->origin_atom_id = origin_atom_id;
  this->target_atom_id = target_atom_id;
  this->bond_type = bond_type;
  this->status_bits = status_bits;
}

Bond::~Bond() {}

string Bond::toString() {
  std::stringstream ss;
  ss << this->bond_id << ":[" << this->origin_atom_id << " -> "
     << this->target_atom_id << "] " << this->bond_type;
  std::string s = ss.str();
  return s;
}
