/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   Atom.cpp
 * Author: Savíns
 *
 * Created on 25 de noviembre de 2015, 16:26
 */

#include "model/Atom.h"

string Atom::periodicTable[] = {
    "",   "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na",
    "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V",
    "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br",
    "Kr", "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag",
    "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr",
    "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
    "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi",
    "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am",
    "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh",
    "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"};

unordered_map<string, double> Atom::defaultVDWR = {
    {"", 0.0},    {"H", 1.20},  {"He", 1.40}, {"Li", 1.82}, {"Be", 2.00},
    {"B", 2.00},  {"C", 1.70},  {"N", 1.55},  {"O", 1.52},  {"F", 1.47},
    {"Ne", 1.54}, {"Na", 2.27}, {"Mg", 1.73}, {"Al", 2.00}, {"Si", 2.10},
    {"P", 1.80},  {"S", 1.80},  {"Cl", 1.75}, {"Ar", 1.88}, {"K", 2.75},
    {"Ca", 2.00}, {"Sc", 2.00}, {"Ti", 2.00}, {"V", 2.00},  {"Cr", 2.00},
    {"Mn", 2.00}, {"Fe", 2.00}, {"Co", 2.00}, {"Ni", 1.63}, {"Cu", 1.40},
    {"Zn", 1.39}, {"Ga", 1.87}, {"Ge", 2.00}, {"As", 1.85}, {"Se", 1.90},
    {"Br", 1.85}, {"Kr", 2.02}, {"Rb", 2.00}, {"Sr", 2.00}, {"Y", 2.00},
    {"Zr", 2.00}, {"Nb", 2.00}, {"Mo", 2.00}, {"Tc", 2.00}, {"Ru", 2.00},
    {"Rh", 2.00}, {"Pd", 1.63}, {"Ag", 1.72}, {"Cd", 1.58}, {"In", 1.93},
    {"Sn", 2.17}, {"Sb", 2.00}, {"Te", 2.06}, {"I", 1.98},  {"Xe", 2.16},
    {"Cs", 2.00}, {"Ba", 2.00}, {"La", 2.00}, {"Ce", 2.00}, {"Pr", 2.00},
    {"Nd", 2.00}, {"Pm", 2.00}, {"Sm", 2.00}, {"Eu", 2.00}, {"Gd", 2.00},
    {"Tb", 2.00}, {"Dy", 2.00}, {"Ho", 2.00}, {"Er", 2.00}, {"Tm", 2.00},
    {"Yb", 2.00}, {"Lu", 2.00}, {"Hf", 2.00}, {"Ta", 2.00}, {"W", 2.00},
    {"Re", 2.00}, {"Os", 2.00}, {"Ir", 2.00}, {"Pt", 1.72}, {"Au", 1.66},
    {"Hg", 1.55}, {"Tl", 1.96}, {"Pb", 2.02}, {"Bi", 2.00}, {"Po", 2.00},
    {"At", 2.00}, {"Rn", 2.00}, {"Fr", 2.00}, {"Ra", 2.00}, {"Ac", 2.00},
    {"Th", 2.00}, {"Pa", 2.00}, {"U", 1.86},  {"Np", 2.00}, {"Pu", 2.00},
    {"Am", 2.00}, {"Cm", 2.00}, {"Bk", 2.00}, {"Cf", 2.00}, {"Es", 2.00},
    {"Fm", 2.00}, {"Md", 2.00}, {"No", 2.00}, {"Lr", 2.00}, {"Rf", 2.00},
    {"Db", 2.00}, {"Sg", 2.00}, {"Bh", 2.00}, {"Hs", 2.00}, {"Mt", 2.00},
    {"Ds", 2.00}, {"Rg", 2.00}};

Atom::Atom() {
  this->atom_id = 0;
  this->old_atom_id = 0;
  this->atom_name = "";
  this->x = 0;
  this->y = 0;
  this->z = 0;
  this->atom_type = "";
  this->subst_id = 0;
  this->subst_name = "";
  this->charge = 0;
  this->status_bit = "";
  this->radiusVanDerWaals = 0;
}

Atom::Atom(int atom_id, string atom_name, double x, double y, double z,
           string atom_type, int subst_id, string subst_name, double charge,
           string status_bit) {
  this->atom_id = atom_id;
  this->old_atom_id = atom_id;
  this->atom_name = atom_name;
  this->x = x;
  this->y = y;
  this->z = z;
  this->atom_type = atom_type;
  this->subst_id = subst_id;
  this->subst_name = subst_name;
  this->charge = charge;
  this->status_bit = status_bit;
}

Atom::~Atom() {}

double Atom::setRadiusVanDerWaals(string atom_type) {

  istringstream f(atom_type);
  string s;
  getline(f, s, '.');
  double radio;
  std::unordered_map<std::string, double>::const_iterator got =
      defaultVDWR.find(s);
  if (got == defaultVDWR.end())
    radio = 1.80;
  else
    radio = got->second;

  return radio;
}

string Atom::toString() {
  std::stringstream ss;
  ss << this->atom_id << ":" << this->atom_name << " [" << this->x << ", "
     << this->y << ", " << this->z << "] : " << this->radiusVanDerWaals;
  std::string s = ss.str();
  return s;
}

void Atom::setVDWR(string file) {

  string linea;
  ifstream infile(file.c_str());

  int i = 1;
  std::string delimiter = ",";
  while (std::getline(infile, linea)) {
    size_t pos = 0;
    std::string token;
    while ((pos = linea.find(delimiter)) != std::string::npos) {
      token = linea.substr(0, pos);
      defaultVDWR[Atom::periodicTable[i]] = atof(token.c_str());
      i = i + 1;
      linea.erase(0, pos + delimiter.length());
    }
  }
  /*
   for (auto& x: defaultVDWR) {
   std::cout << x.first << ": " << x.second << std::endl;
   }*/
}
