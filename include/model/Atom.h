/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Atom.h
 * Author: Savins
 *
 * Created on 25 de noviembre de 2015, 16:26
 */

#ifndef ATOM_H
#define ATOM_H

#include <string>
#include <sstream>
#include <unordered_map>
#include <iostream>
#include <fstream>
using namespace std;

class Atom {
public:
    string atom_name, atom_type, subst_name, status_bit;
    int atom_id, subst_id, old_atom_id;
    double x, y, z, charge, radiusVanDerWaals;
    static unordered_map<string, double> defaultVDWR;
    static string periodicTable[];
    Atom();
    Atom(int atom_id, string atom_name, double x, double y, double z, string
    atom_type, int subst_id, string subst_name, double charge, string status_bit);
    virtual ~Atom();
    /**
     * El radio es el dado por los autores de WEGA en un correo el dia 17/11/2016 Re: Fwd:letras chinas: New WEGA license and version
     */
    static double setRadiusVanDerWaals(string atom_type);
    string toString();
    /**
     * Lee de un archivo los radios de van der waals y modifica los radios por defecto.
     * Los valores de los radios deben de indicarse en una linea separados por comas:
     * Ej.
     * H,He,Na....
     */
    static void setVDWR(string file);
};


#endif /* ATOM_H */

