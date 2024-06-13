/*
 * ManageEspecies.h
 *
 *  Created on: 4 de abr. de 2016
 *      Author: savins
 */

#ifndef SRC_UEGO_MANAGESPECIES_H_
#define SRC_UEGO_MANAGESPECIES_H_
#include "uego.h"
#include <iostream>
using namespace std;

class ManageSpecies {
private:
  SpeciesList *start, *temp1, *temp2, *temp3;

public:
  ManageSpecies();
  virtual ~ManageSpecies();
  void addnode();
  void delnode();
  void display();
  void show();
};

#endif /* SRC_UEGO_MANAGESPECIES_H_ */
