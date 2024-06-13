/*
 * Tanimoto.cpp
 *
 *  Created on: 26 sept. 2016
 *      Author: Savins
 */

#include "functions/Tanimoto.h"
#include <iostream>
using namespace std;
double Tanimoto::calculateTanimotoGeneric(double a, double b, double c) {
  return (c / (a + b - c));
}
