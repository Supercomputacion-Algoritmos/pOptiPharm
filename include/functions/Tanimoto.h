/*
 * Tanimoto.h
 *
 *  Created on: 26 sept. 2016
 *      Author: Savins
 */

#ifndef SRC_FUNCTIONS_TANIMOTO_H_
#define SRC_FUNCTIONS_TANIMOTO_H_

class Tanimoto {
public:
	/**
	 * Calcule the Tanimoto similarity.
	 * a: value of compare first element itself
	 * b: value of compare second element itself
	 * c: value of compare first element with second one.
	 *
	 */
	 double static calculateTanimotoGeneric(double a, double b, double c);
};


#endif /* SRC_FUNCTIONS_TANIMOTO_H_ */
