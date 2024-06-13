/*
 * Parameters.cpp
 *
 *  Created on: 16 jan. 2016
 *      Author: Savins
 */

#include "functions/Parameters.h"

using namespace std;

//----- Read all the argument -----------------------------------------
int Parameters::readArguments(int argc, char **argv, int *indexArray,
                              int lenghtIndexArray, Configuration *c) {

  char msg[250];
  if (argc < 2) {

    sprintf(msg, "Version: %s\n", UEGO_VERSION);
    message(msg, MSG_NOTHING);
    message(DEFAULT_MESSAGE, MSG_NOTHING);
    return -1;
  }

  // Leemos los parametros de entrada y obtenemos las configuraciones
  // correspondientes.

  for (int i = 0; i < lenghtIndexArray; i++) {
    indexArray[i] = -1;
  }

  // Vamos obteniendo el valor de cada variable y lo mostramos
  for (int i = 0; i < argc - 1; i++) {
    string arg(argv[i]);
    string argNext(argv[i + 1]);
    if (arg == "-c") { // configuration file of UEGO
      indexArray[0] = i + 1;
      // cout << "-c: " << argNext << endl;
    } else if (arg == "-q") { // query file
      c->onlyConfig = false;
      indexArray[1] = i + 1;
      // cout << "-q: " << argNext << endl;
    } else if (arg == "-d") { // database molecule with compare the query file
      c->onlyConfig = false;
      indexArray[2] = i + 1;
      // cout << "-d: " << argNext << endl;
    } else if (arg == "-o") { // output folder
      indexArray[3] = i + 1;
      // cout << "-o: " << argNext << endl;
    } else if (arg == "-h") { // consider hidrogens atoms
      indexArray[4] = i + 1;
      if (argNext == "1") {
        c->considerHidrogens = true;
      }
      // cout << "-h: " << argNext << endl;
    } else if (arg == "-w") { // consider same Van der Waals radius
      indexArray[5] = i + 1;
      if (argNext == "0") {
        c->sameRadiusVdW = false;
      }
      // cout << "-w: " << argNext << endl;
    } else if (arg == "-N") { // Max functions evals
      c->updateConfig = true;
      indexArray[6] = i + 1;
      // cout << "-N: " << argNext << endl;
    } else if (arg == "-M") { // Max number of species
      c->updateConfig = true;
      indexArray[7] = i + 1;
      // cout << "-M: " << argNext << endl;
    } else if (arg == "-L") { // Number of levels
      c->updateConfig = true;
      indexArray[8] = i + 1;
      // cout << "-L: " << argNext << endl;
    } else if (arg == "-R") { // Minimum radius last level
      c->updateConfig = true;
      indexArray[9] = i + 1;
      // cout << "-R: " << argNext << endl;
    } else if (arg == "-align") { // Do preprocess PCA script.py
      indexArray[10] = i + 1;
      if (argNext == "0")
        c->align = false;
      // cout << "-pre: " << argNext << endl;
    } else if (arg == "-pz") { // set P1 to (0,0,0)
      indexArray[11] = i + 1;
      if (argNext == "1")
        c->pointOneZero = true;
      // cout << "-pz: " << argNext << endl;
    } else if (arg == "-smol") { // save solution as mol2 file
      indexArray[12] = i + 1;
      if (argNext == "1")
        c->saveMolFile = true;

      // cout << "-smol: " << argNext << endl;
    } else if (arg == "-fvdwr") { // save solution as mol2 file
      indexArray[13] = i + 1;
      Atom::setVDWR(argNext);
      // cout << "-smol: " << argNext << endl;
    } else if (arg == "-lowL") { // LowLimit Score
      indexArray[14] = i + 1;

      // cout << "-smol: " << argNext << endl;
    } else if (arg == "-th") {
      indexArray[15] = i + 1;
      uint32_t th_num = atoi(argNext.c_str());
      if (th_num > 0) {
        c->threads_num = th_num;
      }
    } else if (arg == "-pn") {
      indexArray[16] = i + 1;
      if (argNext == "1") {
        c->pin_threads = true;
      }
    } else if (arg == "-fgpu") {
      indexArray[17] = i + 1;
      if (argNext == "1") {
        c->force_CUDA = true;
      }
    } else if (arg == "-mc") {
      indexArray[18] = i + 1;
      uint32_t mc_limit = atoi(argNext.c_str());
      if (mc_limit > 0) {
        c->min_combined_size = mc_limit;
      }
    }
  }

  return 0;
}
//------ END -----------------------------------------------------------
//-----------------------------------------------------------------------
string Parameters::getNameMol(char *path) {
  const char s[2] = "/";
  char *token;

  /* get the first token */
  token = strtok(path, s);

  /* walk through other tokens */
  string name;

  while (token != NULL) {
    name = token;
    token = strtok(NULL, s);
  }

  token = new char[name.length() + 1];
  strcpy(token, name.c_str());
  name = strtok(token, ".");
  return name;
}

void Parameters::defineDeltaMovement(Ini *ini, double *minX, double *maxX,
                                     double *minY, double *maxY, double *minZ,
                                     double *maxZ) {

  /**
   * Asi lo definimos para que quede la caja mas pequena.
   */
  Point3DDouble maxPointMolFixed;
  Point3DDouble minPointMolFixed;
  VolumeOverlap::boundingBox(*(ini->getMolQuery()), &minPointMolFixed,
                             &maxPointMolFixed);
  // cout << "Limites compuesto query\n";
  // cout << minPointMolFixed.x << ", "<< minPointMolFixed.y << ", "<<
  // minPointMolFixed.z << "\n"; cout << maxPointMolFixed.x << ", "<<
  // maxPointMolFixed.y << ", "<< maxPointMolFixed.z << "\n";

  Point3DDouble maxPointMolVariable;
  Point3DDouble minPointMolVariable;

  VolumeOverlap::boundingBox(*(ini->getMolVariable()), &minPointMolVariable,
                             &maxPointMolVariable);

  // cout << "Limites compuesto target\n";
  // cout << minPointMolVariable.x << ", "<< minPointMolVariable.y << ", "<<
  // minPointMolVariable.z << "\n"; cout << maxPointMolVariable.x << ", "<<
  // maxPointMolVariable.y << ", "<< maxPointMolVariable.z << "\n";

  double deltaXMin = 0, deltaXMax = 0, deltaYMin = 0, deltaYMax = 0,
         deltaZMin = 0, deltaZMax = 0;
  // cout << minPointMolFixed.toString()<< "  "
  // <<maxPointMolFixed.toString()<<'\n'; cout <<
  // minPointMolVariable.toString()<<" "<<maxPointMolVariable.toString()<<'\n';

  deltaXMin = minPointMolVariable.x - minPointMolFixed.x;
  deltaYMin = minPointMolVariable.y - minPointMolFixed.y;
  deltaZMin = minPointMolVariable.z - minPointMolFixed.z;

  deltaXMax = maxPointMolVariable.x - maxPointMolFixed.x;
  deltaYMax = maxPointMolVariable.y - maxPointMolFixed.y;
  deltaZMax = maxPointMolVariable.z - maxPointMolFixed.z;

  *minX = deltaXMin;
  *minY = deltaYMin;
  *minZ = deltaZMin;
  *maxX = deltaXMax;
  *maxY = deltaYMax;
  *maxZ = deltaZMax;
  // cout << "diferencias cajas\n";
  // cout << *minX << ", " << *minY << ", " << *minZ << ", " << *maxX << ", " <<
  // *maxY << ", " << *maxZ <<'\n';
  double X = 0, Y = 0, Z = 0;
  (abs(deltaXMin) < abs(deltaXMax)) ? X = abs(deltaXMax) : X = abs(deltaXMin);
  (abs(deltaYMin) < abs(deltaYMax)) ? Y = abs(deltaYMax) : Y = abs(deltaYMin);
  (abs(deltaZMin) < abs(deltaZMax)) ? Z = abs(deltaZMax) : Z = abs(deltaZMin);
  if (X == 0) {
    X = 0.00001;
  }
  if (Y == 0) {
    Y = 0.00001;
  }
  if (Z == 0) {
    Z = 0.00001;
  }
  *minX = -X;
  *minY = -Y;
  *minZ = -Z;
  *maxX = X;
  *maxY = Y;
  *maxZ = Z;
  // cout << "parametros finales\n";
  // cout << *minX << ", " << *minY << ", " << *minZ << ", " << *maxX << ", " <<
  // *maxY << ", " << *maxZ <<'\n';
}

void Parameters::alignMolecule(Molecule *mol) {
  if (mol->num_atoms < 4)
    return;
  // This process centers the molecule in the center of coordinates.
  double centroidX = 0, centroidY = 0, centroidZ = 0;
  int iteradori = 0;
  unsigned numberOfAtoms = mol->atoms.size();

  for (unsigned i = 0; i < numberOfAtoms; i++) {

    centroidX += mol->atoms[i].x;
    centroidY += mol->atoms[i].y;
    centroidZ += mol->atoms[i].z;
  }
  centroidX /= numberOfAtoms;
  centroidY /= numberOfAtoms;
  centroidZ /= numberOfAtoms;
  for (unsigned i = 0; i < numberOfAtoms; i++) {
    mol->atoms[i].x -= centroidX;
    mol->atoms[i].y -= centroidY;
    mol->atoms[i].z -= centroidZ;
  }

  int totalValue = numberOfAtoms * 3;
  double matrix[totalValue];
  unsigned contador = 0;
  for (int i = 0; i < numberOfAtoms; i++) {
    matrix[i] = mol->atoms[i].x;
    matrix[i + numberOfAtoms] = mol->atoms[i].y;
    matrix[i + numberOfAtoms + numberOfAtoms] = mol->atoms[i].z;
  }

  // This process align the molecule according to the three axis.
  // We convert the array in a eigen MAtrix. We need to give data, number of
  // rows and number of columns
  MatrixXd eigenX = Map<MatrixXd>(matrix, numberOfAtoms, 3);

  // It does stuff
  JacobiSVD<MatrixXd> svd(eigenX, ComputeThinU | ComputeThinV);

  MatrixXd VT = svd.matrixV();

  if (numberOfAtoms != 3) {

    if (VT.determinant() < 0) {
      VT = -1 * VT;
    }
  }
  MatrixXd m(1, 3);
  for (unsigned i = 0; i < numberOfAtoms; i++) {

    m(0, 0) = mol->atoms[i].x;
    m(0, 1) = mol->atoms[i].y;
    m(0, 2) = mol->atoms[i].z;

    m = m * VT;

    mol->atoms[i].x = m(0, 0);
    mol->atoms[i].y = m(0, 1);
    mol->atoms[i].z = m(0, 2);
    //	cout << i << ": " << mol->atoms[i].x << ", " << mol->atoms[i].y << ", "
    //<< mol->atoms[i].z << "\n";
  }
}
