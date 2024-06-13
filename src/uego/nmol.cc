#include "functions/MoveAndRotate.h"
#include "functions/Tanimoto.h"
#include "model/Molecule.h"
#include "model/Point3DDouble.h"
#include "read_write/WriteMolecule.h"
#include "uego/searchsp.h"
#include "uego/uego.h"
#include <random>
////////////////////////////////////////////////////////////
// $Id: nreal.cc,v 2.5 1998/03/17 23:14:54 jelasity Exp $
// nreal.cc
// definitions for class NDimRealElement except the Value()
// function; it is defined in a separate file for the
// sake of clarity
////////////////////////////////////////////////////////////
// modification history:
//	Jelasity 98 02 15 inner repr. is (0,1)^dim
////////////////////////////////////////////////////////////

long NDimMolElement::dim = 10;

// -------------------------------------------------------------------------

NDimMolElement::NDimMolElement(long dimension) {

  dim = dimension;

  FailFlag = ((x = new double[dim * 2]) == NULL);
  if (FailFlag)
    message("No memory for NDimMolElement.", MSG_ERROR);
};

// -------------------------------------------------------------------------

double NDimMolElement::UpdateValue() {

  double *normx = x;

  // --- converting to real coordinates
  x += dim;
  for (long i = 0; i < dim; ++i) {
    x[i] = normx[i] * (INI.Upb(i) - INI.Lowb(i)) + INI.Lowb(i);
  };

  // --- normalizacion punto P2 en esfera
  double Vx, Vy, Vz, mV;
  Vx = x[4] - x[1];
  Vy = x[5] - x[2];
  Vz = x[6] - x[3];
  mV = sqrt(Vx * Vx + Vy * Vy + Vz * Vz);

  if (mV != 0) {
    double newP2x = Vx / mV + x[1];
    double newP2y = Vy / mV + x[2];
    double newP2z = Vz / mV + x[3];

    if ((newP2x <= INI.Upb(4) && newP2x >= INI.Lowb(4)) &&
        (newP2y <= INI.Upb(5) && newP2y >= INI.Lowb(5)) &&
        (newP2z <= INI.Upb(6) && newP2z >= INI.Lowb(6))) {
      x[4] = newP2x;
      x[5] = newP2y;
      x[6] = newP2z;
    } else {
      goto valueNoModifications;
    }

    // --- counting value
    value = Value();

    // --- restoring normalized coordinates
    double *xtemp = x;
    x -= dim;

    for (long i = 0; i < dim; i++) {

      if (INI.Upb(i) - INI.Lowb(i) == 0) {
        x[i] = INI.Upb(i);
      } else {
        x[i] = (xtemp[i] - INI.Lowb(i)) / (INI.Upb(i) - INI.Lowb(i));
      }
    };

  } else { // --- counting value

  valueNoModifications:

    value = Value();
    x = normx; // -->Este valia cuando no se ponia el punto en la superficie de
               // la esfera
  }

  return value;
};

double NDimMolElement::PreciseUpdateValue() {

  double *normx = x;

  // --- converting to real coordinates
  x += dim;
  for (long i = 0; i < dim; ++i) {
    x[i] = normx[i] * (INI.Upb(i) - INI.Lowb(i)) + INI.Lowb(i);
  };

  // --- normalizacion punto P2 en esfera
  double Vx, Vy, Vz, mV;
  Vx = x[4] - x[1];
  Vy = x[5] - x[2];
  Vz = x[6] - x[3];
  mV = sqrt(Vx * Vx + Vy * Vy + Vz * Vz);

  if (mV != 0) {
    double newP2x = Vx / mV + x[1];
    double newP2y = Vy / mV + x[2];
    double newP2z = Vz / mV + x[3];

    if ((newP2x <= INI.Upb(4) && newP2x >= INI.Lowb(4)) &&
        (newP2y <= INI.Upb(5) && newP2y >= INI.Lowb(5)) &&
        (newP2z <= INI.Upb(6) && newP2z >= INI.Lowb(6))) {
      x[4] = newP2x;
      x[5] = newP2y;
      x[6] = newP2z;
    } else {
      goto valueNoModifications;
    }

    // --- counting value
    value = PreciseValue();

    // --- restoring normalized coordinates
    double *xtemp = x;
    x -= dim;

    for (long i = 0; i < dim; i++) {

      if (INI.Upb(i) - INI.Lowb(i) == 0) {
        x[i] = INI.Upb(i);
      } else {
        x[i] = (xtemp[i] - INI.Lowb(i)) / (INI.Upb(i) - INI.Lowb(i));
      }
    };

  } else { // --- counting value

  valueNoModifications:

    value = PreciseValue();
    x = normx; // -->Este valia cuando no se ponia el punto en la superficie de
               // la esfera
  }

  return value;
};
// -------------------------------------------------------------------------

double NDimMolElement::Diameter(Ini *ini) {
  // Euclidian distance of lower and upper bounds

  double sqrsum = 0.0;

  for (long i = 0; i < dim; ++i) {
    sqrsum += (ini->Lowb(i) - ini->Upb(i)) * (ini->Lowb(i) - ini->Upb(i));
  };

  return sqrt(sqrsum);
};

// -------------------------------------------------------------------------

double NDimMolElement::v(double r) {
#define SPEED_CONSTANTS 14
  //  Using that binom(x,(x-1)/2)==.5*binom(x+1,(x+1)/2) and
  //  binom(2n,n) is approx. 2**(2n)/sqrt(pi*n)

  double V[SPEED_CONSTANTS] = {// for r==1.0
                               0.0,
                               .25,
                               M_2_PI / 3.0,
                               3.0 / 16.0,
                               4 * M_2_PI / 15.0,
                               10.0 / 64.0,
                               8 * M_2_PI / 35.0,
                               35.0 / 256.0,
                               64 * M_2_PI / 315.0,
                               63.0 / 512.0,
                               128 * M_2_PI / 693.0,
                               231.0 / 2048.0,
                               512 * M_2_PI / 3003.0,
                               858.0 / 8192.0};

  if (dim < SPEED_CONSTANTS)
    return V[dim] * r;
  else
    return r / sqrt(M_PI * (2 * (dim + 1)));

#undef SPEED_CONSTANTS
};

// -------------------------------------------------------------------------

double NDimMolElement::Distance(SearchSpElement *e1, SearchSpElement *e2) {
  // Euclidian distance in real coordinate space

  double sqrsum, d1, d2;
  NDimMolElement *r1, *r2;

  r1 = (NDimMolElement *)e1;
  r2 = (NDimMolElement *)e2;
  if (r2 == NULL)
    r2 = this; // called with 1 par.

  sqrsum = 0.0;
  for (long i = 0; i < dim; ++i) {
    d1 = r1->x[i] * (INI.Upb(i) - INI.Lowb(i)) + INI.Lowb(i);
    d2 = r2->x[i] * (INI.Upb(i) - INI.Lowb(i)) + INI.Lowb(i);
    sqrsum += (d1 - d2) * (d1 - d2);
  };

  return sqrt(sqrsum);
};

// -------------------------------------------------------------------------

SearchSpElement *NDimMolElement::RandNew() {
  // uniform distribution in search space

  NDimMolElement *result;

  FailFlag = 1 == 1;

  // --- new empty element
  result = new NDimMolElement(dim);
  if (result == NULL || result->Fail()) {
    message("No memory in NDimMolElement::RandNew(1)", MSG_ERROR);
    return (SearchSpElement *)NULL;
  };

  // --- filling the empty element
  for (long i = 0; i < dim; ++i)
    result->x[i] = UegoRand();
  result->UpdateValue();

  FailFlag = 1 == 0;
  return result;
};

// -------------------------------------------------------------------------
/** Este es el mismo que el RandNew() pero que se le ha anadido por parametro
 * una especie ya que no interesa que todos los parametros se generen
 * aleatoriamente.
 */
SearchSpElement *NDimMolElement::generateMutationNew(short radind,
                                                     SearchSpElement *e) {
  // uniform distribution in search space

  NDimMolElement *result, *r;

  FailFlag = 1 == 1;

  r = (NDimMolElement *)e;
  if (r == NULL)
    r = this; // called without par.
              // --- new empty element
  result = new NDimMolElement(dim);
  if (result == NULL || result->Fail()) {
    message("No memory in NDimMolElement::RandNew(1)", MSG_ERROR);
    return (SearchSpElement *)NULL;
  };

  // --- filling the empty element

  for (long i = 0; i < dim; ++i) {
    result->x[i] = r->x[i];
  }
  double rad = INI.R(radind) / INI.R(0) * sqrt(dim);

  // Angulo

  result->x[0] = UegoDoubleRand(0, 1);
  if (result->x[0] < 0.0)
    result->x[0] = 0.0;
  else if (result->x[0] > 1.0)
    result->x[0] = 1.0;

  // Traslacion
  for (int i = 7; i < dim; i++) {
    result->x[i] = UegoDoubleRand(0, 1);
    if (result->x[i] < 0.0)
      result->x[i] = 0.0;
    else if (result->x[i] > 1.0)
      result->x[i] = 1.0;
  }
  result->UpdateValue();

  FailFlag = 1 == 0;
  return result;
};

// -------------------------------------------------------------------------

SearchSpElement *NDimMolElement::RandNew(short radind, SearchSpElement *e) {
  // uniform distribution in given area (a hyper-rectangle)

  NDimMolElement *result, *r;
  double rad;

  if (radind == 0)
    return RandNew(); // default for root

  rad = INI.R(radind) / INI.R(0) * sqrt(dim); // normalized radius

  FailFlag = 1 == 1;

  r = (NDimMolElement *)e;
  if (r == NULL)
    r = this; // called without par.

  // --- new empty element
  result = new NDimMolElement(dim);
  if (result == NULL || result->Fail()) {
    message("No memory in NDimMolElement::RandNew(2)", MSG_ERROR);
    return (SearchSpElement *)NULL;
  };

  // --- filling the empty element
  result->x[0] = UegoDoubleRand(r->x[0] - rad, r->x[0] + rad);
  if (result->x[0] < 0.0 || result->x[0] > 1.0)
    result->x[0] = fmod(((fmod(result->x[0], 1.0)) + 1.0), 1.0);

  for (long i = 1; i < dim; ++i) {
    result->x[i] = UegoDoubleRand(r->x[i] - rad, r->x[i] + rad);
    if (result->x[i] < 0.0)
      result->x[i] = 0.0;
    else if (result->x[i] > 1.0)
      result->x[i] = 1.0;
  };
  result->UpdateValue();

  FailFlag = 1 == 0;
  return result;
};

// -------------------------------------------------------------------------

void NDimMolElement::Add(double *y, signed char sign) {

  /*	printf("Before  ");
   for( long i=0; i < dim; ++i )
   printf(" %lf ", x[i]);
   printf(" \n Sumando \n");*/
  for (long i = 0; i < dim; ++i) {
    x[i] += sign * y[i];
    //	printf(" %lf ", x[i]);
    if (x[i] < 0.0)
      x[i] = 0.0;
    else if (x[i] > 1.0)
      x[i] = 1.0;
  };

  /*printf(" \n After   ");
   for( long i=0; i < dim; ++i )
   printf(" %lf ", x[i]);
   printf(" \n");*/
}

// -------------------------------------------------------------------------

double NDimMolElement::Gauss(double bias, double sigma) {

  double xx, z = 0.0;

  for (int k = 0; k < 12; k++) {
    xx = UegoRand();
    z = z + xx;
  }

  return (bias + sigma * (z - 6.0));
};

// -------------------------------------------------------------------------

SearchSpElement *NDimMolElement::MutateNew(short radind, SearchSpElement *e) {
  // mutation used by several algorithms (SHC, GAS, etc.)

  return RandNew(radind, e);
  // ok, just to make it work quickly, temporarily
};

// -------------------------------------------------------------------------

SearchSpElement *NDimMolElement::BetweenNew(SearchSpElement *e1,
                                            SearchSpElement *e2) {

  NDimMolElement *result, *r1, *r2;

  FailFlag = 1 == 1;

  r1 = (NDimMolElement *)e1;
  r2 = (NDimMolElement *)e2;
  if (r2 == NULL)
    r2 = this; // called with 1 par.

  // --- new empty element
  result = new NDimMolElement(dim);
  if (result == NULL || result->Fail()) {
    message("No memory in NDimMolElement::BetweenNew", MSG_ERROR);
    return (SearchSpElement *)NULL;
  };

  // --- filling the empty element
  for (long i = 0; i < dim; ++i)
    result->x[i] = (r1->x[i] + r2->x[i]) / 2.0;
  result->UpdateValue();

  FailFlag = 1 == 0;
  return result;
};

SearchSpElement *NDimMolElement::MutateBetweenNew(SearchSpElement *e1,
                                                  SearchSpElement *e2) {

  NDimMolElement *result, *r1, *r2;

  FailFlag = 1 == 1;

  r1 = (NDimMolElement *)e1;
  r2 = (NDimMolElement *)e2;
  if (r2 == NULL)
    r2 = this; // called with 1 par.

  // --- new empty element
  result = new NDimMolElement(dim);
  if (result == NULL || result->Fail()) {
    message("No memory in NDimMolElement::BetweenNew", MSG_ERROR);
    return (SearchSpElement *)NULL;
  };

  // --- filling the empty element

  for (long i = 0; i < dim; ++i)
    result->x[i] = r1->x[i];

  result->x[0] = (r1->x[0] + r2->x[0]) / 2.0;
  result->UpdateValue();

  FailFlag = 1 == 0;
  return result;
};

// -------------------------------------------------------------------------

// -------------------------------------------------------------------------

void NDimMolElement::UpdateFrom(SearchSpElement *e) {

  value = e->CurrValue();
  for (long i = 0; i < dim; ++i) {
    x[i] = ((NDimMolElement *)e)->x[i];
  };
};

// -------------------------------------------------------------------------

void NDimMolElement::Save(FILE *stream) {

  double y;

  FailFlag = 1 == 0;

  // --- saving as real coordinates
  for (long i = 0; i < dim; ++i) {
    y = x[i] * (INI.Upb(i) - INI.Lowb(i)) + INI.Lowb(i);
    FailFlag |= fprintf(stream, "%.10le ", y) == EOF;
  };

  FailFlag |= fprintf(stream, "%.10le\n", value) == EOF;
  if (Fail())
    message("Error writing NDimMolElement.", MSG_ERROR);
};
/* SAVINS Arregla esto cuando puedas */
string GetStdoutFromCommand3(string cmd) {

  string data;
  FILE *stream;
  const int max_buffer = 256;
  char buffer[max_buffer];
  cmd.append(" 2>&1");

  stream = popen(cmd.c_str(), "r");
  if (stream) {
    while (!feof(stream))
      if (fgets(buffer, max_buffer, stream) != NULL)
        data.append(buffer);
    pclose(stream);
  }
  return data;
}

void NDimMolElement::Save(FILE *stream, bool best) {

  int iterador = 0;
  // UPDATING COORDINATES
  for (unsigned i = 0; i < INI.getMolQuery()->num_atoms; i++) {
    iterador = i * 3;
    INI.getMolQuery()->atoms[i].x = INI.getMolQuery()->getAtomsXYZ()[iterador];
    INI.getMolQuery()->atoms[i].y =
        INI.getMolQuery()->getAtomsXYZ()[iterador + 1];
    INI.getMolQuery()->atoms[i].z =
        INI.getMolQuery()->getAtomsXYZ()[iterador + 2];
  }
  for (unsigned i = 0; i < INI.getMolVariable()->num_atoms; i++) {
    iterador = i * 3;
    INI.getMolVariable()->atoms[i].x =
        INI.getMolVariable()->getAtomsXYZ()[iterador];
    INI.getMolVariable()->atoms[i].y =
        INI.getMolVariable()->getAtomsXYZ()[iterador + 1];
    INI.getMolVariable()->atoms[i].z =
        INI.getMolVariable()->getAtomsXYZ()[iterador + 2];
  }
  // SAVE BEST RESULT
  double y;
  double xarray[dim];
  FailFlag = 1 == 0;

  // --- saving as real coordinates
  for (long i = 0; i < dim; ++i) {
    y = x[i] * (INI.Upb(i) - INI.Lowb(i)) + INI.Lowb(i);
    xarray[i] = y;
    FailFlag |= fprintf(stream, "%.10le ", y) == EOF;
  };
  // We include tanimoto coefficient before to solution value
  // TRANSFORM THE MOLECULE AND SAVE IN MOL2 FILES FOR SHOWING IT
  Molecule *moleculeRotated = new Molecule();
  Point3DDouble p1, p2;
  p1.x = xarray[1];
  p1.y = xarray[2];
  p1.z = xarray[3];
  p2.x = xarray[4];
  p2.y = xarray[5];
  p2.z = xarray[6];

  MoveAndRotate::RotateMolAccording1Axis(INI.getMolVariable(), xarray[0], p1,
                                         p2, moleculeRotated);
  MoveAndRotate::MolToNewPosition(moleculeRotated, xarray[7], xarray[8],
                                  xarray[9]);
  string path = INI.getPathOutputFiles();

  if (INI.getSaveSolutionMol2Files()) {
    string pathFixed = "", pathVariable = "", pathRotated = "";
    pathFixed = path;
    pathFixed.append("MolQ-" + INI.getMolQuery()->mol_name + "-" +
                     INI.getMolVariable()->mol_name + ".mol2");

    WriteMolecule::writeMol(pathFixed, INI.getMolQuery());

    pathVariable = path;
    pathVariable.append("MolD-" + INI.getMolQuery()->mol_name + "-" +
                        INI.getMolVariable()->mol_name + ".mol2");

    WriteMolecule::writeMol(pathVariable, INI.getMolVariable());

    pathRotated = path;
    pathRotated.append("MolR-" + INI.getMolQuery()->mol_name + "-" +
                       INI.getMolVariable()->mol_name + ".mol2");
    WriteMolecule::writeMol(pathRotated, moleculeRotated);
  }
  /*
   * Here it is calculated the initial tanimoto value. We use the initial
   * position of query molecule, variable molecule and the overlap between them.
   */
  double initialTanimoto = 0;
  if (INI.getExecutableName() == "OPShapeSimilarity" |
      INI.getExecutableName() == "pOPShapeSimilarity") {
    initialTanimoto = Tanimoto::calculateTanimotoGeneric(
        INI.getMolQuery()->tanimoto, INI.getMolVariable()->tanimoto,
        VolumeOverlap::overlapWEGA(INI.getMolQuery()->getAtomsXYZ(),
                                   INI.getMolQuery()->atoms.size(),
                                   INI.getMolQuery()->getWeightAtoms(),
                                   INI.getMolQuery()->getRadiusAtoms(),
                                   INI.getMolVariable()->getAtomsXYZ(),
                                   INI.getMolVariable()->atoms.size(),
                                   INI.getMolVariable()->getWeightAtoms(),
                                   INI.getMolVariable()->getRadiusAtoms(),
                                   INI.getSameVanDerWaalsRadius()));
  } else if (INI.getExecutableName() == "OPElectrostatic") {
    string command =
        "./calc_et " + INI.getPathMolQuery() + " " + INI.getPathMolVariable();

    string ls = GetStdoutFromCommand3(command);

    // string ls = GetStdoutFromCommand("python calc_et
    // "+INI.getPathFixedMolecule()+" "+INI.getPathTempMolecule());

    initialTanimoto = atof(ls.c_str());
  }
  FailFlag |= fprintf(stream, "%.10le ", initialTanimoto) == EOF;

  // Normalize best solution with Tanimoto coefficient
  FailFlag |= fprintf(stream, "%.10le\n", value) == EOF; // Original
  // double tanimoto = Tanimoto::calculateTanimotoGeneric(
  // INI.getMoleculeFixed()->tanimoto,
  // INI.getMoleculeVariable()->tanimoto, value);

  // FailFlag |= fprintf(stream, "%.10le\n", tanimoto) == EOF;
  if (Fail())
    message("Error writing NDimMolElement.", MSG_ERROR);
  delete moleculeRotated;
};
