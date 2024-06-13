#ifndef INI_H
#define INI_H

class SearchSpElement;

#include "model/Molecule.h"
#include <cstdlib>
#include <cstring>
#include <vector>
////////////////////////////////////////////////////////////
// $Id: uegoini.h,v 2.5 1998/03/17 23:14:52 jelasity Exp $
// uegoini.h
// contains the declaration of class Ini
// Ini contains and handles the settings used in an uego
// session.
////////////////////////////////////////////////////////////
// modification history:
//	jelasity 98 02 20 name changed to uegoini.h
////////////////////////////////////////////////////////////

class Ini {
  friend class Configure;

private:
  // --- random numbers --------------------------------------
  unsigned long seed; // random seed

  // --- objective function -------------------------------------
  long type;                  // type of search space
  long fnum;                  // number of obj. funtc. in library
  long dimension;             // dimesion of the search space
  long paramnum;              // length of param[]
  double *param;              // optional parameters fot the obj. funct.
  double *lowb;               // lower bounds of variables for real spaces
  double *upb;                // upper bounds of variables for real spaces
  SearchSpElement *prototype; // points to an instance of the
  // class selected by the user

  // --- uego values -----------------------------------------------
  unsigned long maxevals; // max num. of funct. evals
  long maxspecnum;        // max number of species at a time
  double threshold;       // stability threshold
  short levels;           // max strict level (length of r[])
  double last_r;          // radius of last level
  short level_actual;
  // --- automatic
  double *r;                   // radii of different levels
  unsigned long *evals;        // max funct. evals on  different levels
  unsigned long *newspecevals; // max funct. evals in NewSpecies

  // --- handy things for Save and Ini( FILE* ) ----------------------
  static char SaveVector(FILE *, char *, double *, unsigned long *, long);
  static char ToBuffer(FILE *, char **);
  static char GetVector(char *, char *, double **, unsigned long **, long);
  double GetValue(char *, char *);

  char FailFlag;

  Molecule *molQuery = NULL;
  Molecule *molVariable = NULL;
  Molecule *molTemp = NULL;

  string nameFileMolQuery;
  string nameFileMolVariable;

  string pathFileMolQuery;
  string pathFileMolVariable;

  string pathFileMolTemp;

  string pathOutputFiles;
  bool sameVanDerWaalsRadius;
  bool considerHydrogens;
  bool saveSolutionMol2Files;

  double lowLimit;

  string name_executable;
  string MolQueryToString;

  // Thread stuff
  uint32_t thread_num;
  bool pin_threads;

  // CUDA stuff
  bool force_CUDA;
  uint32_t min_combined_size;

public:
  void setLevelMaster(short level_tem) { level_actual = level_tem; }
  short getLevelMaster() { return level_actual; }

  void setThreadNumber(uint32_t thread_number) {
    this->thread_num = thread_number;
  }
  uint32_t getThreadNumber() { return this->thread_num; }
  void setPinnedThreads(bool pin_th) { this->pin_threads = pin_th; }
  bool isPinnedThreads() { return this->pin_threads; }

  bool isForcingCUDA() { return this->force_CUDA; }
  void setForcingCUDA(bool force) { this->force_CUDA = force; }
  void setCUDAMinCombinedSize(uint32_t min_size) {
    this->min_combined_size = min_size;
  }
  uint32_t getCUDAMinCombinedSize() { return this->min_combined_size; }

  Ini(FILE *);
  Ini(FILE *, int argc, char **argv);
  Ini();
  ~Ini();
  double V(double r);
  unsigned long *Seed() { return &seed; };
  // random seed; modified by UegoRand()
  // that's why it's a pointer

  long Type() { return type; };
  long Fnum() { return fnum; };
  long Dimension() { return dimension; };
  long ParamNum() { return paramnum; };
  double Param(long i) {
    if (i < 0)
      i = 0;
    else if (i >= paramnum)
      i = paramnum - 1;
    return param == NULL ? 0.0 : param[i];
  };
  double Lowb(long i) {
    if (i < 0)
      i = 0;
    else if (i >= dimension)
      i = dimension - 1;
    return lowb == NULL ? 0.0 : lowb[i];
  };
  double Upb(long i) {
    if (i < 0)
      i = 0;
    else if (i >= dimension)
      i = dimension - 1;
    return upb == NULL ? 0.0 : upb[i];
  };

  unsigned long MaxEvals() { return maxevals; };
  long MaxSpecNumber() { return maxspecnum; };
  short Levels() { return levels; };
  double R(long i) {
    if (i < 0)
      i = 0;
    else if (i >= levels)
      i = levels - 1;
    return r == NULL ? 0.0 : r[i];
  };
  unsigned long Evals(long i) {
    if (i < 0)
      i = 0;
    else if (i >= levels)
      i = levels - 1;
    return evals == NULL ? 0 : evals[i];
  };
  unsigned long NewSpecEvals(long i) {

    if (i < 1)
      i = 1;
    else if (i >= levels)
      i = levels - 1;
    return newspecevals == NULL ? 0 : newspecevals[i - 1];
  };

  SearchSpElement *Prototype() { return prototype; };
  void SetPrototype(); // creates prototype of type 'type'

  void Save(FILE *);

  char Fail() { return FailFlag; };

  void setLowb(double *newLowb) {
    for (unsigned i = 0; i < dimension; i++) {
      if (*(newLowb + i) != -19041993) {
        *(lowb + i) = *(newLowb + i);
      }
    }
  }
  void setUpb(double *newUpb) {
    for (unsigned i = 0; i < dimension; i++) {
      if (*(newUpb + i) != -19041993) {
        *(upb + i) = *(newUpb + i);
      }
    }
  }

  Molecule *getMolQuery() { return molQuery; }

  void setPathMolQuery(string path) { this->pathFileMolQuery = path; }

  void setPathMolVariable(string path) { this->pathFileMolVariable = path; }

  string getPathMolQuery() { return this->pathFileMolQuery; }

  string getPathMolVariable() { return this->pathFileMolVariable; }

  Molecule *getMolVariable() { return molVariable; }
  void setMolQuery(Molecule *query) {
    delete this->molQuery;
    this->molQuery = new Molecule(*query);
  }
  void setMolVariable(Molecule *newMoleculev) {
    delete this->molVariable;
    this->molVariable = new Molecule(*newMoleculev);
  }

  string getNameFileMolQuery() { return this->nameFileMolQuery; }

  void setNameFileMolQuery(string name) { this->nameFileMolQuery = name; }

  string getNameFileMolVariable() { return this->nameFileMolVariable; }

  void setNameFileMolVariable(string name) { this->nameFileMolVariable = name; }

  string getPathOutputFiles() { return this->pathOutputFiles; }

  void setPathOutputFiles(string path) { this->pathOutputFiles = path; }

  bool getSameVanDerWaalsRadius() { return this->sameVanDerWaalsRadius; }

  void setSameVanDerWaals(bool same) { this->sameVanDerWaalsRadius = same; }

  bool getConsiderHydrogens() { return this->considerHydrogens; }

  void setConsiderHydrogens(bool consider) {
    this->considerHydrogens = consider;
  }

  bool getSaveSolutionMol2Files() { return this->saveSolutionMol2Files; }
  void setSaveSolutionMol2Files(bool save) {
    this->saveSolutionMol2Files = save;
  }

  string getPathMolTemp() { return this->pathFileMolTemp; }
  void setPathMolTemp(string path) { this->pathFileMolTemp = path; }
  /**
   * Update r, evals, newspeciesevals vector and threshold value. It is
   * necessary because after to obtain the value for each molecule, these values
   * are calculated with automatic process.
   */
  void updateVectorIni();

  void setr(double *newr) { r = newr; }
  void setEvals(unsigned long *newEvals) { this->evals = newEvals; }

  void setNewspecevals(unsigned long *newNewspecevals) {
    this->newspecevals = newNewspecevals;
  }
  void setThreshold(double newThreshold) { this->threshold = newThreshold; }

  double getLastR() { return this->last_r; }
  short getLevels() { return this->levels; }
  void setFailFlag(char newFailFlag) { this->FailFlag = newFailFlag; }

  double getTreshhold() { return this->threshold; }

  void setLastR(double newLastR) { this->last_r = newLastR; }

  void setLowLimit(double newLowLimit) { this->lowLimit = newLowLimit; }
  double getLowLimit() { return this->lowLimit; }

  Molecule *getMoleculeTemp() { return this->molTemp; }

  void setMoleculeTemp(Molecule *molecule) {
    if (this->molTemp != NULL)
      delete this->molTemp;
    this->molTemp = new Molecule(*molecule);
  }
  void setExecutableName(string name) { this->name_executable = name; }
  string getExecutableName() { return name_executable; }
  void setMolQueryToString(string query) { this->MolQueryToString = query; }
  string getMolQueryToString() { return MolQueryToString; }
};

#endif
