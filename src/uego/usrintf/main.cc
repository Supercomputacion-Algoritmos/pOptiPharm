// Este main lee dos moleculas de dos archivos independientes en formato mol2.
// Despues ejecuta UEGO

#include "time.h"

#include "functions/Parameters.h"
#include "model/Molecule.h"
#include "model/Point3DDouble.h"
#include "read_write/ReadMolecule.h"
#include "read_write/WriteMolecule.h"
#include "uego/configur.h"
#include "uego/uego.h"
#include "uego/usrintf.h"
#include <chrono>
#include <iostream>
#include <thread>
#include <typeinfo>
#include <vector>

////////////////////////////////////////////////////////////
// $Id: main.cc,v 2.5 1998/03/17 23:14:55 jelasity Exp $
// main.cc
// main modul of the user interface
////////////////////////////////////////////////////////////
// modification history:
//	jelasity 98 01 17 command line configure
//		98 02 14 non-documented command line passing feature
////////////////////////////////////////////////////////////

Ini *Master::_ini = NULL; // initialization of static member
void setMsgLevel(char);

//--- The following is a non-documented feature, to use it in a modul -----
//--- declare global variables as 'extern' --------------------------------

char **RemArgv = NULL; // command line pars after '--'
int RemArgc = 0;       // num. of command line pars after '--'

void CutCommandLine(int *argc, char **argv) {

  for (int i = 0; i < *argc; ++i) {
    if (strcmp(argv[i], "--") == 0) {
      RemArgc = *argc - (i + 1);
      if (i != *argc - 1)
        RemArgv = argv + (i + 1);
      *argc = i;
      break;
    };
  };
};

//--- These are local to this modul ------------------------------------

unsigned long saveflags; // GetPars() sets it
long repcount;           // GetPars() sets it
char *trace,             // GetPars() sets it
    msg[250];            // working

//-----------------------------------------------------------------------

char GetPars(int argc, char **argv) {
  // sets repcount and trace; returns false on failure

  trace = NULL;
  repcount = 1;
  saveflags = 0;

  for (int i = 0; i < argc; ++i) {
    if (argv[i][0] == '-') {
      switch (argv[i][1]) {
      case 'T':
        if (trace != NULL)
          delete trace;
        trace = new char[strlen(argv[i] + 2) + 1];
        if (trace == NULL) {
          message("No memory for comm. line.", MSG_ERROR);
          return 1 == 0;
        };
        strcpy(trace, argv[i] + 2);
        break;
      case 'r':
        if (sscanf(argv[i] + 2, "%ld", &repcount) != 1) {
          message("Bad repeate count.", MSG_ERROR);
          return 1 == 0;
        };
        break;
      case 'S':
        saveflags |= SHORT_SAVE;
        break;
      };
    }
  };

  return 1 == 1;
};
/**
 * Get name from path
 */
string getFileName(const string &s) {

  char sep = '/';
  size_t i = s.rfind(sep, s.length());
  string name = s;
  if (i != string::npos) {
    name = s.substr(i + 1, s.length() - i);
  }
  size_t lastindex = name.find_last_of(".");
  return (name.substr(0, lastindex));
}
/**
 * Redirect the output to a string
 */
string GetStdoutFromCommand(string cmd) {

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

int main(int argc, char **argv) {

  int errcode = 0;
  Ini *ini = NULL;
  Master *master = NULL;
  FILE *inifile = NULL;
  double tiemporeloj;
  // time_t t1, t2;
  // clock_t clock1, clock2;

  setMsgLevel(MSG_INFORMATION);

  // No se para que sirve, supongo que seria de otros experimentos
  // CutCommandLine(&argc, argv);	// cuts things after '--'

  // 0: -c, 1: -q, 2: -d, 3: -o, 4: -h, 5: -w, 6: -N, 7: -M, 8: -L, 9: -R, 10:
  // -pre, 11: -pz, 12: -smol, 13: -fvdwr, 14: -lowL
  int *indexArray = NULL;

  int lenghtIndexArray = 19;
  indexArray = new int[lenghtIndexArray];

  Parameters::Configuration c;
  if (Parameters::readArguments(argc, argv, indexArray, lenghtIndexArray, &c) ==
      -1) {
    if (master != NULL)
      delete master;

    if (trace != NULL)
      delete trace;

    if (ini != NULL)
      delete ini;

    if (indexArray != NULL)
      delete[] indexArray;
    return -1;
  }

  // Actualizamos el uegoini
  if (c.updateConfig) { // Si es necesario modificar el archivo de configuracion
    sprintf(msg, "%s.%s", UEGO_ININAME, argv[indexArray[0]]);
    Configure::ConfigUpdating(msg, argc, argv);

    if (Configure::Fail())
      return 7;
  }

  if (c.onlyConfig) {
    message("Updated configuration file.", MSG_INFORMATION);
    return 0;
  }

  // Se abre el archivo de configuracion
  sprintf(msg, "Opening ini file %s.%s", UEGO_ININAME, argv[indexArray[0]]);
  message(msg, MSG_INFORMATION);
  sprintf(msg, "%s.%s", UEGO_ININAME, argv[indexArray[0]]);
  inifile = fopen(msg, "rt");

  // Se leen las moleculas
  Molecule *query = new Molecule();
  ReadMolecule::readMol(argv[indexArray[1]], query, c.considerHidrogens,
                        c.sameRadiusVdW);

  vector<Molecule> db;
  ReadMolecule::readDBMol2(argv[indexArray[2]], &db, c.considerHidrogens,
                           c.sameRadiusVdW);
  for (unsigned int imol = 0; imol < db.size(); imol++) {

    // --- Starting optimalization -------------------------------------
    message("Reading ini file.", MSG_INFORMATION);
    ini = new Ini(inifile);

    // Set INI thread stuff
    ini->setThreadNumber(c.threads_num);
    ini->setPinnedThreads(c.pin_threads);

    ini->setForcingCUDA(c.force_CUDA);
    ini->setCUDAMinCombinedSize(c.min_combined_size);

    ini->setExecutableName(getFileName(argv[0]));
    ini->setConsiderHydrogens(c.considerHidrogens);
    ini->setSameVanDerWaals(c.sameRadiusVdW);
    string path1(argv[indexArray[1]]);
    ini->setPathMolQuery(path1);
    string path2(argv[indexArray[2]]);
    ini->setPathMolVariable(path2);
    // ---- Output Files ----------------------------------------------------
    if (indexArray[3] != -1) {
      string path(argv[indexArray[3]]);
      ini->setPathOutputFiles(path);
    }
    // ----- END -----------------------

    // ---- ITERAMOS EN CADA MOLECULA
    ini->setMolQuery(query);
    if (c.align) {
      Parameters::alignMolecule(ini->getMolQuery());
    }
    ini->setMolQueryToString(WriteMolecule::MolToString(ini->getMolQuery()));

    ini->getMolQuery()->setArraysAtoms(query->atoms.size());
    // cout << "Query"<<endl;
    ini->getMolQuery()->ObjectToArray(
        ini->getMolQuery()->getAtomsXYZ(), ini->getMolQuery()->getRadiusAtoms(),
        ini->getMolQuery()->getWeightAtoms(), ini->getSameVanDerWaalsRadius());

    ini->setMolVariable(&(db[imol]));
    if (c.align) {
      Parameters::alignMolecule(ini->getMolVariable());
    }
    ini->setMoleculeTemp(&(db[imol])); // Molecula temporal para rotar
    // Se alinea el compuesto

    ini->getMolVariable()->setArraysAtoms(db[imol].atoms.size());

    // cout << "Target"<<endl;
    ini->getMolVariable()->ObjectToArray(
        ini->getMolVariable()->getAtomsXYZ(),
        ini->getMolVariable()->getRadiusAtoms(),
        ini->getMolVariable()->getWeightAtoms(),
        ini->getSameVanDerWaalsRadius());

    // ---- Getting Tanimoto Value -----------------------------------------
    // if(ini->getExecutableName()=="OPShapeSimilarity"){

    ini->getMolQuery()->setTanimoto(
        ini->getMolQuery()->getAtomsXYZ(), ini->getMolQuery()->atoms.size(),
        ini->getMolQuery()->getWeightAtoms(),
        ini->getMolQuery()->getRadiusAtoms(), ini->getSameVanDerWaalsRadius());

    ini->getMolVariable()->setTanimoto(ini->getMolVariable()->getAtomsXYZ(),
                                       ini->getMolVariable()->atoms.size(),
                                       ini->getMolVariable()->getWeightAtoms(),
                                       ini->getMolVariable()->getRadiusAtoms(),
                                       ini->getSameVanDerWaalsRadius());
    //}

    // ------ END -----------------------------

    // ---- Se copia la molecula variable  para poder modificarla en caso de que
    // sea necesario---
    // ini->setMoleculeTemp(new Molecule(ini->getMolVariable()));

    // -- Se establece el limite inferior de similitud. Si al division de los
    // volumenes de ambas moleculas es inferior a
    // -- a este limite, entonces no se ejecuta el algoritmo. De esta forma se
    // reduce el tiempo de busqueda
    if (indexArray[14] != -1) {
      ini->setLowLimit(atof(argv[indexArray[14]]));
    }
    // Comprobamos que el solapamiento entre ambas moleculas sean mayor que el
    // indicado por parametro
    double minMol =
        min(ini->getMolVariable()->tanimoto, ini->getMolQuery()->tanimoto);
    double maxMol =
        max(ini->getMolVariable()->tanimoto, ini->getMolQuery()->tanimoto);
    if ((minMol / maxMol) < ini->getLowLimit()) {
      cerr << "optipharm: -- Maximum posible score is lower than limit used: "
           << ini->getLowLimit() << "\n";
      continue;
    }

    // --- Get Longest and shortest axis from the two molecules and configuring
    // Upb y Lowb
    double *newLowb = new double[ini->Dimension()];
    double *newUpb = new double[ini->Dimension()];
    *newLowb = 0;
    *newUpb = 6.283185;

    // Set P1 value (0,0,0) or aleatory.
    Point3DDouble lower;
    Point3DDouble upper;
    VolumeOverlap::boundingBox(*(ini->getMolVariable()), &lower, &upper);

    // P1
    *(newLowb + 1) = lower.x;
    *(newLowb + 2) = lower.y;
    *(newLowb + 3) = lower.z;
    *(newUpb + 1) = upper.x;
    *(newUpb + 2) = upper.y;
    *(newUpb + 3) = upper.z;

    if (c.pointOneZero) {
      // P1
      *(newLowb + 1) = -0.5;
      *(newLowb + 2) = -0.5;
      *(newLowb + 3) = -0.5;
      *(newUpb + 1) = 0.5;
      *(newUpb + 2) = 0.5;
      *(newUpb + 3) = 0.5;
      // P2
      *(newLowb + 4) = -0.0001;
      *(newLowb + 5) = -0.0001;
      *(newLowb + 6) = -0.0001;
      *(newUpb + 4) = 0.0001;
      *(newUpb + 5) = 0.0001;
      *(newUpb + 6) = 0.0001;
    } else {
      // P2
      *(newLowb + 4) = lower.x;
      *(newLowb + 5) = lower.y;
      *(newLowb + 6) = lower.z;
      *(newUpb + 4) = upper.x;
      *(newUpb + 5) = upper.y;
      *(newUpb + 6) = upper.z;
    }

    double deltaXMin, deltaXMax, deltaYMin, deltaYMax, deltaZMin, deltaZMax;

    Parameters::defineDeltaMovement(ini, &deltaXMin, &deltaXMax, &deltaYMin,
                                    &deltaYMax, &deltaZMin, &deltaZMax);

    *(newLowb + 7) = deltaXMin;
    *(newLowb + 8) = deltaYMin;
    *(newLowb + 9) = deltaZMin;
    *(newUpb + 7) = deltaXMax;
    *(newUpb + 8) = deltaYMax;
    *(newUpb + 9) = deltaZMax;

    // Lower and Upper Bound
    /*
                            cout << "Lower\t Upper: \n";
                    cout.precision(17);
                            for (int i = 0; i < 10; i++) {
                    cout << *(newLowb + i) << "\t " << *(newUpb + i) << "\n";
                    }
    */
    // Se reconfiguran los limites con los nuevos valores
    ini->setLowb(newLowb);
    ini->setUpb(newUpb);
    delete[] newLowb;
    delete[] newUpb;
    ini->updateVectorIni();

    /*
    cout << "Threshold: "<< ini->getTreshhold()<<"\n";
    //Generacion de nuevas especies
             cout << "Niveles: "<< ini->Levels() << endl;
     cout << "Generacion de nuevas especies:\n";
     for (int i = 0; i < ini->Levels(); i++) {
     cout << ini->NewSpecEvals(i) << ", ";
     }
     cout<< endl;

    //Evaluaciones por nivel
    cout << "Evaluaciones por nivel:\n";
     for (int i = 0; i < ini->Levels(); i++) {
     cout << ini->Evals(i) << ", ";
     }
     cout <<endl;

     //Radios de cada nivel
     cout << "Radios:\n";
     for (int i = 0; i < ini->Levels(); i++) {
     cout << ini->R(i) << ", ";
     }
     cout<< endl;
     */
    // --- END ------------------------------------------------------------
    // ----GET EXECUTABLE NAME -------------------------------------------
    //---------------------------------------------------------------------
    // ---- Save name of files --------------------------------------------
    ini->setNameFileMolQuery(getFileName(argv[indexArray[1]]));
    ini->getMolQuery()->name_file = getFileName(argv[indexArray[1]]);
    ini->setNameFileMolVariable(getFileName(argv[indexArray[2]]));
    ini->getMolVariable()->name_file = getFileName(argv[indexArray[2]]);
    // ----- END -----------------------------------------------------------

    // --- SET save solution in mol2 files----------------------------------
    ini->setSaveSolutionMol2Files(c.saveMolFile);

    // ------ END -----------------------------

    if (ini == NULL || ini->Fail()) {
      message("Could not read ini file.", MSG_ERROR);
      errcode = 1;
    } else {
      message("Creating master.", MSG_INFORMATION);
      if (GetPars(argc - 2, argv + 2))
        master = new Master(ini, trace);

      if (master == NULL || master->Fail()) {
        message("Could not create master.", MSG_ERROR);
        errcode = 2;
      } else {
        // BORRAR
        /*for (int i = 0; i < INI.Levels(); i++) {
         cout << "Level: " << i << endl;
         cout << "Evals:" << ini->Evals(i) << endl;
         cout << "NewEvals: " << ini->NewSpecEvals(i) << endl;
         }*/
        message("Starting optimalization.", MSG_INFORMATION);
        // --- doing repcount experiments -----------------------
        for (long i = 1; i <= repcount; ++i) {

          // t1 = time(NULL);
          // SAIVNS: clock1 = clock();
          auto start = std::chrono::high_resolution_clock::now();
          // sprintf(msg, "Starting experiment %ld.", i);
          message(msg, MSG_INFORMATION);
          master->Go();
          // t2 = time(NULL);
          // SAINVS: clock2 = clock();
          auto end = std::chrono::high_resolution_clock::now();
          std::chrono::duration<double, std::ratio<1>> fp_ms = end - start;
          // tiemporeloj = (double) t2 - t1;
          tiemporeloj = fp_ms.count();
          /*tiemporeloj = (double) (clock2 - clock1)
           / (double) CLOCKS_PER_SEC;
           */
          if (master->Fail()) {
            errcode = 3;
            break;
          } else {
            message("Saving results.", MSG_INFORMATION);
            master->Save(stdout, tiemporeloj, saveflags);
            if (master->Fail()) {
              errcode = 4;
              break;
            };
          };
        };
      }
    };

    if (errcode == 0)
      message("Success.", MSG_INFORMATION);

    if (master != NULL)
      delete master;

    if (trace != NULL)
      delete trace;

    if (ini != NULL)
      delete ini;

    // --- END ITERAMOS EN CADA MOLECULA
  }

  if (inifile != NULL)
    fclose(inifile);
  if (indexArray != NULL)
    delete[] indexArray;
  // --- END I

  return errcode;
};
