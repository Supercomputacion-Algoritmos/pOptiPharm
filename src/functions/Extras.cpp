/*
 * Extras.cpp
 *
 *  Created on: Apr 21, 2018
 *      Author: savins
 */

#include "functions/Extras.h"

string Extras::getName() {
  /*

  string cmd = "git describe --tags";
  string data;
  FILE * stream;
  const int max_buffer = 256;
  char buffer[max_buffer];
  cmd.append(" 2>&1");

  stream = popen(cmd.c_str(), "r");
  if (stream) {
          while (!feof(stream))
                  if (fgets(buffer, max_buffer, stream) != NULL)
                          data.append(buffer);
          pclose(stream);
  }*/
  return "OptiPharm20210112";
}
