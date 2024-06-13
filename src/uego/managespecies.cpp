/*
 * ManageEspecies.cpp
 *
 *  Created on: 4 de abr. de 2016
 *      Author: savins
 */

#include "uegomanagespecies.h"
#include <cstdlib>
ManageSpecies::ManageSpecies() {
  start = NULL;
  temp1 = NULL;
  temp2 = NULL;
  temp3 = NULL;
}

ManageSpecies::~ManageSpecies() {
  delete start;
  delete temp1;
  delete temp2;
  delete temp3;
}

void ManageSpecies::addnode() // adding node
{
  char r;
  temp1 = new SpeciesList();
  temp1->level = rand() % 100;
  cout << "press 's' to add in start,'m' for midd , 'e' for end" << endl;
  cin >> r;
  switch (r) {
  case 's': // add start
    if (start == NULL) {
      start = temp1;
      temp1->next = NULL;
      temp1->prev = NULL;
    } else {
      temp2 = start;
      temp1->next = temp2;
      temp1->prev = NULL;
      start = temp1;
      temp2->prev = temp1;
    }
    break;
  case 'e': // add endif(start==NULL)

    if (start == NULL) {
      start = temp1;
      temp1->next = NULL;
      temp1->prev = NULL;
    } else {
      temp2 = start;
      while (temp2->next != NULL)
        temp2 = temp2->next;
      temp2->next = temp1;
      temp1->prev = temp2;
      temp1->next = NULL;
    }
    break;
  case 'm': // add midint num;
    int num;
    cout << "enter node after which you want to enter" << endl;
    cin >> num;
    temp2 = start;
    for (int i = 0; i < num; i++) {
      if (start == NULL)
        cout << "given node not found" << endl;
      else {
        temp3 = temp2;
        temp2 = temp2->next;
      }
    }
    temp1->next = temp2;
    temp3->next = temp1;
    temp1->prev = temp3;
    temp2->prev = temp1;
    break;
  }
}
void ManageSpecies::display() // displaying
{

  temp3 = start;
  if (start == NULL)
    cout << "no node to display" << endl;
  else {
    while (temp3->next != NULL) {
      cout << "Data stored is " << temp3->level << " at " << temp3 << endl;
      temp3 = temp3->next;
    }
    cout << "Data stored is " << temp3->level << " at " << temp3 << endl;
  }
}

void ManageSpecies::delnode() // deleting
{

  char d;
  cout << "press 's' to delete from start,'m' for midd , 'e' for end" << endl;
  cin >> d;
  switch (d) {
  case 's': // delete startif(start==NULL)

    if (start == NULL) {
      cout << "no node to delete" << endl;
    } else {
      temp1 = start;
      start = start->next;
      start->prev = NULL;
      delete temp1;
    }

    break;
  case 'e': // delete endif(start==NULL)
    if (start == NULL) {
      cout << "no node to delete" << endl;
    } else {
      temp1 = start;
      while (temp1->next != NULL) {
        temp2 = temp1;
        temp1 = temp1->next;
      }
      delete temp1;
      temp2->next = NULL;
    }

    break;
  case 'm': // delete mid
    int num;
    cout << "enter node you want to delete" << endl;
    cin >> num;

    temp1 = start;
    for (int i = 1; i < num; i++) {
      if (start == NULL)
        cout << "given node does not exist" << endl;
      else {
        temp2 = temp1;
        temp1 = temp1->next;
      }
    }
    temp3 = temp1->next;
    temp2->next = temp3;
    temp3->prev = temp2;
    delete temp1;
    break;
  }
}
void ManageSpecies::show() // backward display
{
  cout << "backward display" << endl;
  temp3 = start;
  if (start == NULL)
    cout << "no node to display" << endl;
  else {
    while (temp3->next != NULL) {
      temp3 = temp3->next;
    }
    while (temp3->prev != NULL) {
      cout << "Data stored is " << temp3->level << " at " << temp3 << endl;
      temp3 = temp3->prev;
    }
    cout << "Data stored is " << temp3->level << " at " << temp3 << endl;
  }
}
