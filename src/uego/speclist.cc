#include "uego/uego.h"

////////////////////////////////////////////////////////////
// $Id: speclist.cc,v 2.6 1998/03/29 10:41:12 jelasity Exp $
// speclist.cc
// defines for class SpeciesList
////////////////////////////////////////////////////////////
// modification history:
//
////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------

long SpeciesList::Optimize(long maxevals) { //, short tmplevel){
  long l;
  // level=tmplevel;
  if (maxevals == -1) {
    // maxevals = Evals(level);
    maxevals = Evals(Master::ini().getLevelMaster());

    if (Fail())
      return 0;
  };
  // cout << "LEVEL especie: "<< level <<", LEVEL MASTER " <<
  // Master::ini().getLevelMaster()<< " Maxevals: " << maxevals<<endl;
  l = center->Optimize(level, maxevals);
  // l = center->Optimize(Master::ini().getLevelMaster(), maxevals);

  FailFlag = center->Fail();

  return l;
};

// -----------------------------------------------------------------------

long SpeciesList::Evals(long _level) {

  FailFlag = 1 == 1;

  if (!Master::iniSet()) {
    message("Ini not set in SpeciesList::Evals()", MSG_ERROR);
    return 0;
  };

  if (_level == -1)
    _level = level; // called without param
  if (_level == 0)
    _level = 1;

  FailFlag = 1 == 0;

  return Master::ini().Evals(_level) / (long)Master::M(_level);
};

// -----------------------------------------------------------------------

long SpeciesList::NewSpecies(SpeciesList **result, short newlevel,
                             long maxevals) {

  SpeciesList *head = NULL, *tmp;
  long oldlevel = level, evals = 0, maxev1 = 0, maxev2 = 0;

  if (level == newlevel - 1) // the old way if lowest level
    return _NewSpecies(result, newlevel, maxevals);

  if (maxevals < 6)
    maxev1 = maxevals; // maxev < 3 makes no sense
  else
    maxev1 = maxev2 = maxevals / 2;

  // --- creating species at lowest possible level
  level = newlevel - 1;
  evals += _NewSpecies(result, newlevel, maxev1);
  if (Fail())
    return evals;

  // --- creating species at original level
  level = oldlevel;
  evals += _NewSpecies(&head, newlevel, maxev2);
  if (Fail())
    return evals;

  // --- concatenating
  if (*result == NULL)
    *result = head;
  else {
    for (tmp = *result; tmp->next != NULL; tmp = tmp->next)
      ;
    tmp->next = head;
    if (head != NULL)
      head->prev = tmp;
  };

  return evals;
};
// -----------------------------------------------------------------------

long SpeciesList::Mutate(SpeciesList **result, short newlevel,
                         long nOfMutations) {

  SpeciesList *head = NULL, *tmp;
  long oldlevel = level, evals = 0, maxev1 = 0, maxev2 = 0;

  evals += _Mutate(result, newlevel, nOfMutations);

  // --- concatenating
  if (*result == NULL)
    *result = head;
  else {
    for (tmp = *result; tmp->next != NULL; tmp = tmp->next)
      ;
    tmp->next = head;
    if (head != NULL)
      head->prev = tmp;
  };

  return evals;
};

// -----------------------------------------------------------------------

long SpeciesList::_NewSpecies(SpeciesList **result, short newlevel,
                              long maxevals) {
  // --- length is counted so that approximately
  // ---     length + (length*(length-1))/2 = maxevals
  // --- where (length*(length-1))/2 is the number of pairs in length elements

  SearchSpElement **base, *oldcenter, *tmp;
  long length, evals = 0, i, j;
  char *toresult;
  SpeciesList reshead, // temporary head for result list
      *res;

  *result = NULL;

  FailFlag = 1 == 0;
  // cout << "Maxevals" << maxevals<<endl;
  // getchar();
  if (maxevals <= 0)
    return 0;

  FailFlag = 1 == 1;

  // --- initializing base list with random elements from attr. area
  length = (long)((sqrt(1.0 + 8 * maxevals) - 1.0) / 2.0);
  oldcenter = center;

  // cout << length<<endl;
  // getchar();
  base = new SearchSpElement *[length];
  toresult = new char[length];
  i = 0;
  if (base != NULL)
    for (; i < length; ++i) {
      base[i] = INI.Prototype()->RandNew(level, oldcenter);

      if (base[i] == NULL)
        break;
      ++evals;
      if (base[i]->CurrValue() > center->CurrValue())
        center->UpdateFrom(base[i]);
    };
  if (i != length) {
    message("No memory (NewSpecies,1)", MSG_ERROR);
    goto CLEAN;
  };

  // --- creating new species -------------------------------------
  for (i = 0; i < length; ++i)
    toresult[i] = 1 == 0;

  for (i = 0; i < length - 1; ++i)
    for (j = i + 1; j < length; ++j) {
      if (toresult[i] && toresult[j])
        continue;
      tmp = INI.Prototype()->BetweenNew(base[i], base[j]);
      if (tmp == NULL) {
        message("No memory (NewSpecies,2)", MSG_ERROR);
        goto CLEAN;
      };
      ++evals;
      if (tmp->CurrValue() > center->CurrValue())
        center->UpdateFrom(tmp);
      // else if (tmp->CurrValue() < base[i]->CurrValue()
      //		&& tmp->CurrValue() < base[j]->CurrValue())
      //  --- new species -------------------------------
      toresult[i] = toresult[j] = 1 == 1;
      delete tmp;
    };

  // --- collecting result ---------------------------------------
  res = &reshead;
  for (i = 0; i < length; ++i) {
    if (toresult[i]) {
      res->next = new SpeciesList(base[i], newlevel);
      if (res->next == NULL) {
        message("No memory (NewSpecies,3)", MSG_ERROR);
        for (; i < length; ++i)
          toresult[i] = 1 == 0;
        goto CLEAN;
      };
      res->next->prev = res;
      res = res->next;
    };
  };
  *result = reshead.next;
  if (*result != NULL)
    (*result)->prev = NULL;
  reshead.next = NULL; // to prevent result from destructing

  FailFlag = 1 == 0;

CLEAN:
  if (base != NULL) {
    for (i = 0; i < length && base[i] != NULL; ++i)
      if (!toresult[i])
        delete base[i];
    delete[] base;
  };
  if (toresult != NULL)
    delete[] toresult;

  return evals;
};

// -----------------------------------------------------------------------

long SpeciesList::_Mutate(SpeciesList **result, short newlevel,
                          long nOfMutations) {
  // --- length is counted so that approximately
  // ---     length + (length*(length-1))/2 = maxevals
  // --- where (length*(length-1))/2 is the number of pairs in length elements

  SearchSpElement **base, *oldcenter, *tmp;
  long length, evals = 0, i, j;
  char *toresult;
  SpeciesList reshead, // temporary head for result list
      *res;

  *result = NULL;

  FailFlag = 1 == 0;
  // cout << "Maxevals" << maxevals<<endl;
  // getchar();

  FailFlag = 1 == 1;

  // --- initializing base list with random elements from attr. area
  length = nOfMutations;
  oldcenter = center;

  // cout << length<<endl;
  // getchar();
  base = new SearchSpElement *[length];
  toresult = new char[length];
  i = 0;
  if (base != NULL)
    for (; i < length; ++i) {
      base[i] = INI.Prototype()->generateMutationNew(newlevel, oldcenter);

      if (base[i] == NULL)
        break;
      ++evals;
      // if (base[i]->CurrValue() > center->CurrValue())
      //	center->UpdateFrom(base[i]);
    };
  if (i != length) {
    message("No memory (MutateNew,1)", MSG_ERROR);
    goto CLEAN;
  };

  // --- creating new species -------------------------------------
  for (i = 0; i < length; ++i)
    toresult[i] = 1 == 1;

  // --- collecting result ---------------------------------------
  res = &reshead;
  for (i = 0; i < length; ++i) {
    if (toresult[i]) {
      res->next = new SpeciesList(base[i], newlevel);
      if (res->next == NULL) {
        message("No memory (NewSpecies,3)", MSG_ERROR);
        for (; i < length; ++i)
          toresult[i] = 1 == 0;
        goto CLEAN;
      };
      res->next->prev = res;
      res = res->next;
    };
  };
  *result = reshead.next;
  if (*result != NULL)
    (*result)->prev = NULL;
  reshead.next = NULL; // to prevent result from destructing

  FailFlag = 1 == 0;

CLEAN:
  if (base != NULL) {
    for (i = 0; i < length && base[i] != NULL; ++i)
      if (!toresult[i])
        delete base[i];
    delete[] base;
  };
  if (toresult != NULL)
    delete[] toresult;

  return evals;
};

// -----------------------------------------------------------------------

void SpeciesList::Save(FILE *stream, bool best) {

  // FailFlag = fprintf( stream, "# level: %d\n", level ) == EOF;
  if (!FailFlag) {
    center->Save(stream, best);
    if (!center->Fail())
      return; // FailFlag is properly set!
  }

  FailFlag = 1 == 1;
  message("Error writing species to stream.", MSG_ERROR);
  return;
};

void SpeciesList::Save(FILE *stream) {

  // FailFlag = fprintf( stream, "# level: %d\n", level ) == EOF;
  if (!FailFlag) {
    center->Save(stream);
    if (!center->Fail())
      return; // FailFlag is properly set!
  }

  FailFlag = 1 == 1;
  message("Error writing species to stream.", MSG_ERROR);
  return;
};

// -----------------------------------------------------------------------

void SpeciesList::Absorb(SpeciesList *sl) {

  if (sl->center->CurrValue() > center->CurrValue())
    center->UpdateFrom(sl->center);
};
