#include "in.h"
#include "cstdlib"
#ifndef PINKNUMBER_H
#define PINKNUMBER_H


class PinkNumber
{
private:
  int max_key;
  int key;
  unsigned int white_values[5];
  unsigned int range;
public:
  PinkNumber();
  int GetNextValue();
};
#endif // PINKNUMBER_H
