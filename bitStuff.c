//bitStuff.c
// basic bit manipulations

#include "bitStuff.h"

int setBit(int x, unsigned char position)
{
  int mask = 1 << position;
  return x | mask;
}

int clearBit(int x, unsigned char position)
{
  int mask = 1 << position;
  return x & ~mask;
}

int modifyBit(int x, unsigned char position, bool newState)
{
  int mask = 1 << position;
  int state = (int) newState; // relies on true = 1 and false = 0
  return (x & ~mask) | (-state & mask);
}

int flipBit(int x, unsigned char position)
{
  int mask = 1 << position;
  return x ^ mask;
}
 
bool isBitSet(int x, unsigned char position)
 {
   x >>= position;
   return (x & 1) != 0;
 }