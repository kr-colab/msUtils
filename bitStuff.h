#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

int setBit(int x, unsigned char position);
int clearBit(int x, unsigned char position);
int modifyBit(int x, unsigned char position, bool newState);
int flipBit(int x, unsigned char position);
bool isBitSet(int x, unsigned char position);
