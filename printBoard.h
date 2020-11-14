#ifndef PRINT_BOARD_H
#define PRINT_BOARD_H

#include <ostream>
#include <vector>

void printBoard(std::ostream & stream, const int * queenColIdx, const int n);

void printBoard(std::ostream & stream, const std::vector<int> & queenColIdx);

#endif