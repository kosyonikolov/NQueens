#include "printBoard.h"

void printBoard(std::ostream & stream, const int * queenColIdx, const int n)
{
    for (int y = 0; y < n; y++)
    {
        const int xQueen = queenColIdx[y];
        int colIdx = 0;
        for (; colIdx < xQueen; colIdx++) stream << "- ";
        stream << "X ";
        colIdx++;
        for (; colIdx < n; colIdx++) stream << "- ";
        stream << "\n";
    }
}

void printBoard(std::ostream & stream, const std::vector<int> & queenColIdx)
{
    printBoard(stream, queenColIdx.data(), queenColIdx.size());
}