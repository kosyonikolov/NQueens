#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>

#include "printBoard.h"

bool isSolution(const int * queenColIdx, const int n)
{
    for (int y = 1; y < n; y++)
    {
        const int x = queenColIdx[y];
        
        for (int yPrev = y - 1; yPrev >= 0; yPrev--)
        {
            const int xPrev = queenColIdx[yPrev];

            if (xPrev == x)
            {
                return false; // same row
            } 

            const int dx = x - xPrev; 
            const int dy = y - yPrev; // always positive

            if (std::abs(dx) == dy)
            {
                return false; // diagonal
            } 
        }
    }

    return true;
}

void checkAllPermutations(const int n)
{
    std::vector<int> colIdx(n);
    for (int i = 0; i < n; i++) colIdx[i] = i;

    do
    {
        if (isSolution(colIdx.data(), n))
        {
            printBoard(std::cout, colIdx);
            std::cout << "\n";
        }
    } 
    while (std::next_permutation(colIdx.begin(), colIdx.end()));
    
}

int main(int argc, char ** argv) 
{
    // int test[8] = {5, 3, 6, 0, 7, 1, 4, 2};
    // printBoard(std::cout, test, 8);

    // std::cout << "Solution = " << isSolution(test, 8) << "\n";

    const std::string USAGE_MSG = "Usage: ./NQueens [N]";
    if (argc != 2)
    {
        std::cerr << USAGE_MSG << "\n";
        return 1;
    }

    const int n = std::stoi(argv[1]);

    checkAllPermutations(n);

    return 0;
}