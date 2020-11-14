#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>

#include "MinList.h"
#include "printBoard.h"

bool isSolution(const int * queenColIdx, const int n)
{
    for (int y = 1; y < n; y++)
    {
        const int x = queenColIdx[y];
        
        for (int yPrev = y - 1; yPrev >= 0; yPrev--)
        {
            const int xPrev = queenColIdx[yPrev];

            if (xPrev == x) return false; // same row

            const int dx = x - xPrev; 
            const int dy = y - yPrev; // always positive
            if (std::abs(dx) == dy) return false; // diagonal
        }
    }

    return true;
}

void checkAllPermutations(const int n)
{
    std::vector<int> colIdx(n);
    for (int i = 0; i < n; i++) colIdx[i] = i;

    int solCount = 0;

    do
    {
        if (isSolution(colIdx.data(), n))
        {
            printBoard(std::cout, colIdx);
            std::cout << "\n";
            solCount++;
        }
    } 
    while (std::next_permutation(colIdx.begin(), colIdx.end()));
    
    std::cout << solCount << "\n";
}

void minListTest()
{
    std::vector<int> testVals = { 4, 4, 1, 2, 3, 1, 1, 3, -1, 3, 4, 4, 4, 4, 5 };
    MinList<int, int> minList(testVals.size());
    MinList<int, int> maxList(testVals.size());

    std::default_random_engine dre;

    for (int i = 0; i < testVals.size(); i++)
    {
        minList.Update(testVals[i], i);
        maxList.Update(testVals[i], i, std::greater<int>());

        std::cout << minList.SelectRandom(dre) << "\t" << maxList.SelectRandom(dre) << "\n";
    }
}

int main(int argc, char ** argv) 
{
    minListTest();

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
