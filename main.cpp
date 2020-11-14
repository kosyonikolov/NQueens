#include <iostream>
#include <string>
#include <vector>

#include "printBoard.h"

int main(int argc, char ** argv) 
{
    const std::string USAGE_MSG = "Usage: ./NQueens [N]";
    if (argc != 2)
    {
        std::cerr << USAGE_MSG << "\n";
        return 1;
    }

    const int n = std::stoi(argv[1]);

    std::vector<int> colIdx(n);
    for (int i = 0; i < n; i++) colIdx[i] = i;

    printBoard(std::cout, colIdx);

    return 1;
}
