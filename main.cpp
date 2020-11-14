#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <memory>
#include <chrono>

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

static inline int d0Idx(const int x, const int y)
{
    return y - x;
}

static inline int d1Idx(const int x, const int y, const int n)
{
    return y + x - n;
}

static inline int conflicts(const int x, const int y, const int n,
                                const int * colHist, const int * d0Hist, const int * d1Hist)
{
    return colHist[x] + d0Hist[d0Idx(x, y)] + d1Hist[d1Idx(x, y, n)];
}

template<typename _URNG>
void initBoard(int * outColIdx, const int n, 
               int * colHist, int * d0Hist, int * d1Hist,
               MinList<int, int> & conflictList,
               _URNG & rng)
{
    // prepare histograms
    const int diagonalHistSize = 2 * n - 1;
    std::fill_n(colHist, n, 0);
    std::fill_n(d0Hist - n, diagonalHistSize, 0);
    std::fill_n(d1Hist - n, diagonalHistSize, 0);

    for (int y = 0; y < n; y++)
    {
        // reset list
        conflictList.Reset();

        // Evaluate all positions
        for (int x = 0; x < n; x++)
        {
            conflictList.Update(conflicts(x, y, n, colHist, d0Hist, d1Hist), x);
        }

        // select one of the best positions at random
        const int x = conflictList.SelectRandomValue(rng);
        outColIdx[y] = x;

        // update histograms
        colHist[x]++;
        d0Hist[d0Idx(x, y)]++;
        d1Hist[d1Idx(x, y, n)]++;
    }
}

template<typename _URNG>
bool minimizeConflicts(int * outColIdx, const int n, 
                       int * colHist, int * d0Hist, int * d1Hist,
                       MinList<int, int> & conflictList,
                       _URNG & rng)
{
    const int maxIters = 2 * n;
    for (int i = 0; i < maxIters; i++)
    {
        conflictList.Reset();

        // find queen with max conflicts
        for (int y = 0; y < n; y++)
        {
            conflictList.Update(conflicts(outColIdx[y], y, n, colHist, d0Hist, d1Hist), y, std::greater<int>());
        }

        std::pair<int, int> worstQueen = conflictList.SelectRandomPair(rng);
        const int worstConflicts = worstQueen.first;
        const int worstY = worstQueen.second;
        const int worstX = outColIdx[worstY];

        if (worstConflicts == 3)
        {
            // each queen counts itself three times
            return true;
        }

        // find new position with min conflicts
        conflictList.Reset();
        int x = 0;
        for (; x < worstX; x++) conflictList.Update(conflicts(x, worstY, n, colHist, d0Hist, d1Hist), x);
        // all other positions on the same row do not count the diagonal attacks - subtract them from the source X
        conflictList.Update(conflicts(worstX, worstY, n, colHist, d0Hist, d1Hist) - 2, x);
        x++;
        for (; x < worstX; x++) conflictList.Update(conflicts(x, worstY, n, colHist, d0Hist, d1Hist), x);

        // move queen and update histograms
        const int newX = conflictList.SelectRandomValue(rng);
        outColIdx[worstY] = newX;
        
        colHist[worstX]--;
        colHist[newX]++;

        d0Hist[d0Idx(worstX, worstY)]--;
        d0Hist[d0Idx(newX,   worstY)]++;

        d1Hist[d1Idx(worstX, worstY, n)]--;
        d1Hist[d1Idx(newX,   worstY, n)]++;
    }

    return false;
}

int main(int argc, char ** argv) 
{
    const std::string USAGE_MSG = "Usage: ./NQueens [N]";
    if (argc != 2)
    {
        std::cerr << USAGE_MSG << "\n";
        return 1;
    }

    const int n = std::stoi(argv[1]);

    // =======================
    // *** Allocate memory ***
    // =======================

    // board state
    std::unique_ptr<int[]> ptrQueenColIdx(new int[n]);

    // histograms to keep track of how many queens are where
    std::unique_ptr<int[]> ptrColHistogram(new int[n]);
    std::unique_ptr<int[]> ptrD0Histogram(new int[2 * n + 1]);
    std::unique_ptr<int[]> ptrD1Histogram(new int[2 * n + 1]);

    // extract raw pointers - avoid writing .get() everywhere
    int * queenColIdx = ptrQueenColIdx.get();
    int * colHist = ptrColHistogram.get();

    // make diagonal pointers start at the center of the histogram
    // this allows us to take advantage of negative indices
    int * d0Hist = ptrD0Histogram.get() + n + 1;
    int * d1Hist = ptrD1Histogram.get() + n + 1;

    // list used to choose min/max indices at random
    MinList<int, int> minMaxList(n);

    std::random_device rng;

    auto start = std::chrono::steady_clock::now();

    initBoard(queenColIdx, n, 
              colHist, d0Hist, d1Hist,
              minMaxList, rng);

    auto end = std::chrono::steady_clock::now();
    uint64_t elapsedUs = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    std::cout << "Init board: " << elapsedUs << " us\n";

    start = std::chrono::steady_clock::now();

    bool ok = minimizeConflicts(queenColIdx, n, 
                                colHist, d0Hist, d1Hist,
                                minMaxList, rng);

    end = std::chrono::steady_clock::now();
    elapsedUs = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    std::cout << "Minimize conflicts: " << elapsedUs << " us\n";

    std::cout << "Result: " << ok << "\n";
    std::cout << "Check: " << isSolution(queenColIdx, n) << "\n";

    printBoard(std::cout, queenColIdx, n);

    return 0;
}
