#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <memory>
#include <chrono>

#include "MinList.h"
#include "printBoard.h"

#define PRINT_STATS (1)
#define USE_HORSE_INIT (1)

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
    // split up for debug
    int colContrib = colHist[x];
    int _d0Idx = d0Idx(x, y);
    int d0Contrib = d0Hist[_d0Idx];
    int _d1Idx = d1Idx(x, y, n);
    int d1Contrib = d1Hist[_d1Idx];

    return colContrib + d0Contrib + d1Contrib;

    //return colHist[x] + d0Hist[d0Idx(x, y)] + d1Hist[d1Idx(x, y, n)];
}

void clearHistograms(int * colHist, int * d0Hist, int * d1Hist, const int n)
{
    const int diagonalHistSize = 2 * n + 1;
    std::fill_n(colHist, n, 0);
    std::fill_n(d0Hist - n, diagonalHistSize, 0);
    std::fill_n(d1Hist - n, diagonalHistSize, 0);
}

template<typename _URNG>
void initBoard(int * outColIdx, const int n, 
               int * colHist, int * d0Hist, int * d1Hist,
               MinList<int, int> & conflictList,
               _URNG & rng)
{
#if PRINT_STATS
    auto start = std::chrono::steady_clock::now();
#endif

    // prepare histograms
    clearHistograms(colHist, d0Hist, d1Hist, n);

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

#if PRINT_STATS
    auto end = std::chrono::steady_clock::now();
    const uint64_t elapsedUs = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "[InitBoard] " << elapsedUs << " ms\n";
#endif
}

template<typename _URNG>
void horsesForCourses(int * outColIdx, const int n, 
                      int * colHist, int * d0Hist, int * d1Hist,
                      const float pRandom, const int randomAttemps,
                      MinList<int, int> & conflictList,
                      _URNG & rng)
{
#if PRINT_STATS
    auto start = std::chrono::steady_clock::now();
#endif

    clearHistograms(colHist, d0Hist, d1Hist, n);

    std::uniform_int_distribution<int> xDist(0, n - 1);
    std::uniform_real_distribution<float> pDist(0, 1);

    auto placeQueen = [&](const int x, const int y)
    {
        outColIdx[y] = x;
        colHist[x]++;
        d0Hist[d0Idx(x, y)]++;
        d1Hist[d1Idx(x, y, n)]++;
    };

    int y = 0;
    // place first queen
    {
        const int x = xDist(rng);
        placeQueen(x, y);
        y++;
    }

    // place second queen
    {
        conflictList.Reset();
        
        for (int i = 0; i < randomAttemps; i++)
        {
            const int x = xDist(rng);
            conflictList.Update(conflicts(x, y, n, colHist, d0Hist, d1Hist), x);
        }

        const int xPrev = outColIdx[y - 1];
        if (xPrev >= 2) conflictList.Update(conflicts(xPrev - 2, y, n, colHist, d0Hist, d1Hist), xPrev - 2);
        else conflictList.Update(conflicts(n + xPrev - 2, y, n, colHist, d0Hist, d1Hist), n + xPrev - 2);
        if (xPrev < n - 2) conflictList.Update(conflicts(xPrev + 2, y, n, colHist, d0Hist, d1Hist), xPrev + 2);
        else conflictList.Update(conflicts(xPrev + 2 - n, y, n, colHist, d0Hist, d1Hist), xPrev + 2 - n);

        const int x = conflictList.SelectRandomValue(rng);
        placeQueen(x, y);
        y++;
    }

    for (; y < n; y++)
    {
        conflictList.Reset();

        // positions from queen two lines above
        {
            const int xPrev = outColIdx[y - 2];
            if (xPrev >= 1) conflictList.Update(conflicts(xPrev - 1, y, n, colHist, d0Hist, d1Hist), xPrev - 1);
            else conflictList.Update(conflicts(n + xPrev - 1, y, n, colHist, d0Hist, d1Hist), n + xPrev - 1);
            if (xPrev < n - 1) conflictList.Update(conflicts(xPrev + 1, y, n, colHist, d0Hist, d1Hist), xPrev + 1);
            else conflictList.Update(conflicts(xPrev + 1 - n, y, n, colHist, d0Hist, d1Hist), xPrev + 1 - n);
        }

        // positions from queen one line above
        {
            const int xPrev = outColIdx[y - 1];
            if (xPrev >= 2) conflictList.Update(conflicts(xPrev - 2, y, n, colHist, d0Hist, d1Hist), xPrev - 2);
            else conflictList.Update(conflicts(n + xPrev - 2, y, n, colHist, d0Hist, d1Hist), n + xPrev - 2);
            if (xPrev < n - 2) conflictList.Update(conflicts(xPrev + 2, y, n, colHist, d0Hist, d1Hist), xPrev + 2);
            else conflictList.Update(conflicts(xPrev + 2 - n, y, n, colHist, d0Hist, d1Hist), xPrev + 2 - n);
        }

        // random positions
        if (pDist(rng) <= pRandom)
        {
            for (int i = 0; i < randomAttemps; i++)
            {
                const int x = xDist(rng);
                conflictList.Update(conflicts(x, y, n, colHist, d0Hist, d1Hist), x);
            }
        }

        const int x = conflictList.SelectRandomValue(rng);
        placeQueen(x, y);
    }

#if PRINT_STATS
    auto end = std::chrono::steady_clock::now();
    const uint64_t elapsedUs = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "[HorsesForCourses] " << elapsedUs << " ms\n";
#endif
}

template<typename _URNG>
bool minimizeConflicts(int * outColIdx, const int n, 
                       int * colHist, int * d0Hist, int * d1Hist,
                       MinList<int, int> & conflictList,
                       const int maxIters, _URNG & rng)
{
#if PRINT_STATS
    auto start = std::chrono::steady_clock::now();
#endif

    int i = 0;
    bool ok = false;
    for (; i < maxIters; i++)
    {
        conflictList.Reset();

        // find queen with max conflicts
        for (int y = 0; y < n; y++)
        {
            const int x = outColIdx[y];
            const int currentConflicts = conflicts(x, y, n, colHist, d0Hist, d1Hist);
            conflictList.Update(currentConflicts, y, std::greater<int>());
        }

        std::pair<int, int> worstQueen = conflictList.SelectRandomPair(rng);
        const int worstConflicts = worstQueen.first;
        const int worstY = worstQueen.second;
        const int worstX = outColIdx[worstY];

        if (worstConflicts == 3)
        {
            // each queen counts itself three times
            ok = true;
            break;
        }

        // find new position with min conflicts
        conflictList.Reset();
        int x = 0;
        for (; x < worstX; x++)
        {
            const int candConf = conflicts(x, worstY, n, colHist, d0Hist, d1Hist);
            conflictList.Update(candConf, x);
        } 
        // source position - we will count 3 attacks from ourselves (col and two diagonals)
        const int srcConf = conflicts(worstX, worstY, n, colHist, d0Hist, d1Hist) - 3;
        conflictList.Update(srcConf, x);
        x++;
        for (; x < n; x++)
        {
            const int candConf = conflicts(x, worstY, n, colHist, d0Hist, d1Hist);
            conflictList.Update(candConf, x);
        } 
        // std::cout << "\n";

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

#if PRINT_STATS
    auto end = std::chrono::steady_clock::now();
    const uint64_t elapsedUs = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "[MinConflicts] " << i << " iters / " << elapsedUs << " ms\n";
#endif
    
    return ok;
}

bool solve(int * outColIdx, const int n,
           const int maxStarts, const float maxIterationsRatio)
{
    // **************** Allocate memory ****************

    // histograms to keep track of how many queens are where
    std::unique_ptr<int[]> ptrColHistogram(new int[n]);
    std::unique_ptr<int[]> ptrD0Histogram(new int[2 * n + 1]);
    std::unique_ptr<int[]> ptrD1Histogram(new int[2 * n + 1]);

    // extract raw pointers - avoid writing .get() everywhere
    int * colHist = ptrColHistogram.get();

    // make diagonal pointers start at the center of the histogram
    // this allows ms to take advantage of negative indices
    int * d0Hist = ptrD0Histogram.get() + n;
    int * d1Hist = ptrD1Histogram.get() + n;

    // list msed to choose min/max indices at random
    MinList<int, int> minMaxList(n);

    std::random_device rng;
    const int maxMinConflictIters = std::round(n * maxIterationsRatio);

#if USE_HORSE_INIT
    const float pRandom = 1.0f;
    const float attemptMult = 2.0f;
    const int randomAttempts = std::max(10, static_cast<int>(std::round(attemptMult * std::log2(n))));

    horsesForCourses(outColIdx, n,
                     colHist, d0Hist, d1Hist,
                     pRandom, randomAttempts,
                     minMaxList, rng);
#else   
    initBoard(outColIdx, n, 
              colHist, d0Hist, d1Hist,
              minMaxList, rng);
#endif

    int iStart = 1;
    bool ok = false;
    for (; iStart <= maxStarts; iStart++)
    {
        ok = minimizeConflicts(outColIdx, n, 
                               colHist, d0Hist, d1Hist,
                               minMaxList, 
                               maxMinConflictIters, rng);

        if (ok) break;

#if USE_HORSE_INIT
        horsesForCourses(outColIdx, n,
                        colHist, d0Hist, d1Hist,
                        pRandom, randomAttempts,
                        minMaxList, rng);
#else   
        initBoard(outColIdx, n, 
                colHist, d0Hist, d1Hist,
                minMaxList, rng);
#endif
    }

#if PRINT_STATS
    std::cout << "Starts: " << iStart << "\n";
#endif

    return ok;
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
    // board state
    std::unique_ptr<int[]> ptrQueenColIdx(new int[n]);
    int * queenColIdx = ptrQueenColIdx.get();
    
    auto start = std::chrono::steady_clock::now();

    bool ok = solve(queenColIdx, n, 10, 1.0f);

    auto end = std::chrono::steady_clock::now();
    uint64_t solTimeMs = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Solution time: " << solTimeMs << " ms\n";

    if (ok && n <= 42)
    {
        printBoard(std::cout, queenColIdx, n);
    }

    std::cout << "Result: " << ok << "\n";
    std::cout << "Check: " << isSolution(queenColIdx, n) << "\n";

    return 0;
}
