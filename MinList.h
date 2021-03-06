#ifndef MIN_LIST_H
#define MIN_LIST_H

#include <vector>
#include <functional>
#include <random>
#include <utility>

template<typename _Key, typename _Val>
class MinList
{
private:
    const size_t size;
    size_t insertIdx;
    std::vector<_Key> minKeys;
    std::vector<_Val> minValues;

    template<class _URNG>
    int SelectRandomIndex(_URNG & rng)
    {
        std::uniform_int_distribution<int> idxDistribution(0, insertIdx - 1);
        return idxDistribution(rng);
    }

public:
    MinList(const size_t n) : size(n)
    {
        minKeys = std::vector<_Key>(n);
        minValues = std::vector<_Val>(n);
    }

    void Reset()
    {
        insertIdx = 0;
    }

    template<typename _Comparator = std::less<_Key>>
    void Update(const _Key key, const _Val value, _Comparator comp = _Comparator())
    {
        if (insertIdx == 0)
        {
            // empty list - place into the first element
            minKeys[0] = key;
            minValues[0] = value;
            insertIdx++;
            return;
        }

        if (comp(minKeys[0], key)) return; // key is larger than min

        if (minKeys[0] == key)
        {
            // new key is the same value as the current minimum - add it to the list
            minKeys[insertIdx] = key;
            minValues[insertIdx] = value;
            insertIdx++;
        }
        else // key is smaller than current min
        {
            minKeys[0] = key;
            minValues[0] = value;
            insertIdx = 1;
        }
    }

    template<class _URNG>
    _Val SelectRandomValue(_URNG & rng)
    {
        return minValues[SelectRandomIndex(rng)];
    }

    template<class _URNG>
    std::pair<_Key, _Val> SelectRandomPair(_URNG & rng)
    {
        const auto idx = SelectRandomIndex(rng);
        return std::pair<_Key, _Val>(minKeys[idx], minValues[idx]);
    }
};

#endif