#ifndef MIN_LIST_H
#define MIN_LIST_H

#include <memory>
#include <functional>
#include <random>

template<typename _Key, typename _Val>
class MinList
{
private:
    const size_t size;
    size_t insertIdx;
    std::unique_ptr<_Key[]> minKeys;
    std::unique_ptr<_Val[]> minValues;

public:
    MinList(const size_t n) : size(n)
    {
        minKeys = std::unique_ptr<_Key[]>(new _Key[n]);
        minValues = std::unique_ptr<_Val[]>(new _Val[n]);
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
    _Val SelectRandom(_URNG & rng)
    {
        std::uniform_int_distribution<int> idxDistribution(0, insertIdx - 1);
        auto idx = idxDistribution(rng);

        return minValues[idx];
    }
};

#endif