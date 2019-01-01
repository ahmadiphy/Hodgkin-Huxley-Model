#include "pinknumber.h"

PinkNumber::PinkNumber()
{
    max_key = 0x1f; // Five bits set
    this->range = range;
    key = 0;
    for (int i = 0; i < 5; i++)
        white_values[i] = rand() % (range/5);
}
int PinkNumber::GetNextValue()
{
    int last_key = key;
    unsigned int sum;

    key++;
    if (key > max_key)
        key = 0;
    // Exclusive-Or previous value with current value. This gives
    // a list of bits that have changed.
    int diff = last_key ^ key;
    sum = 0;
    for (int i = 0; i < 5; i++)
    {
        // If bit changed get new random number for corresponding
        // white_value
        if (diff & (1 << i))
            white_values[i] = rand() % (range/5);
        sum += white_values[i];
    }
    return sum;
}
