#pragma once
#include "Typedefs.h"


class SmoothRankine
{
public:
    static double Get(const StressVector& rStress)
    {
        double smoothRankine = rStress[0]; // TODO
        return smoothRankine;
    }

    static Derivative GetDerivative(const StressVector& rStress)
    {
        Derivative smoothRankineDerivative = rStress; // TODO
        return smoothRankineDerivative;
    }
};
