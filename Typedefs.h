#pragma once
#include <eigen3/Eigen/Core>

using StressVector = Eigen::Matrix<double, 6, 1>;
using StressTensor = Eigen::Matrix<double, 3, 3>;
using Derivative = Eigen::Matrix<double, 6, 1>;

typedef struct yieldAndDerivative
{
    double yieldFunction;
    Derivative deriv;
} YieldAndDerivative;
