#include "Typedefs.h"
#include <eigen3/Eigen/Dense>


class SmoothRankineEigen
{
public:
    //! @brief calculates the smooth rankine norm via Eigen::eigenvalues()
    //! @param rStress stress
    //! @return smooth rankine norm
    static double Get(const StressVector& rStress)
    {
        auto eigenvalues = CalculateTensor(rStress).eigenvalues().real();
        double eigenvalueSum = eigenvalues.sum();
        return std::sqrt(eigenvalueSum);
    }

    //! @brief calculates the derivative of the smooth rankine norm with respect to the stress in voigt notation via
    //! central
    //! differences
    //! @param rStress stress
    //! @param rDelta delta of the central differences
    //! @return derivative
    static Derivative GetDerivative(const StressVector& rStress, double rDelta)
    {
        Derivative CDF;
        for (int i = 0; i < rStress.rows(); ++i)
        {
            StressVector sMinus = rStress;
            StressVector sPlus = rStress;
            sMinus[i] -= 0.5 * rDelta;
            sPlus[i] += 0.5 * rDelta;
            CDF[i] = (Get(sPlus) - Get(sMinus)) / rDelta;
        }
        return CDF;
    }

    //! @brief converts the a vector in voigt notation to a tensor
    static StressTensor CalculateTensor(const StressVector& rStress)
    {
        StressTensor s;

        s(0, 0) = rStress[0];
        s(1, 1) = rStress[1];
        s(2, 2) = rStress[2];

        s(1, 0) = rStress[5];
        s(0, 1) = rStress[5];

        s(2, 0) = rStress[4];
        s(0, 2) = rStress[4];

        s(2, 1) = rStress[3];
        s(1, 2) = rStress[3];

        return s;
    }
};
