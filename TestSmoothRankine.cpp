#define BOOST_TEST_MODULE DummyTestName
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <eigen3/Eigen/Dense>
#include "SmoothRankine.h"
#include "SmoothRankineEigen.h"

#include <iostream>

//! @brief compares the smooth rankine norm values (Eigen vs own implementation)
//! @param rStress vector of test stresses
void CheckValues(const std::vector<StressVector>& rStresses)
{
    for (const auto& stress : rStresses)
        BOOST_CHECK_CLOSE(SmoothRankineEigen::Get(stress), SmoothRankine::Get(stress), 1.e-6);
}

//! @brief compares the smooth rankine norm derivatives (Eigen vs own implementation)
//! @param rStress vector of test stresses
void CheckDerivatives(const std::vector<StressVector>& rStresses)
{
    for (const auto& stress : rStresses)
        BOOST_CHECK(SmoothRankineEigen::GetDerivative(stress, 1.e-6).isApprox(SmoothRankine::GetDerivative(stress)));
}

//! @param rNum number of test stresses
//! @return rNum random test stresses
std::vector<StressVector> GetRandomStresses(int rNum)
{
    std::vector<StressVector> stresses(rNum);
    for (auto& stress : stresses)
        stress = StressVector::Random();
    return stresses;
}

//! @return vector of *hand picked* test stresses
std::vector<StressVector> GetCarefullyChosenStresses()
{
    std::vector<StressVector> stresses;
    stresses.push_back((Eigen::VectorXd(6) << 1, 0, 0, 0, 0, 0).finished());
    stresses.push_back((Eigen::VectorXd(6) << 0, 1, 0, 0, 0, 0).finished());
    stresses.push_back((Eigen::VectorXd(6) << 0, 0, 1, 0, 0, 0).finished());
    stresses.push_back((Eigen::VectorXd(6) << 0, 0, 0, 1, 0, 0).finished());
    stresses.push_back((Eigen::VectorXd(6) << 0, 0, 0, 0, 1, 0).finished());
    stresses.push_back((Eigen::VectorXd(6) << 0, 0, 0, 0, 0, 1).finished());
    // please come up with some better ones :)
    return stresses;
}


BOOST_AUTO_TEST_CASE(SmoothRankineValue)
{
    CheckValues(GetCarefullyChosenStresses());
}

BOOST_AUTO_TEST_CASE(SmoothRankineValueRandom)
{
    CheckValues(GetRandomStresses(100));
}

BOOST_AUTO_TEST_CASE(SmoothRankineDerivative)
{
    CheckDerivatives(GetCarefullyChosenStresses());
}

BOOST_AUTO_TEST_CASE(SmoothRankineDerivativeRandom)
{
    CheckDerivatives(GetRandomStresses(100));
}
