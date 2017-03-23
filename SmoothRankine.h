#pragma once
#include "Typedefs.h"
#include <iostream>


const double pi = 3.1415926535897932384626433832795028841971;
const double eps = std::numeric_limits<double>::epsilon();

class SmoothRankine
{
public:
    static double Get(const StressVector& rStress)
    {
//        double smoothRankine = rStress[0]; // TODO

        YieldAndDerivative yad = CalculateYieldFunctionAndDerivative(rStress);
        double smoothRankine = yad.yieldFunction;

        return smoothRankine;
    }

    static Derivative GetDerivative(const StressVector& rStress)
    {
//        Derivative smoothRankineDerivative = rStress; // TODO
        YieldAndDerivative yad = CalculateYieldFunctionAndDerivative(rStress);
        Derivative smoothRankineDerivative = yad.deriv;
        return smoothRankineDerivative;
    }


private:

    //!
    //! \brief calculates the yield surface of rankine type and corresponding derivative with respect to the stress in
    //! voigt notation
    //! \param rStress stress
    //! \return structure containing the yield surface and derivative
    //!
    static YieldAndDerivative CalculateYieldFunctionAndDerivative(const StressVector& rStress)
    {
        StressVector stressPermuted;
        stressPermuted[0] = rStress[0];
        stressPermuted[1] = rStress[1];
        stressPermuted[2] = rStress[2];
        stressPermuted[3] = rStress[5];
        stressPermuted[4] = rStress[3];
        stressPermuted[5] = rStress[4];


        //-------------------- important variables --------------------
        double a = (-1)*(stressPermuted[0] + stressPermuted[1] + stressPermuted[2]);
        double b = stressPermuted[0]*stressPermuted[1] + stressPermuted[0]*stressPermuted[2] + stressPermuted[1]*stressPermuted[2]
                - std::pow(stressPermuted[3],2) - std::pow(stressPermuted[4],2) - std::pow(stressPermuted[5], 2);
        double c = stressPermuted[0]*std::pow(stressPermuted[4], 2) + stressPermuted[1]*std::pow(stressPermuted[5], 2) + stressPermuted[2]*std::pow(stressPermuted[3], 2)
                - stressPermuted[0]*stressPermuted[1]*stressPermuted[2] - 2*stressPermuted[3]*stressPermuted[4]*stressPermuted[5];

        double q = std::pow(a,3)/27 - a*b/6 + c/2;
        double p = b/3 - std::pow(a,2)/9;

        double D = std::pow(p,3) + std::pow(q,2);

        double P = std::sqrt(std::abs(p));
        P = (q < 0 ? (-1)*P : P);
        double beta = 1/3*std::acos(q/(std::pow(P,3)));


        //-------------------- eigenvalues --------------------
        Eigen::Vector3d eigenValues;
        eigenValues[0] = -2*P*cos(beta)-a/3;
        eigenValues[0] = (greaterThan(eigenValues[0], 0.) ? eigenValues[0] : 0);

        eigenValues[1] = 2*P*cos(beta+pi/3) - a/3;
        eigenValues[1] = (greaterThan(eigenValues[1], 0.) ? eigenValues[1] : 0);

        eigenValues[2] = 2*P*cos(beta-pi/3) - a/3;
        eigenValues[2] = (greaterThan(eigenValues[2], 0.) ? eigenValues[2] : 0);

        double f = 0.;


        //-------------------- yield function --------------------
        //just one positive eigenvalue
        if (greaterThan(eigenValues[0], 0.) && equals(eigenValues[1], 0.))
        {
            f = eigenValues[0];
        }
        //two positive, identical eigenvalues
        else if (equals(eigenValues[0], 0.) && greaterThan(eigenValues[1], 0.))
        {
            f = std::sqrt(eigenValues.transpose()*eigenValues);
        }
        //three positive eigenvalues
        else if (greaterThan(eigenValues[0], 0.) && greaterThan(eigenValues[0], 0.))
        {
            f = std::sqrt(std::pow(a,2) - 2*b);
        }


        //-------------------- derivative --------------------
        Derivative Da_Ds;
        Derivative Db_Ds;
        Derivative Dc_Ds;
        Derivative Dp_Ds;
        Derivative Dq_Ds;
        Derivative DP_Ds;
        Derivative Dbeta_Ds;

        double Df_Dq;
        double Df_DP;
        double Df_Da;
        double Df_Db;

        double Dsig_1_Dq;
        double Dsig_1_DP;
        double Dsig_1_Da;

        Derivative Dsig_1_Ds;
        Derivative Dsig_2_Ds;
        Derivative Dsig_3_Ds;

        Derivative Df_Ds;

        Da_Ds << -1, -1, -1, 0, 0, 0;
        Db_Ds << stressPermuted[1]+stressPermuted[2], stressPermuted[0]+stressPermuted[2], stressPermuted[0]+stressPermuted[1], -2*stressPermuted[3], -2*stressPermuted[4], -2*stressPermuted[5];
        Dc_Ds << std::pow(stressPermuted[4],2)-stressPermuted[1]*stressPermuted[2], std::pow(stressPermuted[5],2)-stressPermuted[0]*stressPermuted[2], std::pow(stressPermuted[3], 2)-stressPermuted[0]*stressPermuted[1],
                2*(stressPermuted[2]*stressPermuted[3]-stressPermuted[4]*stressPermuted[5]), 2*(stressPermuted[0]*stressPermuted[4] - stressPermuted[3]*stressPermuted[5]),
                2*(stressPermuted[1]*stressPermuted[5] - stressPermuted[3]*stressPermuted[4]);

        Dp_Ds = -2*a/9*Da_Ds + 1/3*Db_Ds;
        Dq_Ds = (std::pow(a,2)/9 - b/6)*Da_Ds - a/6*Db_Ds + 1/2*Dc_Ds;

        DP_Ds = 1/(2*std::sqrt(-p))*Dp_Ds;
        DP_Ds = (q < 0 ? (-1)*DP_Ds : DP_Ds);

        //D < 0 ?
        if (lessThan(D, 0.))
        {
            Dbeta_Ds = 1/(3*std::pow(P,3)*std::sqrt(1-std::pow(q,2)/(std::pow(P,6))))*Dq_Ds + q/(std::pow(P,4)*std::sqrt(1-std::pow(q,2)/(std::pow(P,6))))*DP_Ds;

            Dsig_1_Ds = 2*P*std::sin(beta)*Dbeta_Ds - 2*std::cos(beta)*DP_Ds - 1/3*Da_Ds;
            Dsig_2_Ds = -2*P*std::sin(beta+pi/3)*Dbeta_Ds + 2*std::cos(beta+pi/3)*DP_Ds - 1/3*Da_Ds;
            Dsig_3_Ds = -2*P*std::sin(beta-pi/3)*Dbeta_Ds + 2*std::cos(beta-pi/3)*DP_Ds - 1/3*Da_Ds;

            Df_Ds = 1/(2*f)*(Dsig_1_Ds + Dsig_2_Ds + Dsig_3_Ds);
        }
        else    //only further possibility is D = 0
        {
            //just one positive eigenvalue
            if (eigenValues[0] > 0 && eigenValues[1] == 0)
            {
                Dsig_1_Dq = -2/(9*std::pow(P,2));
                Dsig_1_DP = -4/3;
                Dsig_1_Da = -1/3;

                Dsig_1_Ds = Dsig_1_Dq*Dq_Ds + Dsig_1_DP*DP_Ds + Dsig_1_Da*Da_Ds;

                Df_Ds = 1/(2*f)*Dsig_1_Ds;
            }
            //two positive, identical eigenvalues
            else if (eigenValues[0] == 0 && eigenValues[1] > 0)
            {
                Df_Dq = -std::sqrt(2)/9*(6*P + a)/(std::abs(3*P - a)*std::pow(P,2));
                Df_DP = std::sqrt(2)/3*(15*P-2*a)/(std::abs(3*P - a));
                Df_Da = std::sqrt(2)/(3);
                Df_Da = ((3*P - a) < 0 ? (-1)*Df_Da : Df_Da);

                Df_Ds = Df_Dq*Dq_Ds + Df_DP*DP_Ds + Df_Da*Da_Ds;
            }
            //three positive eigenvlaues
            else if (eigenValues[0] > 0 && eigenValues[1] > 0)
            {
                Df_Da = a/f;
                Df_Db = -1/f;

                Df_Ds = Df_Da*Da_Ds + Df_Db*Db_Ds;
            }
        }
    }


    template<typename T>
    static bool equals(T num_1, T num_2)
    {
        return (std::abs(num_1 - num_2) < eps);
    }

    template<typename T>
    static bool greaterThan(T num_1, T num_2)
    {
        return (std::abs(num_1 - eps) > num_2);
    }

    template<typename T>
    static bool lessThan(T num_1, T num_2)
    {
        return (greaterThan(num_2, num_1));
    }

    template<typename T>
    static bool greaterOrEqual(T num_1, T num_2)
    {
        return (equals(num_1, num_2) || greaterThan(num_1, num_2));
    }

    template<typename T>
    static bool lessOrEqual(T num_1, T num_2)
    {
        return (equals(num_1, num_2) || lessThan(num_1, num_2));
    }
};
