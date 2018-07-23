#include "JerkBounded.h"

using namespace CNCLite;
using namespace Eigen;

JerkBounded::JerkBounded()
{
    type = POLYNOMIAL;
}

JerkBounded::~JerkBounded()
{

}

double JerkBounded::displacementAcc(double t) const
{
    double jm = kinematicConstraint(2); // alias of maximum jerk.
    double vs = isAcc ? beginVelocity : endVelocity;
    double t1 = timeParameter(0);
    double t2 = timeParameter(1);
    double t_bgn = 0.0;
    double t_end = 0.0;
    /// Out of time range.
    if (t<0.0)
        return 0.0;
    t_bgn = t_end; t_end = t1;
    if (t >= t_bgn && t<t_end) // 0 <= t < t1
        return vs*t + jm*t*t*t/6.0;
    t_bgn = t_end; t_end = t_bgn + t2;
    if (t >= t_bgn && t<t_end) // t1 <= t < t1 + t2
        return vs*t + jm*(t*t*t1*0.5-t*t1*t1*0.5 + t1*t1*t1/6.0);
    t_bgn = t_end; t_end = t_bgn + t1;
    if (t >= t_bgn && t<t_end) // t1 + t2 <= t <t1 + t2 + t3
        return vs*t + jm*(-t*t*t/6.0 + (t1+t2*0.5)*t*t + (-t1*t1-t1*t2-t2*t2*0.5)*t + t1*t1*t1/3.0
                          + t1*t1*t2*0.5 + t1*t2*t2*0.5 + t2*t2*t2/6.0);
    /// Out of time range.
    if (t>=t_end)
        return vs*t + jm*(t1*t1*t1 + 3*t1*t1*t2*0.5 + t1*t2*t2*0.5);
    return 0.0; // for safety.
}

double JerkBounded::velocityAcc(double t) const
{
    double jm = kinematicConstraint(2); // alias of maximum jerk.
    double vs = isAcc ? beginVelocity : endVelocity;
    double t1 = timeParameter(0);
    double t2 = timeParameter(1);
    double t_bgn = 0.0;
    double t_end = 0.0;
    /// Out of time range.
    if (t<0.0)
        return vs;
    t_bgn = t_end; t_end = t1;
    if (t >= t_bgn && t<t_end) // 0 <= t < t1
        return vs + 0.5*jm*t*t;
    t_bgn = t_end; t_end = t_bgn + t2;
    if (t >= t_bgn && t<t_end) // t1 <= t < t1 + t2
        return  vs + jm*t1*(t-t1) + 0.5*jm*t1*t1;
    t_bgn = t_end; t_end = t_bgn + t1;
    if (t >= t_bgn && t<t_end) // t1 + t2 <= t <t1 + t2 + t3
        return  vs - 0.5*jm*(2 * t1 + t2 - t)*(2 * t1 + t2 - t) + jm*(t1*t2 + t1*t1);
    /// Out of time range.
    if (t>=t_end)
        return  vs + jm*(t1*t2+t1*t1);
    return vs; // for safety.
}

double JerkBounded::accelerationAcc(double t) const
{
    double jm = kinematicConstraint(2); // alias of maximum jerk.
    double t1 = timeParameter(0);
    double t2 = timeParameter(1);
    double t_bgn = 0.0;
    double t_end = 0.0;
    /// Out of time range.
    if (t<0.0)
        return 0.0;
    t_bgn = t_end; t_end = t1;
    if (t >= t_bgn && t<t_end) // 0 <= t < t1
        return jm*t;
    t_bgn = t_end; t_end = t_bgn + t2;
    if (t >= t_bgn && t<t_end) // t1 <= t < t1 + t2
        return jm*t1;
    t_bgn = t_end; t_end = t_bgn + t1;
    if (t >= t_bgn && t<t_end) // t1 + t2 <= t <t1 + t2 + t3
        return jm * (t_end-t);
    /// Out of time range.
    if (t>=t_end)
        return 0.0;
    return 0.0; // for safety.
}

double JerkBounded::jerkAcc(double t) const
{
    double jm = kinematicConstraint(2); // alias of maximum jerk.
    double t1 = timeParameter(0);
    double t2 = timeParameter(1);
    double t_bgn = 0.0;
    double t_end = 0.0;
    /// Out of time range.
    if (t<0.0)
        return 0.0;
    t_bgn = t_end; t_end = t1;
    if (t >= t_bgn && t<t_end) // 0 <= t < t1
        return jm;
    t_bgn = t_end; t_end = t_bgn + t2;
    if (t >= t_bgn && t<t_end) // t1 <= t < t1 + t2
        return 0.0;
    t_bgn = t_end; t_end = t_bgn + t1;
    if (t >= t_bgn && t<t_end) // t1 + t2 <= t <t1 + t2 + t3
        return -jm;
    /// Out of time range.
    if (t>=t_end)
        return 0.0;
    return 0.0; // for safety.
}

double JerkBounded::snapAcc(double /*t*/) const
{
    return 0.0;
}

double JerkBounded::lawVV(double vs, double ve)
{
    type = POLYNOMIAL;
    if (fabs(vs-ve) == 0 && kinematicConstraint(1) == 0)
    {
        timeParameter = VectorXd::Zero(order-1);
        duration = 0.0;
        dispm = 0.0;
        return 0.0; // Time law with constant velocity.
    }
    beginVelocity = vs;
    endVelocity = ve;
    isAcc = (ve >= vs);
    double t1 = 0.0, t2 = 0.0; // temp time parameters.
    double am = kinematicConstraint(1); // alias of maximum acceleration.
    double jm = kinematicConstraint(2); // alias of maximum jerk.
    double deltaV = std::fabs(ve - vs);
    double tempArcL = 0.0;
    if (deltaV <= pow(am, 2) / jm)
    {
        tempArcL = (vs + ve)*sqrt(deltaV / jm);

        t1 = sqrt(jm*deltaV) / jm;
        t2 = 0;
    }
    else if (deltaV > pow(am, 2) / jm)
    {
        tempArcL = 0.5*(vs + ve)*(am / jm + deltaV / am);

        t1 = am / jm;
        t2 = deltaV / am - am / jm;
    }
    else
    {
        tempArcL = 0.0;
        t1 = 0.0; t2 = 0.0;
    }

    VectorXd tempTime(order-1);
    tempTime << t1,  t2;
    timeParameter = tempTime;
    duration = 2*t1 + t2;
    dispm = tempArcL;
    return tempArcL;
}

double JerkBounded::lawVD(double vs, double dis, bool acc)
{
    type = POLYNOMIAL;
    if (acc)
        beginVelocity = vs;
    else
        endVelocity = vs;
    isAcc = acc;
    dispm = dis;
    // isAcc = true;
    double am = kinematicConstraint(1); // alias of maximum acceleration.
    if (am == 0.0)
    {
        /// Cruise phase.
        duration = dis / vs;
        VectorXd tempTime(order-1);
        tempTime << 0.0, duration;
        timeParameter = tempTime;
        if (isAcc)
            endVelocity = vs;
        else
            beginVelocity = vs;
        return vs;
    }
    double t1 = 0.0, t2 = 0.0; // temp time parameters.
    bool exit_flag = false; // exit_flag indicates whether the quadratic and cubic equations have real roots.

    double jm = kinematicConstraint(2); // alias of maximum jerk.
    Vector3d solution3(0.0, 0.0, 0.0);
    double tempVe = 0.0; // reachable end velocity
    exit_flag = cubicEquationSolver(1, vs, -pow(vs, 2), -pow(vs, 3) - jm*pow(dis, 2), solution3);
    if (exit_flag)
    {
        tempVe = findMinPosArray(solution3);
        if (tempVe-vs <= pow(am, 2)/jm)
        {
            /// Am can be reached.
            t1 = sqrt(jm * (tempVe-vs) ) / jm;
            t2 = 0.0;
        }
        else
        {
            Vector2d solution2(0.0, 0.0);
            exit_flag = quadraticEquationSolver(jm, pow(am, 2), -jm*pow(vs, 2)
                                                + pow(am, 2)*vs - 2 * am*jm*dis, solution2);
            if(exit_flag)
            {
                tempVe = findMinPosArray(solution2);
                t1 = am / jm;
                t2 = (tempVe-vs) / am - am / jm;
            }
        }
        VectorXd tempTime(order-1);
        tempTime << t1, t2;
        timeParameter = tempTime;
        duration = 2*t1 + t2;
        if (isAcc)
            endVelocity = tempVe;
        else
            beginVelocity = tempVe;
        return tempVe;
    }
    else
        duration = 0.0;
    return vs;
}

void JerkBounded::lawCruise(double v, double dis)
{
    type = POLYNOMIAL;
    beginVelocity = v;
    endVelocity = v;
    dispm = dis;
    isAcc = true;
    for (unsigned int i=1; i<kinematicConstraint.size(); i++)
        kinematicConstraint(i) = 0.0;
    duration = dis / v;
    VectorXd tempTime;
    tempTime =  VectorXd::Zero(order-1);
    /// Assume that the phase corresponding to last time parameter occurs only once.
    tempTime(order-2) = duration;
    timeParameter = tempTime;
}

void JerkBounded::lawCruiseDuration(double v, double t)
{
    type = POLYNOMIAL;
    beginVelocity = v;
    endVelocity = v;
    dispm = v * t;
    isAcc = true;
    for (unsigned int i=1; i<kinematicConstraint.size(); i++)
        kinematicConstraint(i) = 0.0;
    duration = t;
    VectorXd tempTime;
    tempTime =  VectorXd::Zero(order-1);
    /// Assume that the phase corresponding to last time parameter occurs only once.
    tempTime(order-2) = duration;
    timeParameter = tempTime;
}
