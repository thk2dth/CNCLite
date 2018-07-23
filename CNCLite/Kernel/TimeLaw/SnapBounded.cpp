#include "SnapBounded.h"

using namespace CNCLite;
using namespace Eigen;

SnapBounded::SnapBounded()
{
    type = POLYNOMIAL;
}

SnapBounded::~SnapBounded()
{

}


double SnapBounded::displacementAcc(double t) const
{
    double sm = kinematicConstraint(3); // alias of maximum snap.
    double vs = isAcc ? beginVelocity : endVelocity;
    double t1 = timeParameter(0);
    double t2 = timeParameter(1);
    double t3 = timeParameter(2);
    double t_bgn = 0.0;
    double t_end = 0.0;
    /// Out of the time range.
    if (t < 0.0)
        return 0.0;
    t_bgn = t_end; t_end = t1;
    if (t >= t_bgn && t < t_end) // 0 <= t < t1
        return vs*t + sm * t*t*t*t / 24.0;
    t_bgn = t_end; t_end = t_bgn + t2;
    if (t >= t_bgn && t < t_end) // t1 <= t < t1 + t2
        return vs*t +sm*(t*t*t*t1 / 6.0 - t*t*t1*t1*0.25 + t*t1*t1*t1 / 6.0 - t1*t1*t1*t1/24.0);
    t_bgn = t_end; t_end = t_bgn + t1;
    if (t >= t_bgn && t < t_end) // t1 + t2 <= t < 2*t1 + t2
        return vs*t +sm*(-t*t*t*t/24.0 + (t1/3.0+t2/6.0)*t*t*t +(-t1*t1*0.5-t1*t2*0.5-t2*t2*0.25)*t*t + t*t1*t1*t1/3.0
                         + t*t1*t1*t2*0.5 + t*t1*t2*t2*0.5 + t*t2*t2*t2/6.0 - t1*t1*t1*t1/12.0
                         - t1*t1*t1*t2/6.0 - t1*t1*t2*t2*0.25 - t1*t2*t2*t2/6.0 - t2*t2*t2*t2/24.0);
    t_bgn = t_end; t_end = t_bgn + t3;
    if (t >= t_bgn && t < t_end) // 2*t1 + t2 <= t < 2*t1 + t2 + t3
        return vs*t +sm*(t1*(t1+t2)*t*t*0.5 -t*t1*t1*t1 - 3*t*t1*t1*t2*0.5 - t*t1*t2*t2*0.5 + 7*t1*t1*t1*t1/12.0
                         + 7*t1*t1*t1*t2/6.0 + 3*t1*t1*t2*t2*0.25 + t1*t2*t2*t2/6.0);
    t_bgn = t_end; t_end = t_bgn + t1;
    if (t >= t_bgn && t < t_end) // 2*t1 + t2 + t3 <= t < 3*t1 + t2 + t3
        return vs*t +sm * (-t*t*t*t/24.0 + (t1/3.0+t2/6.0+t3/6.0)*t*t*t + (-t1*t1*0.5-t1*t2*0.5-t1*t3
                                                                           - t2*t2*0.25-t2*t3*0.5-t3*t3*0.25)*t*t + (t1*t1*t1/3.0+t1*t1*t2*0.5+2*t1*t1*t3
                                                                                                                     + t1*t2*t2*0.5+2*t1*t2*t3+t1*t3*t3+t2*t2*t2/6.0+t2*t2*t3*0.5+t2*t3*t3*0.5+t3*t3*t3/6.0)*t
                           - t1*t1*t1*t1/12.0 - t1*t1*t1*t2/6.0 - 4*t1*t1*t1*t3/3.0 - t1*t1*t2*t2*0.25 - 2*t1*t1*t2*t3
                           - t1*t1*t3*t3 - t1*t2*t2*t2/6.0 - t1*t2*t2*t3 - t1*t2*t3*t3 - t1*t3*t3*t3/3.0 - t2*t2*t2*t2/24.0
                           - t2*t2*t2*t3/6.0 - t2*t2*t3*t3*0.25 - t2*t3*t3*t3/6.0 - t3*t3*t3*t3/24.0);
    t_bgn = t_end; t_end = t_bgn + t2;
    if (t >= t_bgn && t < t_end) // 3*t1 + t2 + t3 <= t < 3*t1 + 2*t2 + t3
        return vs*t +sm*(-t1*t*t*t/6.0 + (t1*t2+t1*t3*0.5+7*t1*t1*0.25)*t*t + (-25*t1*t1*t1/6.0-4*t1*t1*t2-5*t3*t1*t1*0.5
                                                                               -t1*t2*t2-t3*t1*t2-t3*t3*t1*0.5)*t + 79*t1*t1*t1*t1/24.0 + 13*t1*t1*t1*t2/3.0 + 19*t1*t1*t1*t3/6.0
                         + 2*t1*t1*t2*t2 + 5*t1*t1*t2*t3*0.5 + 5*t1*t1*t3*t3*0.25 + t1*t2*t2*t2/3.0 + t1*t2*t2*t3*0.5
                         + t1*t2*t3*t3*0.5 + t1*t3*t3*t3/6.0);
    t_bgn = t_end; t_end = t_bgn + t1;
    if (t >= t_bgn && t < t_end) // 3*t1 + 2*t2 + t3 <= t < 4*t1 + 2*t2 + t3
        return vs*t +sm*(t*t*t*t/24.0 - (2*t1/3.0+t2/3.0+t3/6.0)*t*t*t + (4*t1*t1+4*t1*t2+2*t1*t3+t2*t2+t2*t3
                                                                          +t3*t3*0.25)*t*t + (-26*t1*t1*t1/3.0-13*t1*t1*t2-7*t1*t1*t3-7*t1*t2*t2
                                                                                              -7*t1*t2*t3-2*t1*t3*t3-4*t2*t2*t2/3.0-2*t2*t2*t3-t2*t3*t3-t3*t3*t3/6.0)*t + 20*t1*t1*t1*t1/3.0
                         + 40*t1*t1*t1*t2/3.0 + 23*t1*t1*t1*t3/3.0 + 11*t1*t1*t2*t2 + 23*t1*t1*t2*t3*0.5 + 7*t1*t1*t3*t3*0.5
                         + 13*t1*t2*t2*t2/3.0 + 13*t1*t2*t2*t3*0.5 + 7*t1*t2*t3*t3*0.5 + 2*t1*t3*t3*t3/3.0 + 2*t2*t2*t2*t2/3.0
                         + 4*t2*t2*t2*t3/3.0 + t2*t2*t3*t3 + t2*t3*t3*t3/3.0 + t3*t3*t3*t3/24.0);
    /// Out of the time range.
    if (t >= t_end) // t > 4*t1 + 2*t2 + t3
        return vs*t +sm*(4*t1*t1*t1*t1 + 8*t1*t1*t1*t2 + 3*t1*t1*t1*t3 + 5*t1*t1*t2*t2 + 4.5*t1*t1*t2*t3
                         + t1*t1*t3*t3*0.5 + t1*t2*t2*t2 + 1.5*t1*t2*t2*t3 + t1*t2*t3*t3*0.5);
    return 0.0; // for safety.
}

double SnapBounded::velocityAcc(double t) const
{
    double sm = kinematicConstraint(3); // alias of maximum snap.
    double vs = isAcc ? beginVelocity : endVelocity;
    double t1 = timeParameter(0);
    double t2 = timeParameter(1);
    double t3 = timeParameter(2);
    double t_bgn = 0.0;
    double t_end = 0.0;
    /// Out of the time range.
    if (t < 0.0)
        return vs;
    t_bgn = t_end; t_end = t1;
    if (t >= t_bgn && t < t_end) // 0 <= t < t1
        return vs + sm * t*t*t / 6.0;
    t_bgn = t_end; t_end = t_bgn + t2;
    if (t >= t_bgn && t < t_end) // t1 <= t < t1 + t2
        return vs + sm*(0.5*(t - t1)*(t - t1)*t1 + 0.5*(t - t1)*t1*t1 + t1*t1*t1 / 6.0);
    t_bgn = t_end; t_end = t_bgn + t1;
    if (t >= t_bgn && t < t_end) // t1 + t2 <= t < 2*t1 + t2
        return vs + sm*(t1*t1*t1 / 3.0 + 0.5*t1*t1*t2 - t1*t1*t + 0.5*t1*t2*t2 - t1*t2*t + t1*t*t + t2*t2*t2 / 6.0 - 0.5*t2*t2*t + 0.5*t2*t*t - t*t*t / 6.0);
    t_bgn = t_end; t_end = t_bgn + t3;
    if (t >= t_bgn && t < t_end) // 2*t1 + t2 <= t < 2*t1 + t2 + t3
        return vs + sm*(-t1*t1*t1 - 1.5*t1*t1*t2 - 0.5*t1*t2*t2 + t1*t2*t + t1*t1*t);
    t_bgn = t_end; t_end = t_bgn + t1;
    if (t >= t_bgn && t < t_end) // 2*t1 + t2 + t3 <= t < 3*t1 + t2 + t3
        return vs + sm * ((t1*t2 + t1*t1)*(t - t_bgn) + (t_bgn - t)*(t_bgn - t)*(t_bgn - t) / 6.0 + t1*t1*t1 + 1.5*t1*t1*t2 + 0.5*t1*t2*t2 + t1*t2*t3 + t1*t1*t3);
    t_bgn = t_end; t_end = t_bgn + t2;
    if (t >= t_bgn && t < t_end) // 3*t1 + t2 + t3 <= t < 3*t1 + 2*t2 + t3
        /// Note that the original equation is incorrect. -Sm*t1*t2^3 should be -Sm*t1*t2*t2.
        return vs + sm*(2 * t1*t2*t + 3.5*t1*t1*t + t1*t3*t - 0.5*t1*t*t - 4 * t1*t1*t2 - t1*t2*t2 - t1*t2*t3 - t1*t1*t1 * 25 / 6.0 - 2.5*t1*t1*t3 - 0.5*t1*t3*t3);
    t_bgn = t_end; t_end = t_bgn + t1;
    if (t >= t_bgn && t < t_end) // 3*t1 + 2*t2 + t3 <= t < 4*t1 + 2*t2 + t3
        return vs + sm*(-26 * t1*t1*t1 / 3.0 + 4 * t1*t3*t - 7 * t1*t2*t3 + 8 * t1*t2*t + 2 * t2*t3*t - 7 * t1*t2*t2 + 8 * t1*t1*t - 2 * t1*t*t - 13 * t1*t1*t2 - 7 * t1*t1*t3 - 2 * t1*t3*t3 - t2*t*t - 0.5*t3*t*t + 2 * t2*t2*t + 0.5*t3*t3*t - 2 * t2*t2*t3 - t2*t3*t3 + t*t*t / 6.0 - 4 * t2*t2*t2 / 3.0 - t3*t3*t3 / 6.0);
    /// Out of the time range.
    if (t >= t_end) // t > 4*t1 + 2*t2 + t3
        return vs + sm * t1 * (t1+t2) * (2*t1 + t2 + t3);
    return vs; // for safety.
}

double SnapBounded::accelerationAcc(double t) const
{
    double sm = kinematicConstraint(3); // alias of maximum snap.
    double t1 = timeParameter(0);
    double t2 = timeParameter(1);
    double t3 = timeParameter(2);
    double t_bgn = 0.0;
    double t_end = 0.0;
    /// Out of the time range.
    if (t < 0.0)
        return 0.0;
    t_bgn = t_end; t_end = t1;
    if (t >= t_bgn && t < t_end) // 0 <= t < t1
        return 0.5*sm*t*t;
    t_bgn = t_end; t_end = t_bgn + t2;
    if (t >= t_bgn && t < t_end) // t1 <= t < t1 + t2
        return sm * (t1*t - 0.5*t1*t1);
    t_bgn = t_end; t_end = t_bgn + t1;
    if (t >= t_bgn && t < t_end) // t1 + t2 <= t < 2*t1 + t2
        return sm*(-0.5*(t_end - t)*(t_end - t) + (t1*t2 + t1*t1));
    t_bgn = t_end; t_end = t_bgn + t3;
    if (t >= t_bgn && t < t_end) // 2*t1 + t2 <= t < 2*t1 + t2 + t3
        return sm*(t1*t2 + t1*t1);
    t_bgn = t_end; t_end = t_bgn + t1;
    if (t >= t_bgn && t < t_end) // 2*t1 + t2 + t3 <= t < 3*t1 + t2 + t3
        return sm * ((t1*t2 + t1*t1) - 0.5*(t_bgn - t)*(t_bgn - t));
    t_bgn = t_end; t_end = t_bgn + t2;
    if (t >= t_bgn && t < t_end) // 3*t1 + t2 + t3 <= t < 3*t1 + 2*t2 + t3
        return sm * (2 * t1*t2 + 3.5*t1*t1 + t1 *t3 - t1*t);
    t_bgn = t_end; t_end = t_bgn + t1;
    if (t >= t_bgn && t < t_end) // 3*t1 + 2*t2 + t3 <= t < 4*t1 + 2*t2 + t3
        return 0.5*sm*(t - t_end) *(t - t_end);
    /// Out of the time range.
    if (t >= t_end) // t > 4*t1 + 2*t2 + t3
        return 0.0;
    return 0.0; // for safety.
}

double SnapBounded::jerkAcc(double t) const
{
    double sm = kinematicConstraint(3); // alias of maximum snap.
    double t1 = timeParameter(0);
    double t2 = timeParameter(1);
    double t3 = timeParameter(2);
    double t_bgn = 0.0;
    double t_end = 0.0;
    /// Out of the time range.
    if (t < 0.0)
        return 0.0;
    t_bgn = t_end; t_end = t1;
    if (t >= t_bgn && t < t_end) // 0 <= t < t1
        return sm * t;
    t_bgn = t_end; t_end = t_bgn + t2;
    if (t >= t_bgn && t < t_end) // t1 <= t < t1 + t2
        return sm * t1;
    t_bgn = t_end; t_end = t_bgn + t1;
    if (t >= t_bgn && t < t_end) // t1 + t2 <= t < 2*t1 + t2
        return sm*(t_end - t);
    t_bgn = t_end; t_end = t_bgn + t3;
    if (t >= t_bgn && t < t_end) // 2*t1 + t2 <= t < 2*t1 + t2 + t3
        return 0.0;
    t_bgn = t_end; t_end = t_bgn + t1;
    if (t >= t_bgn && t < t_end) // 2*t1 + t2 + t3 <= t < 3*t1 + t2 + t3
        return sm * (t_bgn - t);
    t_bgn = t_end; t_end = t_bgn + t2;
    if (t >= t_bgn && t < t_end) // 3*t1 + t2 + t3 <= t < 3*t1 + 2*t2 + t3
        return -sm * t1;
    t_bgn = t_end; t_end = t_bgn + t1;
    if (t >= t_bgn && t < t_end) // 3*t1 + 2*t2 + t3 <= t < 4*t1 + 2*t2 + t3
        return sm * (t - t_end);
    /// Out of the time range.
    if (t >= t_end) // t > 4*t1 + 2*t2 + t3
        return 0.0;
    return 0.0; // for safety.
}

double SnapBounded::snapAcc(double t) const
{
    double sm = kinematicConstraint(3); // alias of maximum snap.
    double t1 = timeParameter(0);
    double t2 = timeParameter(1);
    double t3 = timeParameter(2);
    double t_bgn = 0.0;
    double t_end = 0.0;
    /// Out of the time range.
    if (t < 0.0)
        return 0.0;
    t_bgn = t_end; t_end = t1;
    if (t >= t_bgn && t < t_end) // 0 <= t < t1
        return sm;
    t_bgn = t_end; t_end = t_bgn + t2;
    if (t >= t_bgn && t < t_end) // t1 <= t < t1 + t2
        return 0.0;
    t_bgn = t_end; t_end = t_bgn + t1;
    if (t >= t_bgn && t < t_end) // t1 + t2 <= t < 2*t1 + t2
        return -sm;
    t_bgn = t_end; t_end = t_bgn + t3;
    if (t >= t_bgn && t < t_end) // 2*t1 + t2 <= t < 2*t1 + t2 + t3
        return 0.0;
    t_bgn = t_end; t_end = t_bgn + t1;
    if (t >= t_bgn && t < t_end) // 2*t1 + t2 + t3 <= t < 3*t1 + t2 + t3
        return -sm;
    t_bgn = t_end; t_end = t_bgn + t2;
    if (t >= t_bgn && t < t_end) // 3*t1 + t2 + t3 <= t < 3*t1 + 2*t2 + t3
        return 0.0;
    t_bgn = t_end; t_end = t_bgn + t1;
    if (t >= t_bgn && t < t_end) // 3*t1 + 2*t2 + t3 <= t < 4*t1 + 2*t2 + t3
        return sm;
    /// Out of the time range.
    if (t >= t_end) // t > 4*t1 + 2*t2 + t3
        return 0.0;
    return 0.0; // for safety.
}

double SnapBounded::lawVV(double vs, double ve)
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
    double t1 = 0.0, t2 = 0.0, t3 = 0.0; // temp time parameters.
    double am = kinematicConstraint(1); // alias of maximum acceleration.
    double jm = kinematicConstraint(2); // alias of maximum jerk.
    double sm = kinematicConstraint(3); // alias of maximum snap.
    double vt, pt1, pt2, pt3; // temp variables used to decrease the computation load of the algorithm.
    pt1 = am/sm;
    double dv = std::fabs(ve - vs); // velocity increase.
    if (jm*jm >= sm*am)
    {
        vt = 2*am*sqrt(pt1);
        if (dv > vt)
        {
            /// Jm cannot be reached, yet am can be reached.
            t1 = sqrt(pt1);
            t2 = 0.0;
            t3 = dv/am - 2*t1;
        }
        else
        {
            /// Jm and am cannot be reached.
            t1 = pow(0.5*dv/sm, 1.0/3);
            t2 = 0.0;
            t3 = 0.0;
        }
    }
    else
    {
        pt2 = am*am/jm;
        pt3 = 2*jm*jm*jm / (sm*sm);
        if (dv > pt2 + pt1*jm)
        {
            /// Jm and am can be reached.
            t1 = jm/sm;
            t2 = am/jm - jm/sm;
            t3 = (dv - 2 * sm*pow(t1, 3) - 3 * sm*t1*t1*t2 - sm*t1*t2*t2) / am;
        }
        else if (dv > pt3)
        {
            /// Jm can be reached, yet am cannot be reached.
            t1 = jm / sm;
            t2 = 0.5*( sqrt(t1*t1 + 4 * dv / jm ) - 3 * t1);
            t3 = 0.0;
        }
        else
        {
            /// Jm and am cannot be reached.
            t1 = pow(0.5*dv / sm, 1.0 / 3);
            t2 = 0.0;
            t3 = 0.0;
        }
    }
    VectorXd tempTime(order-1);
    tempTime << t1, t2, t3;
    timeParameter = tempTime;
    duration = 4*t1 + 2*t2 + t3;
    dispm = 0.5*(vs + ve)*duration;
    return dispm;
}

double SnapBounded::lawVD(double vs, double dis, bool acc)
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
        tempTime << 0.0, 0.0, duration;
        timeParameter = tempTime;
        if (isAcc)
            endVelocity = vs;
        else
            beginVelocity = vs;
        return vs;
    }
    double t1 = 0.0, t2 = 0.0, t3 = 0.0; // temp time parameters.
    bool exit_flag = false; // exit_flag indicates whether the quadratic and cubic equations have real roots.
    double vm = kinematicConstraint(0); // alias of maximum velocity.

    double jm = kinematicConstraint(2); // alias of maximum jerk.
    double sm = kinematicConstraint(3); // alias of maximum snap.
    if (jm*jm >= sm*am)
    {
        double d0 = 4 * am*am / sm + 4 * vs * sqrt(am / sm);
        if (dis > d0)
        {
            t1 = sqrt(am /sm);
            t2 = 0.0;
            Vector2d solution(0.0, 0.0);
            exit_flag = quadraticEquationSolver(0.5*sm*t1*t1, 3 * sm*t1*t1*t1+ vs,
                                                4 * sm*pow(t1, 4) + 4 * vs*t1 - dis, solution);
            if(exit_flag)
                t3 = findMinPosArray(solution);
            else
                return -1.0; // no solution available, return a negative velocity.
        }
        else
        {
            t2 = 0.0;
            t3 = 0.0;
            t1 = eq16Newton(vs, dis, vm, sm);
        }
    }
    else
    {
        double d1 = 4 * pow(jm, 4) / pow(sm, 3) + 4 * vs * jm / sm;
        double temp = am / jm + jm / sm;
        double d2 = am*temp*temp + 2*vs*temp;
        if(dis > d2) // dis > d2, Case 2.1.
        {
            t1 = jm / sm;
            t2 = am/jm - jm/sm;
            Vector2d solution(0.0, 0.0);
            exit_flag = quadraticEquationSolver(0.5*sm*(t1*t2 + t1*t1), sm*(3 * t1*t1*t1
                                                                            + 1.5*t1*t2*t2 + 4.5*t1*t1*t2)
                                                + vs, sm*(4 * pow(t1, 4) + 8 * t1*t1*t1*t2 + t1*t2*t2*t2
                                                          + 5 * t1*t1*t2*t2) + vs*(4 * t1 + 2 * t2) - dis, solution);
            if(exit_flag)
                t3 = findMinPosArray(solution);
            else
                return -1.0;
        }
        else if (dis > d1) // d1<dis<=d2, Case 2.2.
        {
            t1 = jm/sm;
            t3 = 0.0;
            Vector3d solution(0.0, 0.0, 0.0);
            exit_flag = cubicEquationSolver(sm*t1, sm * 5 * t1*t1, 8 * sm*t1*t1*t1 + 2 * vs,
                                            4 * sm*pow(t1, 4) + 4 * vs*t1 - dis, solution);
            if(exit_flag)
                t2 = findMinPosArray(solution);
            else
                return -1.0;
        }
        else
        {
            t2 = 0.0;
            t3 = 0.0;
            t1 = eq16Newton(vs, dis, vm, sm);
        }
    }
    VectorXd tempTime(order-1);
    tempTime << t1, t2, t3;
    timeParameter = tempTime;
    duration = 4*t1 + 2*t2 + t3;
    double tempV = vs + sm * t1 * (t1+t2) * (2*t1 + t2 + t3);
    if (isAcc)
        endVelocity = tempV;
    else
        beginVelocity = tempV;
    return tempV;
}

double SnapBounded::eq16Newton(double vs, double dm, double vm, double sm)
{
    double t0 = 0.0;
    double t1 = 0.0;
    /// If vs is zero, the initial denominator of the equation during iteration will aslo be zero.
    if (0.0 == vs)
        t1 = 2 * dm / vm;
    do
    {
        t0 = t1;
        t1 = t0 - (sm * pow(t0, 4) + vs*t0 - 0.25*dm) / (4 * sm*pow(t0, 3) + vs);
    } while (fabs(t0 - t1) > EPS_NUM);
    return t1;
}

void SnapBounded::lawCruise(double v, double dis)
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

void SnapBounded::lawCruiseDuration(double v, double t)
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
