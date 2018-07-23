#include "AccelerationBounded.h"

using namespace CNCLite;
using namespace Eigen;

AccelerationBounded::AccelerationBounded()
{
    type = POLYNOMIAL;
}

AccelerationBounded::~AccelerationBounded()
{

}

double AccelerationBounded::displacementAcc(double t) const
{
    double vs = isAcc ? beginVelocity : endVelocity;
    double a = kinematicConstraint(1); // acceleration
    double t_bgn = 0.0;
    double t_end = timeParameter(0);
    if (t<t_bgn)
        return 0.0;
    else if ( t<=t_end )
        return vs*t + 0.5*a*t*t;
    else
        return vs*t_end + 0.5*a*t_end*t_end;
}

double AccelerationBounded::velocityAcc(double t) const
{
    double vs = isAcc ? beginVelocity : endVelocity;
    double a = kinematicConstraint(1); // acceleration
    double t_bgn = 0.0;
    double t_end = timeParameter(0);
    if (t<t_bgn)
        return vs;
    else if ( t<=t_end )
        return vs + a*t;
    else
        return vs + a*t_end;
}

double AccelerationBounded::accelerationAcc(double /*t*/) const
{
    return kinematicConstraint(1);
}

double AccelerationBounded::jerkAcc(double /*t*/) const
{
    return 0.0;
}

double AccelerationBounded::snapAcc(double /*t*/) const
{
    return 0.0;
}

double AccelerationBounded::lawVV(double vs,double ve)
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
    double a = kinematicConstraint(1);
    VectorXd temp(order-1);
    temp << std::fabs(ve - vs) / a;
    timeParameter = temp;
    duration = timeParameter(0);
    dispm = 0.5 * (ve*ve-vs*vs) / a;
    return dispm;
}

double AccelerationBounded::lawVD(double vs, double dis, bool acc)
{
    type = POLYNOMIAL;
    if (acc)
        beginVelocity = vs;
    else
        endVelocity = vs;
    isAcc = acc;
    dispm = dis;
    double a = kinematicConstraint(1);
    if (a == 0.0)
    {
        /// Cruise phase.
        duration = dis / vs;
        VectorXd tempTime(order-1);
        tempTime << duration;
        timeParameter = tempTime;
        if (isAcc)
            endVelocity = vs;
        else
            beginVelocity = vs;
        return vs;
    }
    double ve = sqrt(2.0*a*dis + vs*vs);
    VectorXd temp(order-1);
    temp << (ve - vs) / a;
    timeParameter = temp;
    duration = timeParameter(0);
    if (isAcc)
        endVelocity = ve;
    else
        beginVelocity = ve;
    return ve;
}

void AccelerationBounded::lawCruise(double v, double dis)
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

void AccelerationBounded::lawCruiseDuration(double v, double t)
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

