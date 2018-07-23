#include "JerkSine.h"

using namespace CNCLite;
using namespace Eigen;

JerkSine::JerkSine()
{
    type = SINE;
}

JerkSine::~JerkSine()
{

}

double JerkSine::displacementAcc(double t) const
{
    if (t <= 0.0)
        return 0.0;
    if (t >= duration)
        return dispm;
    double dt = 0.25*dvAbs*duration / (PI*PI);
    double vs = isAcc ? beginVelocity : endVelocity;
    return dt*cos(2*PI*t / duration) + 0.5*dvAbs*t*t / duration + vs*t - dt;
}

double JerkSine::velocityAcc(double t) const
{
    double vs = isAcc ? beginVelocity : endVelocity;
    if (t <= 0)
        return vs;
    if (t >= duration)
        return vs + dvAbs;
    return -0.5*dvAbs*sin(2*PI*t/duration) / PI + dvAbs*t / duration + vs;
}

double JerkSine::accelerationAcc(double t) const
{
    if (t <= 0 || t >= duration)
        return 0.0;
    else
        return dvAbs/duration * (1 - cos(2*PI*t / duration) );
}

double JerkSine::jerkAcc(double t) const
{
    if (t < 0 || t > duration)
        return 0.0;
    else
    {
        double jt = 2*PI*dvAbs / (duration * duration);
        return jt * sin(2*PI*t / duration);
    }
}

double JerkSine::snapAcc(double t) const
{
    if (t < 0 || t > duration)
        return 0.0;
    else
    {
        double st = 4*PI*PI*dvAbs / (duration*duration*duration);
        return st * cos(2*PI*t / duration);
    }
}

double JerkSine::lawVV(double vs, double ve)
{
    type = SINE;
    timeParameter = VectorXd::Zero(1);
    if (fabs(vs-ve) == 0 && kinematicConstraint(1) == 0)
    {
        duration = 0.0;
        dispm = 0.0;
        return 0.0; // Time law with constant velocity.
    }
    beginVelocity = vs;
    endVelocity = ve;
    isAcc = (ve >= vs);
    dvAbs = fabs(ve - vs);
    /// Acceleration duration.
    duration = transitionTime(dvAbs, kinematicConstraint(1),
                              kinematicConstraint(2) );
    timeParameter << duration;
    dispm = 0.5 * (vs + ve) * duration;
    return dispm;
}

double JerkSine::lawVD(double vs, double dis, bool acc)
{
    type = SINE;
    if (acc)
        beginVelocity = vs;
    else
        endVelocity = vs;
    isAcc = acc;
    dispm = dis;
    timeParameter = VectorXd::Zero(1);
    if (kinematicConstraint(1) == 0.0)
    {
        /// Cruise phase.
        duration = dis / vs;
        timeParameter << duration;
        if (isAcc)
            endVelocity = vs;
        else
            beginVelocity = vs;
        return vs;
    }
    /// Though explicit displacement function exists, the end
    /// velocity cannot be determined by the displacement increase
    /// directly. A binary search is indispensable.
    double high = kinematicConstraint(0);
    double low = vs;
    double mid = 0.0;
    do
    {
        /// lawVV() cannot be used here,
        /// because it will alter some related paramters.
        mid = 0.5 * (high + low);
        dvAbs = mid - vs;
        duration = transitionTime(dvAbs, kinematicConstraint(1),
                                  kinematicConstraint(2) );
        dispm = 0.5 * (mid + vs) * duration;
        if (dispm > dis)
            high = mid;
        else
            low = mid;
    }while (fabs(dispm -dis) > EPS_CAD);
    if (isAcc)
        endVelocity = mid;
    else
        beginVelocity = mid;
    timeParameter << duration;
    return mid;
}

double JerkSine::transitionTime(double dv, double am, double jm)
{
    double temp = 0.5 * PI * am * am / jm;
    if (dv > temp)
        return 2 * dv / am; // Acceleration constraint is active.
    else
        return sqrt(2 * PI * dv / jm); // Jerk constraint is active.
}

void JerkSine::lawCruise(double v, double dis)
{
    type = SINE;
    beginVelocity = v;
    endVelocity = v;
    dispm = dis;
    isAcc = true;
    for (unsigned int i=1; i<kinematicConstraint.size(); i++)
        kinematicConstraint(i) = 0.0;
    dvAbs = 0.0;
    duration = dis / v;
    VectorXd tempTime;
    /// For sine and sine series profile, only one time parameter is needed.
    tempTime = VectorXd::Zero(1);
    tempTime << duration;
    timeParameter = tempTime;
}

void JerkSine::lawCruiseDuration(double v, double t)
{
    type = SINE;
    beginVelocity = v;
    endVelocity = v;
    dispm = v * t;
    isAcc = true;
    for (unsigned int i=1; i<kinematicConstraint.size(); i++)
        kinematicConstraint(i) = 0.0;
    dvAbs = 0.0;
    duration = t;
    VectorXd tempTime;
    /// For sine and sine series profile, only one time parameter is needed.
    tempTime = VectorXd::Zero(1);
    tempTime << duration;
    timeParameter = tempTime;
}

