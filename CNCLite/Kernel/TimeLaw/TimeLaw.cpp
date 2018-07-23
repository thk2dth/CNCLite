#include "TimeLaw.h"

using namespace Eigen;
using namespace CNCLite;

TimeLaw::TimeLaw()
{
    order = 0;
    duration = 0.0;
    dispm = 0.0;
    isAcc = true;
}

TimeLaw::~TimeLaw()
{

}

TimeLaw::TimeLaw(const VectorXd &kineConst, bool acc)
{
    initialize(kineConst, acc);
}

void TimeLaw::initialize(const VectorXd &kineConst, bool acc)
{
    kinematicConstraint = kineConst;
    order = kinematicConstraint.size();
    duration = 0.0;
    isAcc = acc;
}

double TimeLaw::reachableVelocity(double dis)
{
    return lawVD(0.0, dis);
}

double TimeLaw::displacement(double t) const
{
    if (isAcc)
        return displacementAcc(t);
    else
        return dispm - displacementAcc(duration-t);
}

double TimeLaw::velocity(double t) const
{
    if (isAcc)
        return velocityAcc(t);
    else
        return velocityAcc(duration-t);
}

double TimeLaw::acceleration(double t) const
{
    if (isAcc)
        return accelerationAcc(t);
    else
        return -accelerationAcc(duration-t);
}

double TimeLaw::jerk(double t) const
{
    if (isAcc)
        return jerkAcc(t);
    else
        return jerkAcc(duration-t);
}

double TimeLaw::snap(double t) const
{
    if (isAcc)
        return snapAcc(t);
    else
        return -snapAcc(duration-t);
}

void TimeLaw::stretchByV(double v, bool atEnd)
{
    double vs = beginVelocity;
    double ve = endVelocity;
    /// Scale the constraints by identical factor, and the time parameters remain the same.
    /// Scaling factor for kinematic constraints.
    double k = 1.0;
    if (atEnd)
    {
        endVelocity = v; // set end velocity to v.
        /// Scale the trajectory by altering the end velocity.
        if (fabs(v -ve) <= EPS_NUM || duration == 0)
            return;
        k = fabs((v-vs) / (ve - vs) );
    }
    else
    {
        beginVelocity = v; // set begin velocity to v.
        /// Scale the trajectory by altering the start velocity.
        if (fabs(v -vs) <= EPS_NUM || duration == 0)
            return;
        k = fabs((v-ve) / (ve - vs) );
    }
    /// Type of the time law could be changed.
    isAcc = (endVelocity >= beginVelocity);
    dispm = 0.5 * (endVelocity+beginVelocity) * duration;
    /// The maximum velocity does not need to be scaled.
    for (int i=1; i<kinematicConstraint.size(); i++)
        kinematicConstraint(i) *= k;
}

void TimeLaw::stretchByTime(double t)
{
    if (fabs(t-duration) <= EPS_NUM || duration == 0)
        return;
    /// The time parameters and all the kinematic constraints are changed.
    /// Start and end velocities remain the same.
    double k = t / duration;
    timeParameter *= k;
    dispm *= k;
    duration = t;
    /// The maximum velocity does not need to be scaled.
    for (int i=1; i<kinematicConstraint.size(); i++)
        kinematicConstraint(i) /= pow(k, i);
}

void TimeLaw::reset()
{
    order = 0;
    duration = 0.0;
    dispm = 0.0;
    isAcc = true;
}

double TimeLaw::getDuration() const
{
    return duration;
}

double TimeLaw::getDispm() const
{
    return dispm;
}

double TimeLaw::getBeginVelocity() const
{
    return beginVelocity;
}

double TimeLaw::getEndVelocity() const
{
    return endVelocity;
}

bool TimeLaw::getIsAcc() const
{
    return isAcc;
}

void TimeLaw::setIsAcc(bool value)
{
    isAcc = value;
}

TimeLaw::LawType TimeLaw::getType() const
{
    return type;
}

double TimeLaw::powSigned(double x, double y)
{
    if (x>=0)
        return pow(x, y);
    else
        return -pow(fabs(x), y);
}

double TimeLaw::findMinPosArray(const VectorXd &arr)
{
    double temp = 1e6; // some large positive value.
    for (int i=0; i<arr.size(); i++)
    {
        if (arr(i) > 0.0 && arr[i] < temp)
            temp = arr(i);
    }
    /// If there is no positive value in arr, return 0.0.
    return (temp==1e6  ?  0.0 : temp);
}

bool TimeLaw::cubicEquationSolver(double a, double b, double c, double d, Vector3d &solution)
{
    double RootA, RootB, RootC, RootDelta;
    double tempY1, tempY2;

    RootA = b * b - 3 * a * c;
    RootB = b * c - 9 * a * d;
    RootC = c * c - 3 * b * d;
    RootDelta = RootB * RootB - 4 * RootA * RootC;

    if (RootA == 0 && RootB == 0)
    {
        solution(0) = -c / b;
        solution(1) = -solution(0);
        solution(2) = -solution(0);
    }
    else if (RootDelta > 0)
    {
        tempY1 = RootA * b + 3 * a * (-RootB + sqrt(RootDelta)) / 2;
        tempY2 = RootA * b + 3 * a * (-RootB - sqrt(RootDelta)) / 2;
        solution(0) = (-b - (powSigned(tempY1, 1.0/3) + powSigned(tempY2, 1.0/3))) / (3 * a);
        solution(1) = 0;
        solution(2) = 0;
    }
    else if (RootDelta < 0)
    {
        if (RootA > 0)
        {
            tempY1 = 2 * (2 * RootA * b + 3 * a * RootB) / sqrt(RootA * RootA * RootA);

            if (tempY1 > -1 && tempY1 < 1)
            {
                tempY2 = acos(tempY1);
                solution(0) = (-b - 2 * sqrt(RootA) * cos(tempY2 / 3)) / (3 * a);
                solution(1) = (-b + sqrt(RootA) * (cos(tempY2 / 3) + sqrt(3.0) * sin(tempY2 / 3))) / (3 * a);
                solution(2) = (-b + sqrt(RootA) * (cos(tempY2 / 3) - sqrt(3.0) * sin(tempY2 / 3))) / (3 * a);
            }
            else
            {
                return false;
            }
        }
        else
        {
            return false;
        }
    }
    else
    {
        if (RootA != 0)
        {
            tempY1 = RootB / RootA;
            solution(0) = -b / a + tempY1;
            solution(1) = -tempY1 / 2;
            solution(2) = 0;
        }
        else
        {
            std::cerr << "No solution available for the cubic function.\n";
            return false;
        }
    }

    return true;
}

bool TimeLaw::quadraticEquationSolver(double a, double b, double c, Vector2d &solution)
{
    double delta = b*b - 4 * a*c;
    if (delta < 0)
    {
        solution(0) = 0.0;
        solution(1) = 0.0;
        return false; // no real root.
    }
    else
    {
        double temp = sqrt(delta);
        solution(0) = 0.5*(-b + temp) / a;
        solution(1) = 0.5*(-b - temp) / a;
        return true;
    }
}
