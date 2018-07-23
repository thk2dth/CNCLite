#ifndef JERKSINE_H
#define JERKSINE_H
#include "TimeLaw.h"

namespace CNCLite{
class JerkSine : public TimeLaw
{
public:
    JerkSine();
    ~JerkSine();

public:
    double displacementAcc(double t) const;
    double velocityAcc(double t) const;
    double accelerationAcc(double t) const;
    double jerkAcc(double t) const;
    double snapAcc(double t) const;
    double lawVV(double vs, double ve);
    double lawVD(double vs, double dis, bool acc = true);
    void lawCruise(double v, double dis);
    void lawCruiseDuration(double v, double t);

private:
    /// Determine the duration for accelerating.
    double transitionTime(double dv, double am, double jm);

private:
    double dvAbs; // Absolute velocity change.
};

} /// namespace CNCLite

#endif // JERKSINE_H
