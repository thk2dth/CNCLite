#ifndef ACCELERATIONBOUNDED_H
#define ACCELERATIONBOUNDED_H

#include "TimeLaw.h"

namespace CNCLite {
class AccelerationBounded : public TimeLaw
{
public:
    AccelerationBounded();
    ~AccelerationBounded();

public:
    double displacementAcc(double t) const;
    double velocityAcc(double t) const;
    double accelerationAcc(double t) const;
    double jerkAcc(double t) const;
    double snapAcc(double t) const;
    double lawVV(double vs, double ve);
    double lawVD(double vs,double dis, bool acc = true);
    void lawCruise(double v, double dis);
    void lawCruiseDuration(double v, double t);
};

} /// namespace CNCLite


#endif // ACCELERATIONBOUNDED_H
