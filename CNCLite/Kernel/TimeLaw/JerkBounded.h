#ifndef JERKBOUNDED_H
#define JERKBOUNDED_H
#include "TimeLaw.h"

namespace CNCLite {
class JerkBounded : public TimeLaw
{
public:
    JerkBounded();
    ~JerkBounded();

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
};

} /// namespace CNCLite


#endif // JERKBOUNDED_H
