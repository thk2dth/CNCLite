#ifndef SNAPBOUNDED_H
#define SNAPBOUNDED_H

#include "TimeLaw.h"

namespace CNCLite {

/// The algorithm for snap-bounded time law is based on the following paper:
/// IJAMT2012Fan, Interpolation of parametric CNC machining path under confined jounce.
/// The description hereafter is identical to that paper.
class SnapBounded : public TimeLaw
{
public:
    SnapBounded();
    ~SnapBounded();

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
    /// Root of Eq. 16 by Newton method.
    double eq16Newton(double vs, double dm, double vm, double sm);
};

} /// namespace CNCLite

#endif // SNAPBOUNDED_H
