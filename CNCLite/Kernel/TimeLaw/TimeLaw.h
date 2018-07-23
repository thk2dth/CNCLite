#ifndef TIMELAW_H
#define TIMELAW_H

#include "../../ThirdParty/Eigen/Dense"
#include "../../Common/Common.hpp"
#include <cmath>
#include <cassert>
#include <vector>

namespace CNCLite {
/// TimeLaw is the piecewise polynomial of an **acceleration** phase.
/// Note that TimeLaw does not consider the time sampling,
/// hence it needs no sampling period.
/// Also note that the end velocity in TimeLaw should be greater than the start velocity.
class TimeLaw
{
public:
    TimeLaw();
    /// The begin velocity is essential for the time law, yet
    /// end velocity and displacement are not essential.
    /// When tableEnable is true, the constructor will establish a table which stores the distance required s
    /// to accelerate the velocity from 0 to ve, where ve belongs to [0, vm].
    TimeLaw(const Eigen::VectorXd &kineConst, bool acc = true);
    virtual ~TimeLaw();

public:
    /// Types of the time law.
    /// Type POLYNOMIAL includes acceleration-bounded, jerk-bounded, and snap-bounded laws.
    /// Type SINE represents the jerk profile by a single sine function.
    /// Type SINE_SERIES represents the jerk profile by sine series.
    enum LawType{POLYNOMIAL, SINE, SINE_SERIES};

public:
    void initialize(const Eigen::VectorXd &kineConst, bool acc = true);
    /// Assume that the kinematic constraints are symmetric during the ACC or DEC phase.
    double displacement(double t) const;
    double velocity(double t) const;
    double acceleration(double t) const;
    double jerk(double t) const;
    double snap(double t) const;
    /// Determine the time parameters according to the end velocity or the displacement.
    /// lawVV() returns the distance required to accelerate the velocity from vs to ve.
    virtual double lawVV(double vs, double ve) = 0;
    /// lawVD() returns the velocity attained when accelerating from vs with the specified distance dis.
    virtual double lawVD(double vs, double dis, bool acc = true) = 0;
    /// lawCruise() is used to set the law as a cruise phase.
    virtual void lawCruise(double v, double dis) = 0;
    virtual void lawCruiseDuration(double v, double t) = 0;
    /// Stretch time law with specified paramters. The stretched time law is not time-optimal.
    /// The stretch method is used for time synchronization.
    /// We use two methods to stretch the time law:
    /// (a) Stretch thed boundary velocity and the kinematic constraitns; (2) Increase the duration.
    /// When stretching, the time law does not need to be re-calculated.
    /// \param t, duration to be synchronized to;
    /// \param v, velocity to be decreased to.
    /// \param atEnd, if true, the end velocity is stretched; otherwise, the start velocity velocity
    /// will be stretched. By default, the end velocity is decreased.
    /// Scale the acceleration, jerk and snap by k, k*k and k^3, respectively.
    void stretchByV(double v, bool atEnd = true);
    void stretchByTime(double t);
    /// Get the maximum reachable velocity by the specified distance.
    /// This reachable velocity is used to determine the boundary velocity.
    double reachableVelocity(double dis);
    /// Reset the time law to duration = 0, dispm = 0.
    void reset();

public:
    double getDuration() const;
    double getDispm() const;
    double getBeginVelocity() const;
    double getEndVelocity() const;
    bool getIsAcc() const;
    void setIsAcc(bool value);
    LawType getType() const;

protected:
    /// Calculate the information with analytical functions.
    virtual double snapAcc(double t) const = 0;
    virtual double jerkAcc(double t) const = 0;
    virtual double accelerationAcc(double t) const = 0;
    virtual double velocityAcc(double t) const = 0;
    virtual double displacementAcc(double t) const = 0;


protected:
    Eigen::VectorXd kinematicConstraint; // kinematic constraint, symmetric.
    /// Order of the time law, e.g., for acceleration-bounded law, order is 2.
    /// Order is also the number of kinematic constraints.
    unsigned int order;
    Eigen::VectorXd timeParameter; /// number of the time parameter is order-1.
    double beginVelocity; // start velocity
    double endVelocity; // end velocity
    double duration; /// duration of the law.
    double dispm; // displacement.
    bool isAcc; // if the law is an ACC phase.
    LawType type; // type of the trajectory.

public:
    /// If x<0, pow(x,y) may returns the wrong value.
    /// Therefore, when x<0, powSigned returns -pow(fabs(x), y).
    static double powSigned(double x, double y);
    /// Find the minimum positive value in an array.
    static double findMinPosArray(const Eigen::VectorXd &arr);
    /// Solve cubic equation.
    static bool cubicEquationSolver(double a, double b, double c, double d, Eigen::Vector3d &solution);
    /// Solve quadratic equation.
    static bool quadraticEquationSolver(double a, double b, double c, Eigen::Vector2d &solution);
};

} /// namespace CNCLite


#endif // TIMELAW_H
