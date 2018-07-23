#ifndef JERKSINESERIES_H
#define JERKSINESERIES_H
#include "TimeLaw.h"

namespace CNCLite{
class JerkSineSeries : public TimeLaw
{
public:
    JerkSineSeries();
    ~JerkSineSeries();

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
    double optR(double dv) const; // Determine the optimal m.
    double optTheta(double m) const; // Value of theta when h(m,theta) obtains the maximum value.
    inline double gm(double m) const; // maximum value of g w.r.t. m.
    inline double hm(double m) const; // maximum value of h w.r.t. m.
    inline double tam(double dv, double m) const; // Duration determined by the maximum acceleration.
    inline double tjm(double dv, double m) const; // Duration determined by the maximum jerk.

private:
    double r; // Ratio of the geometric sequence.
    double dvAbs; // Absolute velocity change.
    //    /// Whether the optimal r is searched by the map or by binary search.
    //    /// False, by default.
    //    bool isIndexR;

    //public:
    //    /// Establish map for optimal R.
    //    static Eigen::Matrix2Xd establishMapR(const Eigen::VectorXd& con, uint32_t num = 1001);
};

} /// namespace CNCLite

#endif // JERKSINESERIES_H
