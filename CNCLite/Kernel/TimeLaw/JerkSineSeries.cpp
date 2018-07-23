#include "JerkSineSeries.h"

using namespace CNCLite;
using namespace Eigen;
JerkSineSeries::JerkSineSeries()
{
    r = 1;
    //    isIndexR = false;
    type = SINE_SERIES;
}

JerkSineSeries::~JerkSineSeries()
{

}

double JerkSineSeries::displacementAcc(double t) const
{
    std::cout << "No explict displacement function for the sine series profile" << std::endl;
    return 0.0;
}

double JerkSineSeries::velocityAcc(double t) const
{
    double vs = isAcc ? beginVelocity : endVelocity;
    if (t < 0.0)
        return vs;
    if (t >= duration)
        return vs + dvAbs;
    double th = 2 * PI * t / duration;
    double tt = 0.0; // tt \in [0, pi]
    if (th <= PI)
        tt = atan((r+1)*tan(0.5*th)/(r-1) );
    else
        tt = PI + atan((r+1)*tan(0.5*th)/(r-1) );
    return 0.5*dvAbs/PI*(th - (r-1)*(tt - 0.5*th ) ) + vs;
}

double JerkSineSeries::accelerationAcc(double t) const
{
    if (t < 0.0 || t >= duration)
        return 0.0;
    else
    {
        double th = 2 * PI * t / duration;
        return (dvAbs/duration) * (1 - (r-1)*(r*cos(th)-1)/(r*r-2*r*cos(th)+1) );
    }
}

double JerkSineSeries::jerkAcc(double t) const
{
    if (t < 0.0 || t >= duration)
        return 0.0;
    else
    {
        double th = 2 * PI * t / duration;
        double tt = r*r - 2*r*cos(th) + 1;
        return (2*dvAbs*PI*r*(r-1)*(r*r-1)*sin(th)) / (duration*duration*tt*tt);
    }
}

double JerkSineSeries::snapAcc(double t) const
{
    if (t < 0.0 || t >= duration)
        return 0.0;
    else
    {
        double th = 2 * PI * t / duration;
        double tt = r*r - 2*r*cos(th) + 1;
        return 4*dvAbs*PI*PI*r*(r+1)*(r*r-1)*(r*r*cos(th)+2*r*cos(th)*cos(th)-4*r+cos(th) ) /
                pow(duration*tt, 3.0);
    }
}

double JerkSineSeries::lawVV(double vs, double ve)
{
    type = SINE_SERIES;
    /// For trigonometric profile, only one time parameter is needed.
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
    r = optR(dvAbs);
    double t1 = tam(dvAbs, r); // duration determined by acceleration.
    double t2 = tjm(dvAbs, r); // duration determined by herk.
    duration = fmax(t1, t2);
    timeParameter << duration;
    dispm = 0.5 * (vs + ve) * duration;
    return dispm;
}

double JerkSineSeries::lawVD(double vs, double dis, bool acc)
{
    type = SINE_SERIES;
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
    /// No explicit displacement function.
    /// A binary search is needed.
    double high = kinematicConstraint(0);
    double low = vs;
    double mid = 0.0;
    do
    {
        /// lawVV() cannot be used here,
        /// because it will alter some related paramters.
        mid = 0.5 * (high + low);
        dvAbs = mid - vs;
        r = optR(dvAbs);
        duration = fmax(tam(dvAbs, r), tjm(dvAbs, r) );
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

double JerkSineSeries::optR(double dv) const
{
    double low = 1.01;
    double high = 8.0;
    if (tjm(dv, high) >= tam(dv, high) )
        return high;
    else if (tjm(dv, low) <= tam(dv, low) )
        return low;
    else
    {
        //        if (isIndexR)
        //        {
        //            /// Determine optimal r by index.
        //            uint32_t mid = 0;
        //            uint32_t low = 0;
        //            uint32_t high = mapR.cols();
        //            do
        //            {
        //                mid = uint32_t(0.5 * (high + low) );
        //                if (dv >= mapR(0, mid) )
        //                    high = mid;
        //                else
        //                    low = mid;
        //            }while(fabs(high-low) > 1);
        //            return mapR(1, mid);
        //        }
        //        else
        //        {
        /// Do binary search.
        double mid, t1, t2;
        do
        {
            mid = 0.5 * (high + low);
            t1 = tam(dv, mid);
            t2 = tjm(dv, mid);
            if (t1 > t2)
                high = mid;
            else
                low = mid;

        }while (fabs(t1 - t2) > EPS_NUM);
        return mid;
        //        }
    }
}

double JerkSineSeries::optTheta(double m) const
{
    double c = (sqrt(pow(m, 4.0)+1+34*m*m) - 1 - m*m) / (4*m);
    return acos(c); // the positive root is returned.
}

double JerkSineSeries::gm(double m) const
{
    return 2 - 2.0 / (m+1);
}

double JerkSineSeries::hm(double m) const
{
    double theta = optTheta(m); // theta when h attains the maximum value.
    double t = m*m - 2*m*cos(theta) + 1;
    return m*(m-1)*(m*m-1)*sin(theta) / (t*t);
}

double JerkSineSeries::tam(double dv, double m) const
{
    return dv * gm(m) / kinematicConstraint(1);
}

double JerkSineSeries::tjm(double dv, double m) const
{
    return sqrt(2*PI*dv*hm(m) / kinematicConstraint(2) );
}

void JerkSineSeries::lawCruise(double v, double dis)
{
    type = SINE_SERIES;
    beginVelocity = v;
    endVelocity = v;
    dispm = dis;
    isAcc = true;
    for (unsigned int i=1; i<kinematicConstraint.size(); i++)
        kinematicConstraint(i) = 0.0;
    dvAbs = 0.0;
    r = 4.5; // r is irrelevant. An arbitrary value is choosed.
    duration = dis / v;
    VectorXd tempTime;
    /// For sine and sine series profile, only one time parameter is needed.
    tempTime = VectorXd::Zero(1);
    tempTime << duration;
    timeParameter = tempTime;
}

void JerkSineSeries::lawCruiseDuration(double v, double t)
{
    type = SINE_SERIES;
    beginVelocity = v;
    endVelocity = v;
    dispm = v * t;
    isAcc = true;
    for (unsigned int i=1; i<kinematicConstraint.size(); i++)
        kinematicConstraint(i) = 0.0;
    dvAbs = 0.0;
    r = 4.5; // r is irrelevant. An arbitrary value is choosed.
    duration = t;
    VectorXd tempTime;
    /// For sine and sine series profile, only one time parameter is needed.
    tempTime = VectorXd::Zero(1);
    tempTime << duration;
    timeParameter = tempTime;
}

//void JerkSineSeries::establishMapR(uint32_t num)
//{
//    double low = 0.0; // minimum velocity increase
//    double high = kinematicConstraint(1); // maximum velocity increase
//    VectorXd dv = VectorXd::LinSpaced(num, low, high);
//    mapR = Matrix2Xd::Zero(2, num);
//    for (uint32_t i=0; i < num; i++)
//    {
//        mapR(0, i) = dv(i);
//        mapR(1, i) = optR(dv(i) );
//    }
//    isIndexR = true;
//}
