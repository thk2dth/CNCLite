#ifndef CNCLITECOMMON_H
#define CNCLITECOMMON_H
#include <iostream>

namespace CNCLite{
/// distance precision in R^3, in mm.
#define EPS_CAD  1.0e-6
/// angular precision, in deg.
#define EPS_ANG  1.0e-6
/// numeric precision
#define EPS_NUM  1.0e-10
#define PI  3.141592653589793 // pi
/// Since DEG2RAD and RAD2DEG are marcos, they cannot
/// be used as DEG2RAD(x+y), RAD2DEG(x+y), etc.
#define DEG2RAD(x) ( x*3.141592653589793/180.0 )
#define RAD2DEG(x) ( x*180.0/3.141592653589793 )

/// linear interpolation.
template<typename T>
inline T linearInterpolation(double u, double us, double ue,
                             const T &pts, const T &pte)
{
    if (u < us)
    {
        /// std::cerr << "The parameter is less than begin curve parameter. Return the beginning point." << std::endl;
        return pts;
    }
    if(u > ue)
    {
        /// std::cerr << "The parameter is greater than end curve parameter. Return the end point." << std::endl;
        return pte;
    }
    return pts + (pte - pts) * ((u-us)/(ue-us));
}


} /// namespace CNCLite

#endif // CNCLITECOMMON_H
