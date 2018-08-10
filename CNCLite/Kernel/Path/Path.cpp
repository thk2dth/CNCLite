#include "Path.h"

using namespace CNCLite;
using namespace  Eigen;


Path::Path()
{
    isInitialized = false;
}

Path::~Path()
{

}

void Path::initialize()
{
    if(!isInitialized)
        return;
    calculateBeginPoint();
    calculateEndPoint();
    calculateBeginTangent();
    calculateEndTangent();
    calculateTotalLength();
    calculateBeginCurvature();
    calculateEndCurvature();
}

VectorXd Path::getBeginPoint() const
{
    return beginPoint;
}

VectorXd Path::getEndPoint() const
{
    return endPoint;
}

VectorXd Path::getBeginTangent() const
{
    return beginTangent;
}

VectorXd Path::getEndTangent() const
{
    return endTangent;
}

double Path::getLength() const
{
    return length;
}

double Path::getBeginParameter() const
{
    return beginParameter;
}

double Path::getEndParameter() const
{
    return endParameter;
}

bool Path::getIsInitialized() const
{
    return isInitialized;
}

double Path::getBeginCurvature() const
{
    return beginCurvature;
}

double Path::getEndCurvature() const
{
    return endCurvature;
}

Path::CURVE_TYPE Path::getCurveType() const
{
    return curveType;
}

unsigned int Path::getDimension() const
{
    return dimension;
}

unsigned int Path::getNumberControlPoints() const
{
    return numberControlPoints;
}

std::vector<VectorXd> Path::getControlPoints() const
{
    return controlPoints;
}

double Path::calculateNextParameterFirstTaylor(double ds, double us) const
{
    if(us >= endParameter)
    {
        return endParameter;
    }
    VectorXd vec = calculateDer1(us);
    double ut = us + ds / vec.norm();
    return ut >= endParameter ? endParameter : ut;
}

double Path::calculateNextParameterSecondTaylor(double ds, double us) const
{
    if(us >= endParameter)
    {
        return endParameter;
    }
    std::vector<VectorXd> rst(3, VectorXd::Zero(dimension) );
    calculateDerN(us, 2, rst);
    double len1 = rst.at(1).norm();
    if(0 == len1 && ds > 0)
    {
        /*std::cerr << "The arc increasement is non-zero, "
                     "while the first derivative is zero." << std::endl;*/
        return us;
    }
    double ut = us + ds / len1 - 0.5 * rst.at(1).dot(rst.at(2) ) * ds * ds / std::pow(len1, 4);
    if(ut < 0 )
    {
        /*std::cerr << "In second Taylor expansion method, the obtained parameter is negtive.\n "
                     "Thus the first order Taylor expansion method is utlized temporarily." << std::endl;*/
        ut = us + ds / len1;
    }
    return ut >= endParameter ? endParameter : ut;
}

double Path::calculateNextParameterNewton(double ds, double us) const
{
    if(std::fabs(ds) < EPS_NUM)
        return us;
    double u1 = calculateNextParameterSecondRuKuCompensation(ds, us);
    if( u1 >= endParameter)
        return endParameter;

    std::vector<VectorXd> vec0(2, VectorXd::Zero(dimension) );
    VectorXd vec = calculatePoint(us);
    double u0 = u1;
    calculateDerN(u0, 1, vec0);
    VectorXd vec_temp = vec0.at(0) - vec;
    double ds_temp = vec_temp.norm();
    while(std::fabs(ds - ds_temp) /ds > EPS_NUM && std::fabs(ds - ds_temp) > EPS_NUM )
    {
        if(ds_temp == 0.0)
            return u0;
        u1 = u0 - ds_temp * (ds_temp - ds) / vec_temp.dot(vec0.at(1) );
        if(u1 < beginParameter )
            u1 = beginParameter;
        if(u1 > endParameter)
            u1 = endParameter;
        u0 = u1;
        calculateDerN(u0, 1, vec0);
        vec_temp = vec0.at(0) -vec;
        ds_temp = vec_temp.norm();
    }
    return 0.5 * (u0 + u1);
}

double Path::calculateNextParameterSecondRuKuCompensation(double ds, double us) const
{
    if(us >= endParameter)
    {
        /*std::cerr<< "In second order Runge-Kutta and compensation method, the initial parameter is greater than "
                  "endparameter.\nThus the endparameter is returned." << std::endl;*/
        return endParameter;
    }
    VectorXd vec0 = calculateDer1(us);
    double k1 = 1.0 / vec0.norm();
    double u1 = us + ds * k1; /// first order expansion
    if(u1 > endParameter)
        u1 =endParameter;
    else if(u1 < beginParameter)
        u1 = beginParameter;
    VectorXd vec1 = calculateDer1(u1);
    double k2 = 1.0 / vec1.norm();
    double u_temp = us + 0.5 * ds * (k1 + k2);
    if(u_temp > endParameter)
        u_temp = endParameter;
    else if(u_temp < beginParameter)
        u_temp = beginParameter;
    VectorXd pt0 = calculatePoint(us);
    VectorXd pt_temp = calculatePoint(u_temp);
    VectorXd vec_diff = pt_temp - pt0;
    VectorXd vec_temp = calculateDer1(u_temp);
    double A = vec_temp.norm();
    double B = 2.0 * vec_temp.dot(vec_diff);
    double D = vec_diff.dot(vec_diff) - ds*ds;
    double det = B*B - 4.0*A*D;
    double du = 0.0;
    if(det >= 0)
    {
        du = 0.5 * (-B + std::sqrt(det) ) / A;
    }
    else
        du = 0.0;
    double ue = du + u_temp;
    if(ue >= endParameter)
    {
        /*std::cerr << "In second order Runge-Kutta and compensation method, the obtained parameter is greater than "
                     "endparameter.\nThus the endparameter is returned." << std::endl;*/
        return endParameter;
    }
    else if(ue <= beginParameter)
    {
        /*std::cerr << "In second order Runge-Kutta and compensation method, the obtained parameter is less than "
                     "beginParameter.\nThus the beginParameter is returned." << std::endl;*/
        return beginParameter;
    }
    else
        return ue;
}

double Path::calculateNextParameterPCI(double v, double ts, double up1,
                                      double up2, double up3, double tol,
                                      double beta, unsigned int num) const
{
    if(0.0 == v)
        return up1;
    double ui = 3* (up1 - up2) + up3; /// Initial parameters in PCI.
    if(ui == 0 && v > 0)
        return calculateNextParameterFirstTaylor(v*ts, up1);
    double vf; /// velocity fluctuation
    VectorXd prePoint = calculatePoint(up1);
    double vr, alpha;
    unsigned int ite = 0;
    vr = (calculatePoint(ui) - prePoint).norm() / ts;
    vf = std::fabs(vr - v) / v;
    while(vf > tol && std::fabs(vr-v) > tol)
    {
        alpha = v / (-beta * (vr-v) + v);
        ui = -alpha * (up1-ui) + up1;
        ui = ui >= 0 ? ui : up1;
        vr = (calculatePoint(ui) - prePoint).norm() / ts;
        vf = std::fabs(vr- v) / v;
        ite++;
        if(ite > num)
        {
            std::cerr << "Interpolation of the PCI method do not converge. Use the first-order Tarloy expansion instead." << std::endl;
            return calculateNextParameterFirstTaylor(v*ts, up1);
        }
    }
    return ui;
}


void Path::calculateLengthIntervalSimpson(double us, double ue, double *pLength, unsigned int *pCount) const
{
    double mid, len1, len2;
    mid = 0.5 * (us + ue);
    len1 = simpsonFormulation(us, ue);
    len2 = simpsonFormulation(us, mid) + simpsonFormulation(mid, ue);
    if( std::fabs(len1 -len2) > EPS_CAD )
    {
        calculateLengthIntervalSimpson(us, mid, pLength, pCount);
        calculateLengthIntervalSimpson(mid, ue, pLength, pCount);
        return; /// return is essential!
    }
    *pLength = *pLength + 0.5 * (len1 + len2);
    *pCount = *pCount + 1;
    return;
}

double Path::simpsonFormulation(double us, double ue) const
{
    double mid = 0.5 * (us + ue);
    double du = 0.0;
    if(us > ue)
    {
        du = us - ue;
    }
    else
    {
        du = ue - us;
    }
    VectorXd vec1 = calculateDer1(us);
    VectorXd vec2 = calculateDer1(mid);
    VectorXd vec3 = calculateDer1(ue);
    return du * (vec1.norm() + 4.0 * vec2.norm() + vec3.norm() ) / 6.0;
}

double Path::calculateCurvatureByDerivatives(double u) const
{
    std::vector<VectorXd> vec(3, VectorXd::Zero(3) );
    calculateDerN(u, 2, vec);
    /// If dimension != 3
    Vector3d vec1(0.0, 0.0, 0.0), vec2(0.0, 0.0, 0.0);
    for(unsigned int i=0; i<dimension; i++)
    {
        vec1(i) = vec.at(1)(i);
        vec2(i) = vec.at(2)(i);
    }
    Vector3d vecB = vec1.cross(vec2);
    double t = vec1.norm();
    if( t == 0)
    {
        std::cerr << "The tangential vector is zero." << std::endl;
        return 0.0;
    }
    return vecB.norm() / (t * t * t);
}

double Path::calculateThresholdCurvature(double vm, double am, double jm,
                                         double ts, double ce)
{
    /// Curvature determined by chord error.
    double tc = 8 * ce / (vm*vm*ts*ts + 4*ce*ce);
    /// Curvature determined by kinematic constraints.
    double tv = am / (vm * vm);
    double t = std::fmin(tc, tv);
    tv = sqrt(jm / pow(vm, 3.0) );
    t = std::fmin(t, tv);
    return t;
}

double Path::calculateTraverseVelocityR3(double cur, double vm, double am,
                                                     double jm, double ce, double ts)
{
    if (cur <= 0)
        return 0.0;
    double v1, v2; // velocities determined by chord error and kinematic constraints.
    double rou = 1 / cur;
    if (rou > 0.5 * ce)
        v1 = (2.0 / ts) * sqrt(2*rou*ce - ce*ce); // velocity determined by chord error.
    else
        v1 = 0.0;
    if (v1 > vm)
        v1 = vm;
    v2 = sqrt(am * rou); // velocity determined by maximum acceleration.
    if (v1 > v2)
        v1 = v2;
    v2 = pow(jm * rou * rou, 1.0 / 3); // velocity determined by maximum jerk.
    if (v1 > v2)
        v1 = v2;
    /*********************************************************************
     * Add the code for snap constraint here.
     * ******************************************************************/
    return v1;
}

double Path::calculateTraverseVelocityS2(double cur, double wm,
                                                      double am, double bm,
                                                      double ce, double ts)
{
    if (cur <= 0)
        return 0.0;
    if (cur <= 1)
        return wm; /// If the curvature is less than 1, return the maximum velocity.
    double w1, w2; // velocities determined by chord error and kinematic constraints.
    double cur_geo = sqrt(cur*cur - 1); // geodeesic curvature.
    /// Angular velocity determined by chord error.
    if (cur_geo == 0)
        w1 = wm;
    else
    {
        w1 = (2.0/ts) * acos(cur_geo / (cur_geo*cos(ce) + sin(ce) ) );
        if (w1 > wm)
            w1 = wm;
    }
    w2 = sqrt(am / cur_geo); // Angular velocity determined by maximum acceleration.
    if (w1 > w2)
        w1 = w2;
    w2 = pow(bm/(cur_geo*cur_geo), 1.0/3); // Angular velocity determined by maximum jerk.
    if (w1 > w2)
        w1 = w2;
    /*********************************************************************
     * Add the code for snap constraint here.
     * ******************************************************************/
    return w1;
}
