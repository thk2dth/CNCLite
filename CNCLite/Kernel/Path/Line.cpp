#include "Line.h"

using namespace CNCLite;
using namespace Eigen;

Line::Line()
{
    curveType = LINE;
    isInitialized = false;
}

Line::Line(const Eigen::VectorXd &bgnPt, const Eigen::VectorXd &endPt, double bgnPara, double endPara)
{
    initializeFromPoints(bgnPt, endPt, bgnPara, endPara);
}

Line::~Line()
{

}

void Line::initializeFromPoints(const std::vector<Eigen::VectorXd> &ctrlPts, double bgnPara, double endPara)
{
    initializeFromPoints( ctrlPts.front(), ctrlPts.back(), bgnPara, endPara);
}

void Line::initializeFromPoints(const Eigen::VectorXd &bgnPt, const Eigen::VectorXd &endPt, double bgnPara, double endPara)
{
    beginPoint = bgnPt;
    endPoint = endPt;
    beginParameter = bgnPara;
    endParameter = endPara;
    dimension = bgnPt.size();
    controlPoints.push_back(beginPoint);
    controlPoints.push_back(endPoint);
    numberControlPoints = 2;
    curveType = LINE;
    isInitialized = true;
    initialize();
}

VectorXd Line::calculatePoint(double u) const
{
    return linearInterpolation(u, beginParameter, endParameter, beginPoint, endPoint);
}

VectorXd Line::calculateDer1(double /*u*/) const
{
    return (endPoint - beginPoint) / (endParameter - beginParameter);
}

VectorXd Line::calculateDer2(double /*u*/) const
{
    return VectorXd::Zero(dimension);
}

void Line::calculateDerN(double u, unsigned int order, std::vector<VectorXd> &rst) const
{
    rst[0] = calculatePoint(u);
    if (order >= 1)
        rst[1] = calculateDer1(u);
    if (order >= 2)
    {
        VectorXd temp = VectorXd::Zero(dimension);
        for (unsigned int i=2; i<=order; i++)
        {
            rst[i] = temp;
        }
    }
}

double Line::calculateLengthInterval(double us, double ue) const
{
    return length * (ue - us) / (endParameter - beginParameter);
}

void Line::scaleByParemeter(double r)
{
    endParameter = r * (endParameter - beginParameter);
    beginTangent /= r;
    endTangent /= r;
}

double Line::calculateNextParameter(double ds, double us, const INTERPOLATION_METHOD /*method*/) const
{
    return us + (endParameter - beginParameter) * (ds / length);
}

double Line::calculateCurvature(double /*u*/) const
{
    return 0.0;
}

void Line::calculateBeginPoint()
{
    return;
}

void Line::calculateEndPoint()
{
    return;
}

void Line::calculateBeginTangent()
{
    beginTangent = calculateDer1(0.0);
}

void Line::calculateEndTangent()
{
    endTangent = calculateDer1(1.0);
}

void Line::calculateBeginCurvature()
{
    beginCurvature = 0.0;
}

void Line::calculateEndCurvature()
{
    endCurvature = 0.0;
}

void Line::calculateTotalLength()
{
    length = (endPoint - beginPoint).norm();
}

double Line::calculateCornerVelocity(const Vector3d &pt0, const Vector3d &pt1,
                                    const Vector3d &pt2, double vm, double am,
                                     double ts, double ce)
{
    Vector3d vec0 = pt1 - pt0;
    Vector3d vec1 = pt2 - pt1;
    double len0 = vec0.norm();
    double len1 = vec1.norm();
    double cos_alpha = 0.0; /// angle between vec0 and vec1.
    /// If pt0=pt1, or pt1=pt2, let alpha = pi/2;
    if(len0 <= EPS_CAD || len1 <= EPS_CAD)
    {
        cos_alpha = 0.0;
    }
    else
    {
        cos_alpha = vec0.dot(vec1) / (len0 *len1);
        if(cos_alpha >= 1)
            cos_alpha = 1; /// Sometimes, due to nemerical issues, temp can be greater than 1.
        else if(cos_alpha <= -1)
            cos_alpha = -1; /// Also, temp can be less than -1.
    }
    if(cos_alpha == 1)
    {
        return vm; /// The three points are collinear.
    }
    else
    {
        double cos_half_alpha = sqrt( 0.5 * (1+cos_alpha) );
        double R = ce * cos_half_alpha / ( 1-cos_half_alpha );
        double vc = sqrt( am* R );
        double sin_half_alpha = sqrt(0.5*(1-cos_alpha) );
        double vt = 0.5*am*ts / sin_half_alpha;
        return fmin(vm, fmin(vc, vt) );
    }
}
