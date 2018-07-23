#include "BsplineAngularBlend.h"

using namespace CNCLite;
using namespace Eigen;

BsplineAngularBlend::BsplineAngularBlend()
{
    curveType = BSPLINEANGULARBLEND;
    isInitialized = false;
}

BsplineAngularBlend::BsplineAngularBlend(const std::vector<VectorXd> &ctrlPts)
{
    initializeFromPoints(ctrlPts, 0.0, 1.0);
}

BsplineAngularBlend::~BsplineAngularBlend()
{

}

void BsplineAngularBlend::initializeFromPoints(const std::vector<VectorXd> &ctrlPts,
                                               double bgnPara, double endPara)
{
    if (ctrlPts.size() > 4)
    {
        std::cout << "The number of points is more than 4. Use the first four points "
                     "by default." << std::endl;
    }
    double mr, cr, ce; /// Three parameters for blending.
    if (ctrlPts.size() == 3)
    {
        /// If there are only three points, use default values for the transition parameters.
        mr = 1 / 3.0;
        cr = 0.25;
        ce = 0.1; // in radius.
    }
    else
    {
        mr = ctrlPts.at(3)(0);
        cr = ctrlPts.at(3)(1);
        ce = ctrlPts.at(3)(2);
    }
    initializeFromPoints(ctrlPts.at(0), ctrlPts.at(1), ctrlPts.at(2), mr,cr, ce,
                         bgnPara, endPara);
}

void BsplineAngularBlend::initializeFromPoints(const Vector3d &p0, const Vector3d &p1,
                                               const Vector3d &p2, double mr, double cr,
                                               double ce, double bgnPara, double endPara)
{
    assert(3 == p0.size() );
    double radius = p0.norm();
    if (fabs(p1.norm() - radius) >= EPS_CAD
            || fabs(p2.norm() - radius) >= EPS_CAD)
        std::cerr << "The lengths of the three vectors are not the same." << std::endl;
    /// Construct the Bspline.
    mRatio = mr;
    cRatio = cr;
    blendError = ce; // in radius.
    beginParameter = bgnPara;
    endParameter = endPara;
    /// By default, the Bspline is of degree 3, and has five control points.
    /// Its knot vector is [0, 0, 0, 0, 0.5, 1, 1, 1, 1].

    std::vector<double> knotVector(9);
    for (uint8_t i=0; i<4; i++)
    {
        knotVector[i] = beginParameter;
        knotVector[i+5] = endParameter;
    }
    knotVector[4] = 0.5 * (beginParameter + endParameter); // Knot vector.

    Vector3d vec1 = p1 - p0;
    Vector3d vec2 = p2 - p1;
    double len1 = vec1.norm();
    double len2 = vec2.norm();
    vec1 = vec1 / len1; // unit vector along the first line.
    vec2 = vec2 / len2; // unit vector along the second line.
    Vector3d tempD = vec2 - vec1;
    double k1_ori = 0.25 * tempD.dot(p1);
    double k2_ori = tempD.squaredNorm() / 16.0;
    double cce2 = cos(ce) * cos(ce);
    double para2_ori = k1_ori * k1_ori - k2_ori * cce2; // Eq. (15).
    double d2_geo_up, d2_blend_up;
    d2_geo_up = std::min(len1, len2) * (1 - mr) * 0.5; // Eq. (9).
    if (fabs(para2_ori) <= EPS_NUM)
    {
        /// The parameter for the second order item is 0.
        if (k1_ori >= 0)
            d2 = d2_geo_up; /// The transition length is constrained by the line length.
        else
        {
            d2_blend_up = -0.5 / k1_ori;
            d2 = std::min(d2_blend_up, d2_geo_up);
        }
    }
    else if (para2_ori > 0)
        d2 = d2_geo_up; /// The transition length is constrained by the line length.
    else
    {
        double delta_ori = 4 * (1 - cce2 ) * (k2_ori - k1_ori*k1_ori) * cce2;
        /// The quadratic equation have one positive and one negative roots.
        /// The positive one is adopted.
        d2_blend_up = 0.5 * (-2 * (1-cce2) * k1_ori - sqrt(delta_ori) ) / para2_ori;
        d2 = std::min(d2_geo_up, d2_blend_up);
    }

    numberControlPoints = 5;
    VectorXd temp(3);
    for (uint8_t i=0; i<numberControlPoints; i++)
        controlPoints.push_back(temp);
    controlPoints[2] = p1;
    controlPoints[1] = controlPoints[2] - vec1 * d2;
    controlPoints[0] = controlPoints[1] - vec1 * (cRatio * d2);
    controlPoints[3] = controlPoints[2] + vec2 * d2;
    controlPoints[4] = controlPoints[3] + vec2 * (cRatio * d2);
    nurbs.initializeFromPoints(controlPoints, knotVector, false);

    /// Initialize the other parameters.
    normDirection = vec1.cross(vec2);
    normDirection.normalize();
    curveType = BSPLINEANGULARBLEND;
    dimension = p0.size();
    isInitialized = true;
    initialize();
}

VectorXd BsplineAngularBlend::calculatePoint(double u) const
{
    VectorXd p = nurbs.calculatePoint(u);
    return p.normalized();
}

VectorXd BsplineAngularBlend::calculateDer1(double u) const
{
    VectorXd v0 = nurbs.calculatePoint(u);
    VectorXd v1 = nurbs.calculateDer1(u);
    double t0 = v0.norm();
    double t01 = v0.dot(v1);
    return v1 / t0 - v0 * (t01 / pow(t0, 3.0) );
}

VectorXd BsplineAngularBlend::calculateDer2(double u) const
{
    VectorXd v0 = nurbs.calculatePoint(u);
    VectorXd v1 = nurbs.calculateDer1(u);
    VectorXd v2 = nurbs.calculateDer2(u);
    double t0 = v0.norm();
    double t1 = v1.norm();
    double t01 = v0.dot(v1);
    double t02 = v0.dot(v2);

    return v2 / t0 - v1 * (2*t01 / pow(t0, 3.0) )
            + v0 * ((3*t01*t01 - t1*t1*t0*t0 - t02*t0*t0) / pow(t0, 5.0) );
}

void BsplineAngularBlend::calculateDerN(double u, unsigned int order,
                                            std::vector<VectorXd> &rst) const
{
    VectorXd v0 = nurbs.calculatePoint(u);
    double t0 = v0.norm();
    rst[0] = v0 / t0;
    if (order >= 1)
    {
        VectorXd v1 = nurbs.calculateDer1(u);
        double t01 = v0.dot(v1);
        rst[1] = v1 / t0 - v0 * (t01 / pow(t0, 3.0) );
        if (order >= 2)
        {
            VectorXd v2 = nurbs.calculateDer2(u);
            double t1 = v1.norm();
            double t02 = v0.dot(v2);
            rst[2] = v2 / t0 - v1 * (2*t01 / pow(t0, 3.0) )
                    + v0 * ((3*t01*t01 - t1*t1*t0*t0 - t02*t0*t0) / pow(t0, 5.0) );
        }
        if (order >= 3)
        {
            std::cout << "The third and higher order derivatives are not "
                         "available yet." << std::endl;
            for (uint8_t i=3; i<=order; i++)
                rst[i] = VectorXd::Zero(dimension); // Use zero vector.
        }
    }
}

double BsplineAngularBlend::calculateLengthInterval(double us, double ue) const
{
    double len = 0.0;
    unsigned int count = 0;
    calculateLengthIntervalSimpson(us, ue, &len, &count);
    return len;
}

double BsplineAngularBlend::calculateNextParameter(double ds, double us,
                                                   const INTERPOLATION_METHOD method) const
{
    double ue = beginParameter;
    switch (method) {
    case FIRSTTE:
        ue = calculateNextParameterNewton(ds, us);
        break;
    case SECONDTE:
        ue = calculateNextParameterSecondTaylor(ds, us);
        break;
    case NEWTON:
        ue = calculateNextParameterNewton(ds, us);
        break;
    case SECONDRUKUCOM:
        ue = calculateNextParameterSecondRuKuCompensation(ds, us);
        break;
    default:
        std::cerr << "Undefined interpolation method: " << method << std::endl;
        break;
    }
    return ue;
}

double BsplineAngularBlend::calculateCurvature(double u) const
{
    return calculateCurvatureByDerivatives(u);
}

void BsplineAngularBlend::scaleByParemeter(double /*r*/)
{
    std::cerr << "Scaling method for spherical curve is not available yet." << std::endl;
    return;
}

void BsplineAngularBlend::calculateBeginPoint()
{
    beginPoint = controlPoints.at(0).normalized();
}

void BsplineAngularBlend::calculateEndPoint()
{
    endPoint = controlPoints.at(4).normalized();
}

void BsplineAngularBlend::calculateBeginTangent()
{
    beginTangent = calculateDer1(beginParameter);
}

void BsplineAngularBlend::calculateEndTangent()
{
    endTangent = calculateDer1(endParameter);
}

void BsplineAngularBlend::calculateTotalLength()
{
    length = calculateLengthInterval(beginParameter, endParameter);
}

void BsplineAngularBlend::calculateBeginCurvature()
{
    beginCurvature = calculateCurvature(beginParameter);
}

void BsplineAngularBlend::calculateEndCurvature()
{
    endCurvature = calculateCurvature(endParameter);
}

Path* BsplineAngularBlend::angularBlendingNurbs(const std::vector<VectorXd> &pts)
{
    double mr, cr, ce, d2;
    if (pts.size() == 3)
    {
        /// default value
        mr = 1 / 3.0;
        cr = 0.25;
        ce = 0.1;
    }
    else
    {
        mr = pts.at(3)(0);
        cr = pts.at(3)(1);
        ce = pts.at(3)(2);
    }

    Vector3d p0 = pts.at(0);
    Vector3d p1 = pts.at(1);
    VectorXd p2 = pts.at(2);
    Vector3d vec1 = p1 - p0;
    Vector3d vec2 = p2 - p1;
    double len1 = vec1.norm();
    double len2 = vec2.norm();
    vec1 = vec1 / len1; // unit vector along the first line.
    vec2 = vec2 / len2; // unit vector along the second line.
    Vector3d tempD = vec2 - vec1;
    double k1_ori = 0.25 * tempD.dot(p1);
    double k2_ori = tempD.squaredNorm() / 16.0;
    double cce2 = cos(ce) * cos(ce);
    double para2_ori = k1_ori * k1_ori - k2_ori * cce2; // Eq. (15).
    double d2_geo_up, d2_blend_up;
    d2_geo_up = std::min(len1, len2) * (1 - mr) * 0.5; // Eq. (9).
    if (fabs(para2_ori) <= EPS_NUM)
    {
        /// The parameter for the second order item is 0.
        if (k1_ori >= 0)
            d2 = d2_geo_up; /// The transition length is constrained by the line length.
        else
        {
            d2_blend_up = -0.5 / k1_ori;
            d2 = std::min(d2_blend_up, d2_geo_up);
        }
    }
    else if (para2_ori > 0)
        d2 = d2_geo_up; /// The transition length is constrained by the line length.
    else
    {
        double delta_ori = 4 * (1 - cce2 ) * (k2_ori - k1_ori*k1_ori) * cce2;
        /// The quadratic equation have one positive and one negative roots.
        /// The positive one is adopted.
        d2_blend_up = 0.5 * (-2 * (1-cce2) * k1_ori - sqrt(delta_ori) ) / para2_ori;
        d2 = std::min(d2_geo_up, d2_blend_up);
    }
    std::vector<double> knots(9, 0.0);
    for (uint8_t i = 5; i < 9; i++)
        knots[i] = 1.0;
    knots[4] = 0.5;

    std::vector<VectorXd> ctrlPts(5, VectorXd::Zero(3) );
    ctrlPts[2] = p1;
    ctrlPts[1] = ctrlPts[2] - vec1 * d2;
    ctrlPts[0] = ctrlPts[1] - vec1 * (cr * d2);
    ctrlPts[3] = ctrlPts[2] + vec2 * d2;
    ctrlPts[4] = ctrlPts[3] + vec2 * (cr * d2);
    Path* nrb = new Nurbs(ctrlPts, knots);
    return nrb;
}

