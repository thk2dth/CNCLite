#include "BsplineLinearBlend.h"

using namespace CNCLite;
using namespace Eigen;

BsplineLinearBlend::BsplineLinearBlend()
{
    curveType = BSPLINELINEARBLEND;
    isInitialized = false;
}

BsplineLinearBlend::BsplineLinearBlend(const std::vector<VectorXd> &ctrlPts)
{
    initializeFromPoints(ctrlPts, 0.0, 1.0);
}


BsplineLinearBlend::~BsplineLinearBlend()
{

}

void BsplineLinearBlend::initializeFromPoints(const std::vector<VectorXd> &ctrlPts,
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
        ce = 0.1; // in mm.
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

void BsplineLinearBlend::initializeFromLines(const Line &line0, const Line &line1,
                                             double mr, double cr, double ce,
                                             double bgnPara, double endPara)
{
    /// Check whether the end point of the first line is identical to the start point of
    /// the second line.
    if ((line0.getEndPoint() - line1.getBeginPoint() ).norm() > EPS_CAD)
    {
        std::cout << "The end points of the first line is not identical to the begin point "
                     "of the second line." << std::endl;
    }
    initializeFromPoints(line0.getBeginPoint(), line0.getEndPoint(), line1.getEndPoint(),
                         mr, cr, ce, bgnPara, endPara);
}

void BsplineLinearBlend::initializeFromPoints(const VectorXd &p0, const VectorXd &p1,
                                              const VectorXd &p2, double mr, double cr,
                                              double ce, double bgnPara, double endPara)
{
    mRatio = mr;
    cRatio = cr;
    blendError = ce;
    beginParameter = bgnPara;
    endParameter = endPara;
    /// By default, the Bspline is of degree 3, and has five control points.
    /// Its knot vector is [0, 0, 0, 0, 0.5, 1, 1, 1, 1].
    degree = 3;
    numberControlPoints = 5;
    knotVector.resize(9);
    for (uint8_t i=0; i<4; i++)
    {
        knotVector[i] = beginParameter;
        knotVector[i+5] = endParameter;
    }
    knotVector[4] = 0.5 * (beginParameter + endParameter); // Knot vector.
    controlPoints.resize(5);
    isHomogeneous = false;
    curveType = BSPLINELINEARBLEND;
    dimension = p0.size();
    VectorXd vec1 = p1 - p0;
    VectorXd vec2 = p2 - p1;
    double len1 = vec1.norm();
    double len2 = vec2.norm();
    vec1 = vec1 / len1; // unit vector along the first line.
    vec2 = vec2 / len2; // unit vector along the second line.
    /// Maximum transition length determined by the transition error.
    double d2_blend_up = 0.0;
    /// Maximum transition length determined by the two line lengths.
    double d2_geo_up = std::min(len1, len2) * (1 - mr) * 0.5; // Eq. (9).
    double beta = 0.5 * acos(vec1.dot(vec2) );
    d2_blend_up = 2 * blendError / sin(beta); // Eq. (9).
    d2 = std::min(d2_blend_up, d2_geo_up);
    controlPoints[2] = p1;
    controlPoints[1] = controlPoints[2] - vec1 * d2;
    controlPoints[0] = controlPoints[1] - vec1 * (cRatio * d2);
    controlPoints[3] = controlPoints[2] + vec2 * d2;
    controlPoints[4] = controlPoints[3] + vec2 * (cRatio * d2);
    isInitialized = true;
    initialize();
}


double BsplineLinearBlend::calculateMaximumCurvature() const
{
    return 4*sin(beta) / (3*d2*cos(beta)*cos(beta) ); // Eq. (4).
}

double BsplineLinearBlend::calculateTransitionLength() const
{
    return (1 + cRatio) * d2; // Transition length.
}

