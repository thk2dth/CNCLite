#include "Arc.h"

using namespace CNCLite;
using namespace Eigen;

Arc::Arc()
{
    curveType = ARC;
    isInitialized = false;
}

Arc::Arc(const VectorXd &bgnPt, const VectorXd &endPt,
         const Eigen::VectorXd &ctPt, double bgnPara, double endPara)
{
    initializeFromPoints(bgnPt, endPt, ctPt, bgnPara, endPara);
}

Arc::~Arc()
{

}

void Arc::initializeFromPoints(const std::vector<VectorXd> &ctrlPts, double bgnPara, double endPara)
{
    initializeFromPoints(ctrlPts.at(0), ctrlPts.at(1), ctrlPts.at(2), bgnPara, endPara);
}


void Arc::initializeFromPoints(const Vector3d &bgnPt, const Vector3d &endPt, const Vector3d &ctPt,
                               double bgnPara, double endPara, bool isMinorArc)
{
    assert( 3 == bgnPt.size() );
    Vector3d vec0 = bgnPt - ctPt;
    Vector3d vec1 = endPt - ctPt;
    double len0 = vec0.norm();
    double len1 = vec1.norm();
    if ( fabs(len0 - len1) >= EPS_CAD )
    {
        std::cout << "It seems the input points do not form an arc." << std::endl;
    }
    radius = len0;
    normDirection = vec0.cross(vec1).normalized();
    angle = acos( vec0.dot(vec1) / (len0 * len1) );
    beginPoint = bgnPt;
    endPoint = endPt;
    centerPoint = ctPt;
    beginParameter = bgnPara;
    endParameter = endPara;
    this->isMinorArc = isMinorArc;
    if (!isMinorArc)
    {
        normDirection = -normDirection; // For major arc, change the norm direction of the plane.
        angle = 2*PI - angle;
    }
    dimension = bgnPt.size();
    controlPoints.push_back(beginPoint);
    controlPoints.push_back(endPoint);
    controlPoints.push_back(centerPoint);
    numberControlPoints = 3;
    curveType = ARC;
    isInitialized = true;
    initialize();
}

VectorXd Arc::calculatePoint(double u) const
{
    Vector3d vec0 = (beginPoint - centerPoint) / radius;
    // Vector3d vec1 = (endPoint - centerPoint) / radius;
    double t = angle * u;
    Vector3d temp = vec0 * cos(t) + normDirection.cross(vec0) * sin(t);
    return temp;
}

VectorXd Arc::calculateDer1(double u) const
{
    Vector3d pt = calculatePoint(u);
    Vector3d vec = pt - centerPoint; // norm of vec equals to the radius.
    Vector3d temp = normDirection.cross(vec); // norm of temp equals to the radius.
     // Norm of the first order derivative equals to angle*r/du.
    return temp * (angle / (endParameter-beginParameter) );
}

VectorXd Arc::calculateDer2(double u) const
{
    Vector3d pt = calculatePoint(u);
    Vector3d vec = pt - centerPoint; // norm of vec equals to the radius.
     // norm of the second derivative equals to r*du/angle.
    return vec * ( (beginParameter - endParameter)/angle );
}

void Arc::calculateDerN(double u, unsigned int order, std::vector<VectorXd> &rst) const
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

double Arc::calculateLengthInterval(double us, double ue) const
{
    return length * (ue - us) / (endParameter - beginParameter);
}

void Arc::scaleByParemeter(double r)
{
    endParameter = r * (endParameter - beginParameter);
    beginTangent /= r;
    endTangent /= r;
}

double Arc::calculateNextParameter(double ds, double us, const INTERPOLATION_METHOD /*method*/) const
{
    return us + (endParameter - beginParameter) * (ds / length);
}

double Arc::calculateCurvature(double /*u*/) const
{
    return 1.0 / radius;
}

void Arc::calculateBeginPoint()
{
    return;
}

void Arc::calculateEndPoint()
{
    return;
}

void Arc::calculateBeginTangent()
{
    beginTangent = calculateDer1(0.0);
}

void Arc::calculateEndTangent()
{
    endTangent = calculateDer1(1.0);
}

void Arc::calculateBeginCurvature()
{
    beginCurvature = 1.0 / radius;
}

void Arc::calculateEndCurvature()
{
    endCurvature = 1.0 / radius;
}

void Arc::calculateTotalLength()
{
    length = angle * radius;
}

double Arc::calculateCornerVelocityS2(const Vector3d &O0, const Vector3d &O1,
                                     const Vector3d &O2, double wm, double am,
                                     double ts, double ae)
{
    Vector3d n0 = O1.cross(O0);
    Vector3d n1 = O2.cross(O1);
    /// phi is the dihedral angle of planes O_0OO_1 and O_1OO_2.
    double cos_phi = n0.dot(n1) / (n0.norm() * n1.norm() );
    if(cos_phi >= 1)
        cos_phi = 1;
    else if(cos_phi <= -1)
        cos_phi = -1;
    if(1 == cos_phi)
        return wm; /// Ot0, Ot1, Ot2 is along a great arc.
    else
    {
        double phi_half = acos(cos_phi) * 0.5;
        double k = 1.0 / (sin(ae) * pow( tan(phi_half), 2.0) );
        double tan_alpha = ( cos(ae) + 1.0 / sin(phi_half) ) / (k + sin(ae));
        double wc = sqrt( am * tan_alpha ); /// in rad.
        double wt = 0.5*am*ts / sin(phi_half);
        return fmin(wm, fmin(wc, wt) );
    }
}

double Arc::getRadius() const
{
    return radius;
}

Eigen::VectorXd Arc::getCenterPoint() const
{
    return centerPoint;
}
