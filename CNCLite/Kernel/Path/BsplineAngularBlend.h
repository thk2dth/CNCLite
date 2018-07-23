#ifndef BSPLINEANGULARBLEND_H
#define BSPLINEANGULARBLEND_H
#include "Path.h"
#include "Nurbs.h"

namespace CNCLite{

class BsplineAngularBlend : public Path
{
public:
    /// The spherical curve blends three points in S^2 for angular movement.
    /// First, a Bspline is constructed to blend three points in R^3.
    /// Then, the constructed Bspline is projected onto unit sphere  S^2.
    /// According to the two steps, the spherical curve is actually on a plane in R^3.
    /// The implementation is based on the algorithm in following article:
    /// "Real-time local smoothing for five-axis linear toolpath considering
    /// smoothing error constraints", J Huang, X Du, and LM Zhu, IJMTM, 2018.
    BsplineAngularBlend();
    /// Construct a spherical curve to blend the corner of three given points on S^2.
    /// The constructed Bspline is of degree 3, and has five control points. Its knot
    /// vector is always [0, 0, 0, 0, 0.5, 1, 1, 1, 1]. All the five weights are 1.
    /// In addition to the three points, a fourth point can be used to specify the
    /// parameters of the length margin ratio (i.e., lm/min(l0, l1) ), the cRatio,
    /// and the blending error. The blending error is in radius.
    /// If only three points are available, the default parameters are 1/3, 0.25,
    /// and 0.1, respectively.
    BsplineAngularBlend(const std::vector<Eigen::VectorXd > &ctrlPts);
    ~BsplineAngularBlend();

public:
    void initializeFromPoints(const std::vector<Eigen::VectorXd> &ctrlPts,
                              double bgnPara = 0.0, double endPara = 1.0);
    /// ce is in degree.
    void initializeFromPoints(const Eigen::Vector3d &p0, const Eigen::Vector3d &p1,
                              const Eigen::Vector3d &p2, double mr, double cr, double ce,
                              double bgnPara = 0.0, double endPara = 1.0);
    Eigen::VectorXd calculatePoint(double u) const;
    Eigen::VectorXd calculateDer1(double u) const;
    Eigen::VectorXd calculateDer2(double u) const;
    void calculateDerN(double u, unsigned int order, std::vector<Eigen::VectorXd> &rst) const;
    double calculateLengthInterval(double us, double ue) const;
    void scaleByParemeter(double r);
    double calculateNextParameter(double ds, double us, const INTERPOLATION_METHOD method) const;
    double calculateCurvature(double u) const;


private:
    void calculateBeginPoint();
    void calculateEndPoint();
    void calculateBeginTangent();
    void calculateEndTangent();
    void calculateTotalLength();
    void calculateBeginCurvature();
    void calculateEndCurvature();

private:
    /// The constructed Bspline.
    Nurbs nurbs;
    /// The spherical curve is on a plane in R^3. There is a norm direction of the
    /// spherical curve.
    Eigen::Vector3d normDirection;
    double d2, cRatio; // The transition length is (1+cRatio)*d2. Eq. (3).
    /// Maximum transition error, achieved at the middle of the spline.
    /// In radius for linear movement.
    double blendError;
    double beta; // Half of the angle between the two lines. In radian. Fig. 1.
    /// Length margin. lm/min(l0, l1), by default, 1/3.
    double mRatio;

public:
    /// The nurbs curve after blending.
    /// This function is used when the nurbs curve other than the spherical curve is needed.
    static Path* angularBlendingNurbs(const std::vector<Eigen::VectorXd>& pts);

};
} // namespace CNCLite

#endif // BSPLINEANGULARBLEND_H
