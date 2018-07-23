#ifndef BSPLINELINEARBLEND_H
#define BSPLINELINEARBLEND_H
#include "Nurbs.h"
#include "Line.h"

namespace CNCLite{

class BsplineLinearBlend : public Nurbs
{
public:
    /// The Bspline blends three points in R^3 for linear movement.
    /// The implementation is based on the algorithm in following article:
    /// "Real-time local smoothing for five-axis linear toolpath considering
    /// smoothing error constraints", J Huang, X Du, and LM Zhu, IJMTM, 2018.
    BsplineLinearBlend();
    /// Construct a Bspline to blend the corner of three given points.
    /// The Bspline is of degree 3, and has five control points. Its knot vector
    /// is always [0, 0, 0, 0, 0.5, 1, 1, 1, 1]. All the five weights are 1.
    /// In addition to the three points, a fourth point can be used to specify the
    /// parameters of the length margin ratio (i.e., lm/min(l0, l1) ), the cRatio,
    /// and the blending error.
    /// If only three points are available, the default parameters are 1/3, 0.25,
    /// and 0.1, respectively.
    BsplineLinearBlend(const std::vector<Eigen::VectorXd > &ctrlPts);
    ~BsplineLinearBlend();

public:
    virtual void initializeFromPoints(const std::vector<Eigen::VectorXd> &ctrlPts,
                                      double bgnPara = 0.0, double endPara = 1.0);
    void initializeFromLines(const Line &line0, const Line &line1, double mr, double cr,
                             double ce, double bgnPara = 0.0, double endPara = 1.0);
    void initializeFromPoints(const Eigen::VectorXd &p0, const Eigen::VectorXd &p1,
                              const Eigen::VectorXd &p2, double mr, double cr, double ce,
                              double bgnPara = 0.0, double endPara = 1.0);

    double calculateMaximumCurvature() const;
    double calculateTransitionLength() const;



private:
    double d2, cRatio; // The transition length is (1+cRatio)*d2. Eq. (3).
    /// Maximum transition error, achieved at the middle of the spline.
    /// In mm for linear movement.
    double blendError;
    double beta; // Half of the angle between the two lines. In radian. Fig. 1.
    /// Length margin. lm/min(l0, l1), by default, 1/3.
    double mRatio;

};
} /// namespace CNCLite

#endif // BSPLINELINEARBLEND_H
