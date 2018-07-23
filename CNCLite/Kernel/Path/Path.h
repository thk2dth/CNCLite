#ifndef PATH_H
#define PATH_H
#include "../../ThirdParty/Eigen/Dense"
#include "../../Common/Common.hpp"
#include <cmath>
#include <cassert>
#include <vector>

namespace CNCLite {

class Path
{
public:
    Path();
    virtual ~Path();

public:
    /// Curve types.
    enum CURVE_TYPE {UNDEFINED, LINE, ARC, BEZIER, NURBS, BSPLINELINEARBLEND,
                    BSPLINEANGULARBLEND};
    /// Interpolation methods.
    /// FIRST order Taylor Expansion, SECOND order Taylor Expansion, PCI, NEWTON,
    /// SECOND RUnge-KUnta with COMpensation, respectively.
    enum INTERPOLATION_METHOD {FIRSTTE, SECONDTE, PCI, NEWTON, SECONDRUKUCOM};

public:
    Eigen::VectorXd getBeginPoint() const;
    Eigen::VectorXd getEndPoint() const;
    Eigen::VectorXd getBeginTangent() const;
    Eigen::VectorXd getEndTangent() const;
    double getLength() const;
    double getBeginParameter() const;
    double getEndParameter() const;
    bool getIsInitialized() const;
    double getBeginCurvature() const;
    double getEndCurvature() const;
    Path::CURVE_TYPE getCurveType() const;
    unsigned int getDimension() const;
    unsigned int getNumberControlPoints() const;
    std::vector<Eigen::VectorXd> getControlPoints() const;

public:
    virtual void initialize();
    virtual void initializeFromPoints(const std::vector<Eigen::VectorXd> &ctrlPts,
                                      double bgnPara = 0.0, double endPara = 1.0) = 0;

    virtual Eigen::VectorXd calculatePoint(double u) const = 0;
    virtual Eigen::VectorXd calculateDer1(double u) const = 0;
    virtual Eigen::VectorXd calculateDer2(double u) const = 0;
    /// Calculate the derivatives of the curve at u up to N order.
    /// The N+1 results are stored in rst.
    virtual void calculateDerN(double u, unsigned int order,
                               std::vector<Eigen::VectorXd> &rst) const = 0;
    virtual double calculateLengthInterval(double us, double ue) const = 0;

    /// Add virtual function for reparameterize curve parameter by scaling.
    /// The new parameter v=r*(u-us)+us, where u is the original curve parameter and
    /// r is the scaling factor.
    virtual void scaleByParemeter(double r) = 0;

    virtual double calculateNextParameter(double ds, double us,
                                         const Path::INTERPOLATION_METHOD method = Path::FIRSTTE) const = 0;
    virtual double calculateCurvature(double u) const = 0;

public:
    double calculateNextParameterFirstTaylor(double ds, double us) const;
    double calculateNextParameterSecondTaylor(double ds, double us) const;
    double calculateNextParameterNewton(double ds, double us) const;
    /// This method is a combination of second order Runge-Kuntta and parameter compensation method.
    /// It is based on the work of Jia et al., IJMTM 2016.
    double calculateNextParameterSecondRuKuCompensation(double ds, double us) const;
    /// PCI (Predictor-Corrector Interpolator) is based on Tsai and Cheng, JMSE 2003.
    /// v, velocity; ts, interpolation period; up1, previous parameter;
    /// up2, two steps ahead parameter; up3, three steps ahead parameter;
    /// tol, velocity flucation tolerance; num, maximum iteration time.
    /// beta, correction coefficient.
    double calculateNextParameterPCI(double v, double ts, double up1, double up2,
                                    double up3, double tol=0.01, double beta=0.1, unsigned int num=10) const;


protected:
    virtual void calculateBeginPoint()= 0 ;
    virtual void calculateEndPoint()= 0;
    virtual void calculateBeginTangent() = 0;
    virtual void calculateEndTangent() = 0;
    virtual void calculateTotalLength() = 0;
    virtual void calculateBeginCurvature() = 0;
    virtual void calculateEndCurvature() = 0;

protected:
    /// Calculate the length bewteen the parameter interval [us, ue] by using adaptive Simpson method.
    void calculateLengthIntervalSimpson(double us, double ue,
                                        double *pLength, unsigned int *pCount) const;
    double simpsonFormulation(double us, double ue) const;
    double calculateCurvatureByDerivatives(double u) const;


protected:
    /// Dimension of points of curve.
    unsigned int dimension;
    /// curve type
    CURVE_TYPE curveType;
    bool isInitialized;
    double beginParameter, endParameter;
    double length;
    /// Each column of the matrix is a control point.
    /// The meaning of the control points is dependent on the type of the curve.
    std::vector<Eigen::VectorXd> controlPoints;
    unsigned int numberControlPoints;
    Eigen::VectorXd beginPoint, endPoint;
    Eigen::VectorXd beginTangent, endTangent;
    double beginCurvature, endCurvature;


public:
    /// Calculate the curvature threshold according to the kinematic constraints.
    static double calculateThresholdCurvature(double vm, double am, double jm, double ts,
                                              double ce);

    /// Calculate the maximum allowable speed for the inserted transition Bspline.
    static double calculateTraverseVelocityR3(double cur, double vm, double am, double jm,
                                            double ce, double ts);
    /// Calculate the maximum allowable speed for the inserted spherical curve.
    /// wm, am, and bm are maximum angular velocity, acceleration, and jerk, respectively.
    /// They are in radius.
    /// ce is the chord error, in radius.
    /// ts is the interpolation period, in second.
    static double calculateTraverseVelocityS2(double cur, double wm, double am, double bm,
                                            double ce, double ts);
};
} /// namespace CNCLite

#endif // PATH_H
