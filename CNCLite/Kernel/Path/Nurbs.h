#ifndef NURBS_H
#define NURBS_H
#include "Path.h"
#include <string>

namespace CNCLite {

class Nurbs : public Path
{
public:
    Nurbs();
    /// If homogeneous, control points are formatted as [w*x, w*y, w*z, w]. The last value of
    /// each control point is the weight. Otherwise, elements of the weights are all 1.
    /// By default, the spline is supposed to be a B-spline, ie, all weights are 1.
    Nurbs(const std::vector<Eigen::VectorXd> &ctrlPts, const std::vector<double> &knotVec,
          bool isHomo = false);
    Nurbs(const std::vector<Eigen::VectorXd> &ctrlPts, double bgnPara=0.0, double endPara=1.0);
    ~Nurbs();

public:
    /// This initialization funciton is derived from class Path, and the knot vector in this function is
    /// supposed to be in the form [u0,u0,u0, u1, u1, u1]. In this case, the spline is a Bezier curve.
    virtual void initializeFromPoints(const std::vector<Eigen::VectorXd> &ctrlPts, double bgnPara=0.0, double endPara=1.0);
    void initializeFromPoints(const std::vector<Eigen::VectorXd> &ctrlPts, const std::vector<double> &knotVec,
                              bool isHomo = false);
    /// Two additional initialization methods for Nurbs.
    /// r, scaling ratio of the curve.
    /// If initialize from control point and knot vector files, the curve is supposed to be a B-spline.
    void initializeFromFiles(const std::string &ctrlPtsFileName, const std::string &knotVecFileName,
                             bool isHomo = false, double r = 1.0);
    /// Initialize from json files.
    void initializeFromJson(const std::string &fileName, double r = 1.0, bool isValidate = false);
    /// Write Nurbs data to json file.
    bool toJsonFile(const std::string &fileName) const;

    Eigen::VectorXd calculatePoint(double u) const;
    Eigen::VectorXd calculateDer1(double u) const;
    Eigen::VectorXd calculateDer2(double u) const;
    void calculateDerN(double u, unsigned int order, std::vector<Eigen::VectorXd> &rst) const;
    double calculateLengthInterval(double us, double ue) const;
    void scaleByParemeter(double r);
    double calculateNextParameter(double ds, double us, const INTERPOLATION_METHOD method) const;
    double calculateCurvature(double u) const;

    /// Knot insertion algorithm.
    /// Refer the Nurbs Book - Algorithm 5.1.
    void knotVectorInsert(double u, uint8_t r = 1); // insert u for n times.

private:
    void calculateBeginPoint();
    void calculateEndPoint();
    void calculateBeginTangent();
    void calculateEndTangent();
    void calculateTotalLength();
    void calculateBeginCurvature();
    void calculateEndCurvature();

    /// The following functions are used to accelerate the calculation of nurbs curve.
private:
    unsigned int findSpan(double u) const;
    Eigen::MatrixXd alphaMatrix(double u, unsigned int uIndex) const;
    Eigen::VectorXd tempIterative(unsigned int index) const;
    Eigen::VectorXd tempIterativeN(unsigned int index, unsigned int order) const;

public:
    /// Additional member functions.
    unsigned int getDegree() const;
    std::vector<double> getKnotVector() const;
    bool getIsHomogeneous() const;

public:
    static unsigned int findSpan(double u, unsigned int deg, const std::vector<double> &knotVec);
    static unsigned int findMultiplicity(double u, const std::vector<double> &knotVec);
    static Eigen::VectorXd basisFunction(double u, unsigned int index, unsigned int deg,
                                         const std::vector<double> &knotVec);
    static Eigen::MatrixXd basisMatrix(const Eigen::VectorXd &u, unsigned int deg,
                                       const std::vector<double> &knotVec);
    /// Affine a N-D homogenous point to (N-1)-D point.
    /// [w*x, w*y, w*z, w] ==> [x, y, z].
    static inline Eigen::VectorXd affineHomogenous(const Eigen::VectorXd &homoCoor);

protected:
    /// Additional member variables of Nurbs curve.
    unsigned int degree; // degree of the bezier curve. degree = order -1.
    /// Use std::vector to facilitate knot insertion.
    std::vector<double> knotVector;
    bool isHomogeneous; /// Whether the control points are in homogenous coordinates.
};
} /// namespace CNCLite


#endif // NURBS_H
