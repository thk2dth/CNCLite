#ifndef BEZIER_H
#define BEZIER_H
#include "Path.h"

namespace CNCLite {
class Bezier : public Path
{
public:
    Bezier();
    Bezier(const std::vector<Eigen::VectorXd> &ctrlPts, double bgnPara = 0.0, double endPara = 1.0);
    ~Bezier();

public:
    void initializeFromPoints(const std::vector<Eigen::VectorXd> &ctrlPts,
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

    /// Additional member functions.
    unsigned int getDegree() const;
    double bernsteinPolynomial(unsigned int degree, unsigned int index, double u) const;
    unsigned int factorialInterval(unsigned int end, unsigned int start = 0) const;

protected:
    /// Additional member variables of bezier curve.
    unsigned int degree; // degree of the bezier curve. degree = order -1.
};
} /// namespace CNCLite

#endif // BEZIER_H
