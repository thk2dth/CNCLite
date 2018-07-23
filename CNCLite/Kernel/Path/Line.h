#ifndef LINE_H
#define LINE_H
#include "Path.h"

namespace CNCLite{

class Line : public Path
{
public:
    Line();
    /// By default, the begin and end parameters of the line are 0 and 1.
    Line(const Eigen::VectorXd &bgnPt, const Eigen::VectorXd &endPt,
         double bgnPara = 0.0, double endPara = 1.0);
    ~Line();
public:
    void initializeFromPoints(const std::vector<Eigen::VectorXd> &ctrlPts,
                              double bgnPara = 0.0, double endPara = 1.0);
    /// Additional initialization method for line.
    void initializeFromPoints(const Eigen::VectorXd &bgnPt, const Eigen::VectorXd &endPt,
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

public:
    /// Calculate corner velocity. The corner is defined by three points.
    /// The corner velocity is constrianed by the maximum acceleration and maximum velocity.
    /// \ref Luo2007AMT, Lee2011CAD.
    /// Traditional method is related to the interpolation period, which is problematic if a
    /// short interpolation period is used.
    /// The algorithm in Grbl is adopted to determine the corner velocity. This algorithm
    /// uses the centripetal acceleration and a chord error to constrain the velocity.
    /// \ref https://onehossshay.wordpress.com/2011/09/24/improving_grbl_cornering_algorithm/
    static double calculateCornerVelocity(const Eigen::Vector3d &pt0, const Eigen::Vector3d &pt1, const Eigen::Vector3d &pt2,
                                         double vm, double am,
                                          double ts = 0.001, double ce = 0.005);
};

} /// namespace CNCLite

#endif // LINE_H
