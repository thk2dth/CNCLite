#ifndef ARC_H
#define ARC_H
#include "Path.h"

namespace CNCLite {
class Arc : public Path
{
public:
    Arc();
    /// By default, the begin and end parameters of the arc are 0 and 1.
    /// The dimension of the arc is always 3.
    Arc(const Eigen::VectorXd &bgnPt, const Eigen::VectorXd &endPt,
        const Eigen::VectorXd &ctPt, double bgnPara = 0.0, double endPara = 1.0);
    ~Arc();
public:
    /// By default, this is a minor arc, i.e., rotation angle <= pi.
    void initializeFromPoints(const std::vector<Eigen::VectorXd> &ctrlPts,
                              double bgnPara = 0.0, double endPara = 1.0);
    /// Additional initialization method for arc. The dimension of the arc is 3.
    void initializeFromPoints(const Eigen::Vector3d &bgnPt, const Eigen::Vector3d &endPt, const Eigen::Vector3d &ctPt,
                              double bgnPara = 0.0, double endPara = 1.0, bool isMinorArc = true);
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
    double getRadius() const;
    Eigen::VectorXd getCenterPoint() const;

protected:
    /// Additional member variables of arc.
    double radius; // radius
    double angle; // rotation angle of the arc, in radius.
    Eigen::VectorXd centerPoint; // center point of the circle.
    Eigen::Vector3d normDirection; // norm direction of the arc plane.
    bool isMinorArc;

public:
    /// Rodrigues' rotation formula. Make sure that the beginning and end points are unit vectors.
    /// The point returned is also a unit vector.
    /// \ref https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
    /// \param t is the angular angle between the returned vector and v0.
    static Eigen::Vector3d sphericalInterpolation(const Eigen::Vector3d &v0, const Eigen::Vector3d &v1,
                                                  double u, bool isMinorArc = true);

    /// Calculate the corner angular velocity when blending on S^2. The method is similar to the planar situation.
    /// The three points O_i-1, O_i, and O_i+1 are unit vectors.
    /// Note that on S^2, the geodesic curvature rather than the ambient curvature is  used to constrain the velocity.
    /// \param ae, angle error, in rad, 0.005 by default;
    /// \param wm, maximum angular velocity, in rad/s;
    /// \param am, maximum angular acceleration, in rad/s^2.
    static double calculateCornerVelocityS2(const Eigen::Vector3d &O0, const Eigen::Vector3d &O1, const Eigen::Vector3d &O2,
                                           double wm, double am,
                                           double ts = 0.001, double ae = 0.005);

};
} /// namespace CNCLite


#endif // ARC_H
