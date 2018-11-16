#include "Kinematic.h"

using namespace CNCLite;
using namespace Eigen;

Kinematic::Kinematic()
{

}

Kinematic::Kinematic(const VectorXd &param, const VectorXd &lb,
                     const VectorXd &ub)
    : parameters(param), lowerBound(lb), upperBound(ub)
{
    dof = lb.size();
}

Kinematic::Kinematic(const VectorXd &param, const VectorXd &lb,
                     const VectorXd &ub, const MatrixXd &kineConstraint)
    : parameters(param), lowerBound(lb), upperBound(ub),
      kinematicConstraint(kineConstraint)
{
    dof = lb.size();
}

Kinematic::~Kinematic()
{

}

MachineTool::MachineTool()
{

}

MachineTool::MachineTool(const VectorXd &param, const VectorXd &lb,
                         const VectorXd &ub)
{
    parameters = param;
    lowerBound = lb;
    upperBound = ub;
    dof = lb.size();
}

MachineTool::MachineTool(const VectorXd &param, const VectorXd &lb,
                         const VectorXd &ub, const MatrixXd &kineConstraint)
{
    parameters = param;
    lowerBound = lb;
    upperBound = ub;
    dof = lb.size();
    kinematicConstraint = kineConstraint;
}

MachineTool::~MachineTool()
{

}

VectorXd MachineTool::forwardKinematic(const VectorXd &axialCoor,
                                       const VectorXd & /*axialCoorPrevious*/) const
{
    Vector3d trans = axialCoor.head(3);
    Vector2d rot = axialCoor.tail(2);
    Vector3d ori = orientationCoor(rot);
    Vector3d pos = positionCoor(trans, rot);
    VectorXd fkt(6);
    fkt << pos(0), pos(1), pos(2),
            ori(0), ori(1), ori(2);
    return fkt;
}

VectorXd MachineTool::inverseKinematic(const VectorXd &cartesianCoor,
                                       const VectorXd &axialCoorPrevious) const
{
    Vector3d pos = cartesianCoor.head(3);
    Vector3d ori = cartesianCoor.tail(3);
    ori.normalize();
    Vector2d rotPre = axialCoorPrevious.tail(2);
    Vector2d rot = rotaryCoor(ori, rotPre);
    Vector3d trans = translationalCoor(pos, rot);
    VectorXd ikt(5);
    ikt << trans(0), trans(1), trans(2),
            rot(0), rot(1);
    return ikt;
}

Vector2d MachineTool::rotaryCoor(const Vector3d &ori, const Vector2d &rotPre) const
{
    std::vector<Vector2d> rotSeries, rotSeriesPruned;
    // Eigen::Vector3d oriNormalized = ori.normalized();
    allRot(ori, rotSeries);
    removeRotByBound(rotSeries, rotSeriesPruned);
    Vector2d rotTemp = selectRot(rotSeriesPruned, rotPre);
    return rotTemp;
}

void MachineTool::removeRotByBound(const std::vector<Vector2d> &rotSeries,
                                   std::vector<Vector2d>& rotSeriesPruned) const
{
    for (uint8_t i = 0; i < rotSeries.size(); i++)
    {
        // Suppose the limits of the rotary axes are the last two in the lowerer and upper bounds.
        if ( rotSeries.at(i)(0) >= lowerBound(dof-2) && rotSeries.at(i)(0) <= upperBound(dof-2) &&
             rotSeries.at(i)(1) >= lowerBound(dof-1) && rotSeries.at(i)(1) <= upperBound(dof-1) )
            rotSeriesPruned.push_back(rotSeries.at(i) );
    }
}

Vector2d MachineTool::selectRot(const std::vector<Vector2d> &rotSeries,
                                const Vector2d& rotPre) const
{
    double dis;
    double disMin = std::numeric_limits<double>::max();
    Vector2d rot = Vector2d::Zero();
    for (uint16_t i =0; i < rotSeries.size(); i++)
    {
        dis = (rotSeries.at(i) - rotPre).norm();
        if (dis < disMin)
        {
            disMin = dis;
            rot = rotSeries.at(i);
        }
    }
    return rot;
}

VMC_C20::VMC_C20()
{

}

VMC_C20::VMC_C20(const VectorXd &param, const VectorXd &lb, const VectorXd &ub)
{
    parameters = param;
    lowerBound = lb;
    upperBound = ub;
    dof = lb.size();
}

VMC_C20::VMC_C20(const VectorXd &param, const VectorXd &lb, const VectorXd &ub,
                 const MatrixXd &kineConstraint)
{
    parameters = param;
    lowerBound = lb;
    upperBound = ub;
    dof = lb.size();
    kinematicConstraint = kineConstraint;
}

VMC_C20::~VMC_C20()
{

}

Vector3d VMC_C20::translationalCoor(const Vector3d &pos,
                                    const Vector2d &rot) const
{
    double A = DEG2RAD(rot(0) );
    double C = DEG2RAD(rot(1) );
    double Sa = sin(A);
    double Ca = cos(A);
    double Sc = sin(C);
    double Cc = cos(C);
    /// Suppose the first three parameters are related to
    /// the translational information.
    Vector3d delta = pos - parameters.head(3);
    Matrix3d mat;
    mat << Cc, Sc, 0,
            Ca*Sc, -Ca*Cc, -Sa,
            Sa*Sc, -Sa*Cc, Ca;
    Vector3d temp;
    temp << parameters(0), -parameters(1), parameters(2);
    return  mat * delta + temp;
}

Vector3d VMC_C20::orientationCoor(const Vector2d &rot) const
{
    double A = DEG2RAD(rot(0) );
    double C =DEG2RAD(rot(1) );
    Vector3d temp;
    temp << sin(A)*sin(C), -sin(A)*cos(C), cos(A);
    return temp;
}

Vector3d VMC_C20::positionCoor(const Vector3d &trans,
                               const Vector2d &rot) const
{
    double A = DEG2RAD(rot(0) );
    double C = DEG2RAD(rot(1) );
    double Sa = sin(A);
    double Ca = cos(A);
    double Sc = sin(C);
    double Cc = cos(C);
    Vector3d delta;
    delta << trans(0) - parameters(0),
            trans(1) + parameters(1),
            trans(2) - parameters(2);
    Matrix3d mat;
    mat << Cc, Ca*Sc, Sa*Sc,
            Sc, -Ca*Cc, -Sa*Cc,
            0, -Sa, Ca;
    return mat*delta + parameters.head(3);
}

Vector2d VMC_C20::oneRot(const Vector3d &ori) const
{
    double A, C;
    double i = ori(0);
    double j = ori(1);
    double k = ori(2);
    A = acos(k);
    assert(A>=0 && A <= PI);
    if ((fabs(i)<EPS_NUM && fabs(j)<EPS_NUM) || fabs(fabs(k)-1)<EPS_NUM)
    {
        C = 0;
    }
    else
    {
        /****C= -atan2(i, j) is not correct!**********/
        C = atan2(i, -j); // ori = [Sa*Sc, -Sa*Cc, Ca]
        if (C < 0)
            C += 2*PI; // ensure C belongs to [0, 2*PI]
    }
    return Vector2d(RAD2DEG(A), RAD2DEG(C) );
}

void VMC_C20::allRot(const Vector3d &ori,
                     std::vector<Vector2d> &rotSeries) const
{
    assert(ori.norm() > 0);
    Vector2d rot = oneRot(ori);
    double A = rot(0);
    double C = rot(1);
    Vector2d temp;
    /// All the possible solutions for the rotational coordinates are listed.
    /// (A, C)
    rotSeries.push_back(rot);
    /// (-A, C+PI)
    temp << -A, C+180;
    rotSeries.push_back(temp);
    /// (A, C+360)
    temp << A, C+360;
    rotSeries.push_back(temp);
    /// (A, C-360)
    temp << A, C-360;
    rotSeries.push_back(temp);
    /// (-A, C+540)
    temp << -A, C+540;
    rotSeries.push_back(temp);
    /// (-A, C-180)
    temp << -A, C-180;
    rotSeries.push_back(temp);
}

Eigen::MatrixXd VMC_C20::jacobian(const VectorXd &axialCoor) const
{
    double X = axialCoor(0);
    double Y = axialCoor(1);
    double Z = axialCoor(2);
    double A = DEG2RAD(axialCoor(3) );
    double C = DEG2RAD(axialCoor(4) );
    double Sa = sin(A);
    double Ca = cos(A);
    double Sc = sin(C);
    double Cc = cos(C);
    MatrixXd mat(6, 5); // \frac{\partial x}{\partial X}, etc.
    Vector3d v;
    v << X-parameters(0), Y+parameters(1), Z-parameters(2);
    mat << Cc, Ca*Sc, Sa*Sc, Ca*Sc*v(2)-Sa*Sc*v(1), Cc*Sa*v(2)-Sc*v(0)+Ca*Cc*v(1),
            Sc, -Ca*Cc, -Cc*Sa, Cc*Sa*v(1)-Ca*Cc*v(2), Cc*v(0)+Sa*Sc*v(2)+Ca*Sc*v(1),
            0, -Sa, Ca, -Ca*v(1)-Sa*v(2), 0,
            0, 0, 0, Ca*Sc, Cc*Sa,
            0, 0, 0, -Ca*Cc, Sa*Sc,
            0, 0, 0, -Sa, 0;
    return mat;
}

Eigen::MatrixXd VMC_C20::inverseJacobian(const VectorXd &axialCoor) const
{
    MatrixXd mat = jacobian(axialCoor);
    /// The Moore-Penrose inverse is derieved by the SVD decomposition.
    /// M = UDV^T
    /// M^+ = VDU^T.
    JacobiSVD<MatrixXd> svd(mat, ComputeThinU | ComputeThinV);
    uint8_t rank = svd.singularValues().size();
    MatrixXd Dp = MatrixXd::Zero(5, 5);
    for (uint8_t i=0; i<rank; i++)
        Dp(i, i) = 1.0 / svd.singularValues()(i);
    MatrixXd mi(5, 6);
    mi = svd.matrixV() * Dp * svd.matrixU().transpose();
    return mi;
}

void VMC_C20::kinematicConstraintMapping(const VectorXd &axialCoor,
                                         const VectorXd &axialCoorNext,
                                         VectorXd &posConstraint,
                                         VectorXd &oriConstraint) const
{
    /// Axial movement.
    /// Unit of axial movement is not relevant.
    VectorXd delta = axialCoorNext - axialCoor;
    /// Eq. (14)
    double posTemp = delta.head(3).norm();
    double A = DEG2RAD(axialCoor(3) );
    double to = sin(A)*delta(4);
    double oriTemp = sqrt(delta(3)*delta(3) + to*to);
    double pm;
    double om;
    double t1; // temporary variable.
    for (uint8_t i=0; i<kinematicConstraint.cols(); i++)
    {
        /// Reset pm and om to some large values.
        pm = std::numeric_limits<double>::max();
        om = std::numeric_limits<double>::max();
        for (uint8_t j=0; j<3; j++)
        {
            double da = fabs(delta(j) );
            if (da == 0)
                t1 = kinematicConstraint(j, i);
            else
                t1 = kinematicConstraint(j, i) / da;
            if (t1 < pm)
                pm = t1;
        }
        posConstraint(i) = pm * posTemp;
        for (uint8_t j=3; j<5; j++)
        {
            double da = fabs(delta(j) );
            if (da == 0)
                t1 = kinematicConstraint(j, i);
            else
                t1 = kinematicConstraint(j, i) / da;
            if (t1 < om)
                om = t1;
        }
        oriConstraint(i) = om * oriTemp;
    }
}
