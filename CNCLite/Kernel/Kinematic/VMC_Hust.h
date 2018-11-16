#ifndef VMC_HUST_H
#define VMC_HUST_H
#include "Kinematic.h"

namespace CNCLite{
class VMC_Hust : public MachineTool
{
public:
    VMC_Hust();
    VMC_Hust(const Eigen::VectorXd& param, const Eigen::VectorXd& lb,
             const Eigen::VectorXd& ub);
    VMC_Hust(const Eigen::VectorXd &param, const Eigen::VectorXd &lb,
             const Eigen::VectorXd &ub, const Eigen::MatrixXd& kineConstraint);
    ~VMC_Hust();

public:
    Eigen::Vector3d translationalCoor(const Eigen::Vector3d &pos,
                                      const Eigen::Vector2d &rot) const;
    Eigen::Vector3d orientationCoor(const Eigen::Vector2d &rot) const;

    Eigen::Vector3d positionCoor(const Eigen::Vector3d &trans,
                                 const Eigen::Vector2d &rot) const;
    Eigen::MatrixXd jacobian(const Eigen::VectorXd &axialCoor) const;
    Eigen::MatrixXd inverseJacobian(const Eigen::VectorXd &axialCoor) const;

    void kinematicConstraintMapping(const Eigen::VectorXd &axialCoor,
                                    const Eigen::VectorXd &axialCoorNext,
                                    Eigen::VectorXd &posConstraint,
                                    Eigen::VectorXd &oriConstraint) const;
protected:
    Eigen::Vector2d oneRot(const Eigen::Vector3d &ori) const;
    void allRot(const Eigen::Vector3d &ori, std::vector<Eigen::Vector2d> &rotSeries) const;
};

} // End of namespace CNCLite
#endif // VMC_HUST_H
