#ifndef KINEMATIC_H
#define KINEMATIC_H
#include "../../../ThirdParty/Eigen/Dense"
#include "../../../Common/Common.hpp"
#include <cmath>
#include <cassert>
#include <vector>

namespace CNCLite {
/**********************************************
 * Abstract class Kinematic.
 * ********************************************/
class Kinematic
{
public:
    Kinematic();
    /// \param param parameters of the kinematic model;
    /// \param lb lower bound of the axial strides;
    /// \param up upper bound of the axial strides.
    /// \param kineConstraint axial kinematic constraints, each column is
    ///         the same order constraints for all axes.
    Kinematic(const Eigen::VectorXd& param, const Eigen::VectorXd& lb,
              const Eigen::VectorXd& ub);
    Kinematic(const Eigen::VectorXd &param, const Eigen::VectorXd &lb,
              const Eigen::VectorXd &ub, const Eigen::MatrixXd& kineConstraint);
    virtual ~Kinematic();

public:
    /// Forward kinematic transformation: from joint space to task space.
    virtual Eigen::VectorXd forwardKinematic(const Eigen::VectorXd& axialCoor,
                                             const Eigen::VectorXd& axialCoorPrevious) const = 0;
    virtual Eigen::VectorXd inverseKinematic(const Eigen::VectorXd& cartesianCoor,
                                             const Eigen::VectorXd& axialCoorPrevious) const = 0;
    virtual Eigen::MatrixXd jacobian(const Eigen::VectorXd& axialCoor) const = 0;
    virtual Eigen::MatrixXd inverseJacobian(const Eigen::VectorXd& axialCoor) const = 0;
    virtual void kinematicConstraintMapping(const Eigen::VectorXd& axialCoor,
                                            const Eigen::VectorXd& axialCoorNext,
                                            Eigen::VectorXd& posConstraint,
                                            Eigen::VectorXd& oriConstraint) const = 0;

protected:
    Eigen::VectorXd parameters;
    Eigen::VectorXd lowerBound, upperBound;
    Eigen::MatrixXd kinematicConstraint;
    unsigned int dof;
};

/// Kinematic transformations of machine tools.
/// Machine tools usually consist of three translational axes and two rotary axes.
/// The task space is represented by position and orientation.
class MachineTool : public Kinematic
{
public:
    MachineTool();
    MachineTool(const Eigen::VectorXd& param, const Eigen::VectorXd& lb,
                const Eigen::VectorXd& ub);
    MachineTool(const Eigen::VectorXd &param, const Eigen::VectorXd &lb,
                const Eigen::VectorXd &ub, const Eigen::MatrixXd& kineConstraint);
    virtual ~MachineTool();
public:
    virtual Eigen::VectorXd forwardKinematic(const Eigen::VectorXd &axialCoor,
                                             const Eigen::VectorXd &axialCoorPrevious) const;
    virtual Eigen::VectorXd inverseKinematic(const Eigen::VectorXd &cartesianCoor,
                                             const Eigen::VectorXd &axialCoorPrevious) const;
    virtual Eigen::MatrixXd jacobian(const Eigen::VectorXd &axialCoor) const = 0;
    virtual Eigen::MatrixXd inverseJacobian(const Eigen::VectorXd &axialCoor) const = 0;
    virtual void kinematicConstraintMapping(const Eigen::VectorXd &axialCoor,
                                            const Eigen::VectorXd &axialCoorNext,
                                            Eigen::VectorXd &posConstraint,
                                            Eigen::VectorXd &oriConstraint) const = 0;

public:
    /// Calculate the translational axes by the determined rotaty axes.
    /// This function acts like the RTCP function.
    virtual Eigen::Vector3d translationalCoor(const Eigen::Vector3d& pos,
                                              const Eigen::Vector2d& rot) const = 0;
    /// Map the coordinates of the rotary axes to an orientation in task space.
    virtual Eigen::Vector3d orientationCoor(const Eigen::Vector2d& rot) const = 0;

    /// Calculate the tool position by the machine coordinates.
    /// This function is similar to function forwardKinematic.
    virtual Eigen::Vector3d positionCoor(const Eigen::Vector3d& trans,
                                         const Eigen::Vector2d& rot) const = 0;

    /// Calculate the coordinates of rotary axes.
    /// Return values in degree.
    virtual Eigen::Vector2d rotaryCoor(const Eigen::Vector3d& ori,
                                       const Eigen::Vector2d& rotPre) const;

protected:
    virtual Eigen::Vector2d oneRot(const Eigen::Vector3d& ori) const = 0;
    virtual void allRot(const Eigen::Vector3d& ori,
                        std::vector<Eigen::Vector2d>& rotSeries ) const = 0;
    /// Remove angles which exceeds the bounds of the rotational axes.
    void removeRotByBound(const std::vector<Eigen::Vector2d> &rotSeries,
                          std::vector<Eigen::Vector2d>& rotSeriesPruned) const;
    /// Select an rotational coordinate which is nearest to the rotPre coordinate.
    /// The distance is defined by the Euclidean distance.
    /// This function is neccessary if other selection method is required.
    /// This function is usually called after calling the the removeRotByBound() function.
    /// If there is no limits on the axial strides, it should be called directly.
    Eigen::Vector2d selectRot(const std::vector<Eigen::Vector2d> &rotSeries,
                              const Eigen::Vector2d &rotPre) const;
};


/// Specific machine tool.
class VMC_C20 : public MachineTool
{
public:
    VMC_C20();
    VMC_C20(const Eigen::VectorXd& param, const Eigen::VectorXd& lb,
            const Eigen::VectorXd& ub);
    VMC_C20(const Eigen::VectorXd &param, const Eigen::VectorXd &lb,
            const Eigen::VectorXd &ub, const Eigen::MatrixXd& kineConstraint);
    ~VMC_C20();
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
    Eigen::Vector2d oneRot(const Eigen::Vector3d& ori) const;
    void allRot(const Eigen::Vector3d& ori,
                std::vector<Eigen::Vector2d>& rotSeries ) const;
};

} // End of namespace CNCLite.

#endif // KINEMATIC_H
