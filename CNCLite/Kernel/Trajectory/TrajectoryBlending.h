#ifndef TRAJECTORYBLENDING_H
#define TRAJECTORYBLENDING_H
#include "Trajectory.h"

/// TrajectoryBlending is used for the inserted parametric curves during curve blending.
/// Unlike the usual one, TrajectoryBlending uses a CRU trajectory during planning, and
/// a DEC-CRU-ACC trajectory during time synchronizing.
/// The start and end velocities have been limited according to the traversing distance.
/// Therefore, the two boundary velocities are always reachable.
namespace CNCLite{
class TrajectoryBlending : public Trajectory
{
public:
    TrajectoryBlending();
    TrajectoryBlending(double vs, double ve, const Eigen::VectorXd &kineConst,
                       Path *path, TimeLaw::LawType type = TimeLaw::POLYNOMIAL,
                       double us = 0.0, double ue = 1.0);
    virtual ~TrajectoryBlending();

public:
    /// Only two virtual member functions are redefined.
    /// The acceleration and higher order kinematic quantities are not zero.
    /// However, always use a constant velocity profile to plan the trajectory.
    virtual bool plan();
    /// Extend the trajectory duration to tg with a DEC-CRU-ACC trajectory.
    virtual bool synchronize(double tg);
};

} /// Namespace CNCLite

#endif // TRAJECTORYBLENDING_H
