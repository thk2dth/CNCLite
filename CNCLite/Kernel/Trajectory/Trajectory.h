#ifndef TRAJECTORY_H
#define TRAJECTORY_H
#include <vector>
#include "../../Path/Path.h"
#include "../../TimeLaw/TimeLaw.h"
#include "../../Common/Common.hpp"
#include "../../TimeLaw/AccelerationBounded.h"
#include "../../TimeLaw/JerkBounded.h"
#include "../../TimeLaw/SnapBounded.h"
#include "../../TimeLaw/JerkSineSeries.h"
#include "../../TimeLaw/JerkSine.h"

namespace CNCLite{

const unsigned int maxNumLaw = 3; /// maximum number of laws.
/// A trajectory consists of a path and at most three successive time laws.
class Trajectory
{
public:
    Trajectory();
    /// With parameters us and ue, a portion of the path can be selected. This is useful for trajectory
    /// planning of parametric curves.
    Trajectory(double vs, double ve, const Eigen::VectorXd &kineConst, Path *path,
               TimeLaw::LawType type = TimeLaw::LawType::POLYNOMIAL,double us = 0.0, double ue = 1.0);
    virtual ~Trajectory();

public:
    /// All the possible trajectory types, including the types after synchronization.
    /// Trajectories of type DEC_ACC are only used for blending curves, including linear and angular
    /// blending curves.
    enum TrajType{UNDEFINED, ACC, DEC, ACC_DEC, ACC_CRU_DEC, CRU, ACC_CRU, CRU_DEC, DEC_CRU, CRU_ACC,
                 DEC_ACC, DEC_CRU_ACC, ACC_CRU_ACC};

public:
    void initialize(double vs, double ve, const Eigen::VectorXd &kineConst, Path *path,
                    TimeLaw::LawType type = TimeLaw::LawType::POLYNOMIAL,double us = 0.0, double ue = 1.0);
    /// Return true if the boundary velocity is decreased.
    virtual bool plan(); // Use binary search method to plan the trajetory.
    /// Synchronize the trajectory with a specifed duration.
    /// Return true if the boundary velocity is decreased during synchronization.
    virtual bool synchronize(double tg);
    double displacement(double t) const;
    double velocity(double t) const;
    double acceleration(double t) const;
    double jerk(double t) const;
    double snap(double t) const;
    /// For most controller, axial position and velocity information is sufficient.
    /// Function axialPosition() will update current parameter uc.
    /// Therefore, result of axialVelocity() will be influenced by axialPosition().
    Eigen::VectorXd axialPosition(double t);
    Eigen::VectorXd axialVelocity(double t);
    /// When there is no explicit displacement function, the arc length increase is
    /// necessary for determination of the axial position and velocity.
    /// Note that unlike the above two functions, the following two functions do not
    /// update the parameters uc and dc.
    /// ds is the arc length increase.
    /// u is the curve parameter in the previous sampling period.
    /// f is the magnitude of the current tangential velocity.
    Eigen::VectorXd axialPosition(double ds, double u);
    Eigen::VectorXd axialVelocity(double ds, double u, double f);

    Eigen::VectorXd getKinematicConstraint() const;
    double getDuration() const;
    Path* getPath() const;
    double getVs() const;
    double getVe() const;
    double getVf() const;
    double getLength() const;
    void setVs(double v);
    void setVe(double v);
    /// Current parameter.
    double getUc() const;

    TrajType getTrajType() const;

    double getUs() const;

    double getUe() const;

protected:
    /// Path.
    Path *path;
    double us; // begin parameter of the path.
    double ue; // end parameter of the path.
    double length; // length of the path between us and ue.
    /// uc and dc are used for interpolation.
    double uc; // current parameter.
    double dc; // current displacement.
    double duration; // total duration of the time law.
    double vs; // begin velocity.
    double ve; // end velocity.
    double vf; // reachable velocity.
    double sa, sc, sd; // displacements during ACC, CRU and DEC phases.
    double ta, tc, td; // durations during ACC, CRU and DEC phases.
    Eigen::VectorXd kinematicConstraint;
    TrajType trajType; // Trajectory type.
    /// Time law. To facilitate the sychronization process, the trajectory always has three time laws.
    std::vector<TimeLaw* > law;
    bool isPlanned;
    TimeLaw::LawType lawType;
};

} /// namespace CNCLite


#endif // TRAJECTORY_H
