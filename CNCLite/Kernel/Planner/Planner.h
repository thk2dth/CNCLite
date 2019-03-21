#ifndef PLANNER_H
#define PLANNER_H
#include "../../ThirdParty/Eigen/Dense"
#include "../../Trajectory/Trajectory.h"
#include "../../Path/Path.h"
#include <deque>
#include <vector>

namespace CNCLite {

class Planner
{
public:
    Planner(const Planner &) = delete;
    Planner& operator = (const Planner &) = delete;
    /// Explicitly construct the buffer with its size and the feedrate function.
    /// size, buffer size; t, sampling period; order, ord of trajectory.
    explicit Planner(unsigned int size, double t, unsigned int ord);

public:
    /// Append the new boundary feedrate values to the buffer head.
    bool write(Trajectory *traj);
    /// Take away feedrate profile at the tail of the buffer.
    /// If no new boundary feedrate value is availabel, take away all the profiles after updating.
    bool read(Trajectory*  traj);
    void setIsWriteCompleted(bool value); /// Set by the interpreter when finishing interpreting.

    /// Set by the interpreter when bufferSize bound values have been fed into the buffer.
    void setIsBoundPoolFull(bool value);
    bool updateBuffer();

private:
    /// Step 1. if vs>ve, backward scanning: set the backward and after backward feedrates in the
    /// boundFeed buffer.
    bool backwardScanning();
    /// Step 2. if ve>vs, forward scanning: set the final feedrates in the boundFeed buffer.
    bool forwardScanning();

private:
    /// Buffer size of the boundary feedrates. It should be greater than 2.
    const unsigned int bufferSize;
    const double ts; /// sampling period, in second.
    const unsigned int order;
    /// The trajectoryPool stores the profiles to be taken away. It is used in backward and forward
    /// scanning steps as well.
    /// Each path consists of bufferSize trajectories.
    std::deque< Trajectory* >  trajectoryPool;
    /// If no data to be written, the feedrate in the buffer should be all taken away.
    bool isWriteCompleted;
    bool isPoolFull;
    /// The boundary feedrates from th tail to this index need no updating.
    unsigned int noUpdateIndex;
    /// When the new bound value is fed in, the feedPool needs update.
    bool isNeedUpdate;
    /// After updating, set this value to true, and it is now ready for read.
    bool isUpdateCompleted;
    /// Index of current output profile.
    unsigned int currentOutputIndex;

public:
    /// The time synchronization method.
    /// To synchronize with any duration, the trajectories to be synchorinzed should be checked under
    /// the compatible conditions.
    /// Depend on whether blending trajectory exists, there are two synchronization methods.
    /// When there is no blending trajectory, the first one is adopted. It is the default method.
    /// When the blending trajectory exists, the synchronization method shall also limit the
    /// velocity of the blending trajectory.
    /// If there is no blending trajectory, SLATP is used.
    static bool SLATP(std::vector< std::deque<Trajectory* > > trajList);
    /// If there exists blending trajectories, SLATPEx is used.
    /// In SLATPEx, the trajectories of the remaining lines are planned and synchronized first,
    /// and then the trajectories of the inserted blending curves are planned.
    /// The boundary velocities of the trajectories of the remaining lines may be decreased.
    /// The boundary velocities of the trajectories of the inserted curves are reserved.
    static bool SLATPEx(std::vector< std::deque<Trajectory* > > trajList);

    /// The traditional bi-directional scanning method.
    /// return true if the overall start velocity of the trajectory pool is changed.
    static bool BDS(std::deque<Trajectory *> trajPool);
};

} /// namespace CNCLite

#endif // PLANNER_H
