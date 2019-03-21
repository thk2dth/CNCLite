#include "Planner.h"

using namespace Eigen;
using namespace CNCLite;

Planner::Planner(unsigned int size, double t, unsigned int ord) :
    bufferSize(size), ts(t), order(ord)
{
    isPoolFull = false;
    isWriteCompleted = false;
    noUpdateIndex = 0;
    isNeedUpdate = false;
    isUpdateCompleted = false;
    Trajectory* traj;
    currentOutputIndex = 0;
    /// Initially, each pool consists of bufferSize-1 trajectories.
    for (unsigned int i=0; i<bufferSize-1; i++)
    {
        traj = new Trajectory();
        trajectoryPool.push_back(traj);
    }
}

bool Planner::read(Trajectory *traj)
{
    if (isUpdateCompleted)
    {
        *traj = *trajectoryPool.front();
        trajectoryPool.pop_front();
        isPoolFull = false;
        currentOutputIndex++;
        return true;
    }
    else
    {
        std::cerr << "Buffer updating is not completed.\n";
        return false;
    }
}

bool Planner::write(Trajectory* traj)
{
    /// Write the trajectory data to the head.
    trajectoryPool.push_back(traj);
    isUpdateCompleted = false;
    if(trajectoryPool.size() == bufferSize)
        isPoolFull = true;
    /// The buffer needs update in the following two conditions.
    ///     (1) Buffer is full and write has not completed.
    ///     (2) Write is completed. When write is completed, function write() will not be called.
    ///         Thus, the following situation will happen only once.
    if( (isPoolFull && !isWriteCompleted) || isWriteCompleted )
        isNeedUpdate = true;
    return true;
}

bool Planner::updateBuffer()
{
    /// If the boundPool is full, or no new trajectory data is available, then update the feedPool.
    /// The feedPool is not a circular buffer, and it always start from index 0.
    if (isNeedUpdate)
    {
        backwardScanning();
        forwardScanning();
        isUpdateCompleted = true;
    }
    return true;
}

bool Planner::backwardScanning()
{
    for (int i = bufferSize; i>=1; i--)
    {
        /// Backward scanning if vs > ve
        if (trajectoryPool.at(i)->getVs() > trajectoryPool.at(i)->getVe() )
        {
            if ( trajectoryPool[i]->plan() )
            {
                trajectoryPool[i-1]->setVe(trajectoryPool.at(i)->getVs() );
            }
        }
    }
    return true;
}

bool Planner::forwardScanning()
{
    for (int i=0; i<bufferSize; i++)
    {
        /// Forward scanning if ve > vs
        if (trajectoryPool.at(i)->getVe() > trajectoryPool.at(i)->getVs() )
        {
            if (trajectoryPool[i]->plan() )
            {
                trajectoryPool[i+1]->setVs(trajectoryPool.at(i)->getVe() );
            }
        }
    }
    return true;
}

bool Planner::SLATP(std::vector<std::deque<Trajectory *> > trajList)
{
    /// Number of paths to be synchronized.
    unsigned int np = trajList.size();
    /// Number of trajectories each path consisting of.
    unsigned int nt = trajList.front().size();
    int i = 0;
    int j = 0;
    double tmax = 0.0;
    /// Backward scanning with time synchronization.
    for (i=nt-1; i>=0; i--)
    {
        tmax = 0.0;
        for (j=0; j<np; j++)
        {
            if (trajList[j][i]->plan() )
            {
                /// The end velocity of the previous trajectory should be decreased too.
                if (i > 0 &&
                        trajList[j][i]->getVs() < trajList[j][i-1]->getVe() )
                    trajList[j][i-1]->setVe(trajList[j][i]->getVs() );
                /// The start velocity of the next trajectory should be decreased too.
                if (i < nt -1  &&
                        trajList[j][i]->getVe() < trajList[j][i+1]->getVs() )
                    trajList[j][i+1]->setVs(trajList[j][i]->getVe() );

            }
            if (trajList[j][i]->getDuration() > tmax)
                tmax = trajList[j][i]->getDuration();
        }
        /// Time synchronization.
        for (j=0; j<np; j++)
        {
            /// If vs>ve, synchronize the trajectory during backward scanning.
            /// Otherwise, synchronize the trajectory during forward scanning.
            if (tmax > trajList[j][i]->getDuration() &&
                    trajList[j][i]->getVs() > trajList[j][i]->getVe() )
            {
                if (trajList[j][i]->synchronize(tmax) && i > 0)
                {
                    trajList[j][i-1]->setVe(trajList[j][i]->getVs() );
                }
            }
        }
    }

    /// Forward scanning with time synchronization.
    for (i=0; i<nt; i++)
    {
        tmax = 0.0;
        for (j=0; j<np; j++)
        {
            if (trajList[j][i]->plan() )
            {
                /// The start velocity of the next trajectory should be decreased too.
                if (i < nt -1 &&
                        trajList[j][i]->getVe() < trajList[j][i+1]->getVs() )
                    trajList[j][i+1]->setVs(trajList[j][i]->getVe() );
                /// The end velocity of the previous trajectory should be decreased too.
                if (i > 0 &&
                        trajList[j][i]->getVs() < trajList[j][i-1]->getVe() )
                {
                    trajList[j][i-1]->setVe(trajList[j][i]->getVs() );
                    std::cerr << "During forward scanning, the end velocity of the previous "
                                 "trajectory is decreased. However, this is not allowed." << std::endl;
                }
            }
            if (trajList[j][i]->getDuration() > tmax)
                tmax = trajList[j][i]->getDuration();
        }
        /// Time synchronization.
        for (j=0; j<np; j++)
        {
            /// During forward scanning, time synchronization is always applied,
            /// regardless of the relationship between the boundary velocities.
            if (tmax > trajList[j][i]->getDuration() )
            {
                if (trajList[j][i]->synchronize(tmax) )
                {
                    if (i < nt -1 &&
                            trajList[j][i]->getVe() < trajList[j][i+1]->getVs() )
                        trajList[j][i+1]->setVs(trajList[j][i]->getVe() );
                    /// The end velocity of the previous trajectory should be decreased too.
                    if (i > 0 &&
                            trajList[j][i]->getVs() < trajList[j][i-1]->getVe() )
                    {
                        trajList[j][i-1]->setVe(trajList[j][i]->getVs() );
                        std::cerr << "During forward scanning, the end velocity of the previous "
                                     "trajectory is decreased. However, this is not allowed." << std::endl;
                    }
                }
            }
        }
    }
    return true;
}

bool Planner::SLATPEx(std::vector<std::deque<Trajectory *> > trajList)
{
    /// Number of paths to be synchronized.
    uint32_t np = trajList.size();
    /// Number of trajectories each path consisting of.
    uint8_t nt = trajList.front().size();
    uint32_t i = 0;
    uint8_t j = 0;
    double tmax = 0.0;
    /// Plan the trajectories of the remaining lines.
    for (i = 0; i < nt; i++)
    {
        j = 0; // reset j to 0.
        if (trajList[j][i]->getPath()->getCurveType() != Path::BSPLINELINEARBLEND
                && trajList[j][i]->getPath()->getCurveType() != Path::BSPLINEANGULARBLEND)
        {
            tmax = 0.0;
            for (j = 0; j < np; j++)
            {
                if (trajList[j][i]->plan() )
                {
                    /// The start boundary velocity is decreased.
                    if (i > 0 &&
                            trajList[j][i-1]->getVe() > trajList[j][i]->getVs() )
                        trajList[j][i-1]->setVe(trajList[j][i]->getVs() );
                    if (i < nt - 1 &&
                            trajList[j][i+1]->getVs() > trajList[j][i]->getVe() )
                        trajList[j][i+1]->setVs(trajList[j][i]->getVe() );
                }
                if (trajList[j][i]->getDuration() > tmax)
                    tmax = trajList[j][i]->getDuration();
            }
            /// Time Synchronization of the trajectories.
            for (j = 0; j < np; j++)
            {
                if (trajList[j][i]->getDuration() < tmax)
                {
                    if (trajList[j][i]->synchronize(tmax) )
                    {
                        /// The start boundary velocity is decreased.
                        if (i > 0 &&
                                trajList[j][i-1]->getVe() > trajList[j][i]->getVs() )
                            trajList[j][i-1]->setVe(trajList[j][i]->getVs() );
                        if (i < nt - 1 &&
                                trajList[j][i+1]->getVs() > trajList[j][i]->getVe() )
                            trajList[j][i+1]->setVs(trajList[j][i]->getVe() );
                    }
                }
            }
        }
    }
    /// Plan the trajectories of the inserted curves.
    for (i = 0; i < nt; i++)
    {
        j = 0; // reset j to 0.
        if (trajList[j][i]->getPath()->getCurveType() == Path::BSPLINELINEARBLEND
                || trajList[j][i]->getPath()->getCurveType() == Path::BSPLINEANGULARBLEND)
        {
            tmax = 0.0;
            for (j = 0; j < np; j++)
            {
                if (trajList[j][i]->plan() )
                    std::cout << "The boundary velocity of the inserted curve is decreased "
                                 "during planning, which is not allowed.\n";
                if (trajList[j][i]->getDuration() > tmax)
                    tmax = trajList[j][i]->getDuration();
            }
            /// Time Synchronization of the trajectories.
            for (j = 0; j < np; j++)
            {
                if (trajList[j][i]->getDuration() < tmax)
                {
                    if (trajList[j][i]->synchronize(tmax) )
                        std::cout << "The boundary velocity of the inserted curve is decreased "
                                     "during synchronization, which is not allowed.\n";
                }
            }
        }
    }
    return true;
}


bool Planner::BDS(std::deque<Trajectory *> trajPool)
{
    /// Number of trajectories.
    int nt = trajPool.size();
    int i = 0;
    double vs = trajPool.front()->getVs();
    /// Backward scanning.
    for (i = nt-1; i >= 0; i--)
    {
        /// If vs > ve
        if (trajPool.at(i)->getVs() > trajPool.at(i)->getVe())
            if ( trajPool[i]->plan() && i > 0)
            {
                trajPool[i-1]->setVe(trajPool.at(i)->getVs() );
                /*if (trajPool[i-1]->getKinematicConstraint()(1) == 0.0)
                {
                    /// Constant velocity
                    trajPool[i-1]->setVs(trajPool.at(i)->getVs() );
                    if (i >= 2)
                        trajPool[i-2]->setVe(trajPool.at(i)->getVs() );
                }*/
            }
    }
    /// Forward scanning
    for (i=0; i <= nt-1; i++)
    {
        /// In all situations, plan the trajectory.
        if ( trajPool[i]->plan() && i < nt-1)
        {
            trajPool[i+1]->setVs( trajPool.at(i)->getVe() );
            /*if (trajPool[i+1]->getKinematicConstraint()(1) == 0.0)
            {
                /// Constant velocity
                trajPool[i+1]->setVe(trajPool.at(i)->getVe() );
                if (i <= nt-3)
                    trajPool[i+2]->setVs(trajPool.at(i)->getVe() );
            }*/
        }
    }

    return trajPool.front()->getVs() == vs;
}
