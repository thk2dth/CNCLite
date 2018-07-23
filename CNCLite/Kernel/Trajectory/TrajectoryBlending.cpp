#include "TrajectoryBlending.h"

using namespace CNCLite;
using namespace Eigen;

TrajectoryBlending::TrajectoryBlending()
{

}

TrajectoryBlending::TrajectoryBlending(double vs, double ve,
                                       const VectorXd &kineConst, Path *path,
                                       TimeLaw::LawType type,
                                       double us, double ue)
    : Trajectory(vs, ve, kineConst, path, type, us, ue)
{

}

TrajectoryBlending::~TrajectoryBlending()
{
    for(unsigned int i=0; i<maxNumLaw; i++)
    {
        if (law.at(i) != nullptr)
        {
            delete law[i];
            law[i] = nullptr;
        }
    }
}

bool TrajectoryBlending::plan()
{
    unsigned int order = kinematicConstraint.size(); // order of the time law.
    for (int i=0; i<maxNumLaw; i++)
    {
        switch (order)
        {
        case 2:
            law[i] = new AccelerationBounded();
            break;
        case 3:
            law[i] = new JerkBounded();
            break;
        case 4:
            law[i] = new SnapBounded();
            break;
        default:
            law[i] = new SnapBounded();
        }
    }
    law[0]->initialize(kinematicConstraint, true);
    law[1]->initialize(kinematicConstraint, true); /// Cruise phase.
    law[2]->initialize(kinematicConstraint, false);
    if (vs == ve)
    {
        law[1]->lawCruise(vs, length);
        ta = 0;
        td = 0;
        tc = length / vs;
        sa = 0;
        sd = 0;
        sc = length;
        vf = vs;
        duration = ta + tc + td;
        trajType = CRU;
        isPlanned = true;
    }
    else if (vs < ve)
    {
        sa = law[0]->lawVV(vs, ve);
        ta = law[0]->getDuration();
        sc = length - sa;
        law[1]->lawCruise(ve, sc);
        tc = law[1]->getDuration();
        vf = ve;
        sd = 0.0;
        td = 0.0;
        duration = ta + tc + td;
        trajType = ACC_CRU;
        isPlanned = true;
    }
    else
    {
        sd = law[2]->lawVV(vs, ve);
        td = law[2]->getDuration();
        sc = length - sd;
        law[1]->lawCruise(vs, sc);
        tc = law[1]->getDuration();
        vf = vs;
        sa = 0.0;
        sd = 0.0;
        duration = ta + tc + td;
        trajType = CRU_DEC;
        isPlanned = true;
    }
    return false;
}

bool TrajectoryBlending::synchronize(double tg)
{
    if (!isPlanned)
        plan();
    if (tg <= duration)
    {
        std::cout << "Cannot synchronize the traversing time to a shorter duration.\n";
        return true;
    }
    /*double T0 = 2 * length / (vs + ve); // Duration threshold.
    if (Tg < T0)
    {
        double Tl = 0.0;
        if (ve > )
    }*/
    /*******************************************************************
     * Binary search is not robust because the totol duration is not
     * monotonically increasing with respect to vf.*/
    /// Use binary search to determine vf. vf belongs to [0, high].
    /// Before searching, the existence of the cruise phase should be guaranteed.
    double high = fmax(vs, ve);
    double low = 0.0;
    do
    {
        vf = 0.5 * (high + low); // vf can be decreased to 0.
        /// Note law[0] is not necessary an acceleration phase.
        sa = law[0]->lawVV(vs, vf);
        sd = law[2]->lawVV(vf, ve);
        ta = law[0]->getDuration();
        td = law[2]->getDuration();
        tc = (length - sa - sd) / vf;
        duration = ta + tc + td;
        if (duration > tg && tc >= 0)
            /// Total duration is too long. vf should be increased.
            low = vf;
        else
        {
            /*if (ta + td > tg)
                low = vf; // Increase vf to shorten ta and td.
            else
                high = vf;*/
            high = vf;
        }
        // vf = 0.5 * (low + high);
    }while (fabs(duration - tg) >= EPS_NUM && vf >= EPS_NUM);
    sa = law[0]->lawVV(vs, vf);
    ta = law[0]->getDuration();
    sd = law[2]->lawVV(vf, ve);
    td = law[2]->getDuration();
    sc = length - sa - sd;
    if (vf >= EPS_NUM && fabs(sc) >= EPS_CAD)
    {
        law[1]->lawCruise(vf, sc);
        tc = law[1]->getDuration();
    }
    else
    {
        /// Cruise phase with zero velocity.
        tc = tg - ta -td;
        law[1]->lawCruiseDuration(vf, tc);
    }
    if (vs < vf)
        trajType = DEC_CRU_ACC;
    else
        trajType= ACC_CRU_ACC;
    duration = ta + tc + td;
    /* Binary search is not robust because the totol duration is not
     * monotonically increasing with respect to vf.
    ***************************************************************************/
    return false;
}
