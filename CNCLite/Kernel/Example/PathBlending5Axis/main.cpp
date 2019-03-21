#include <QCoreApplication>
#include "../../../Kinematic/Kinematic.h"
#include "../../../Kinematic/AptParser.h"
#include "../../../Planner/Planner.h"
#include "../../Path/Line.h"
#include "../../Path/Arc.h"
#include "../../Path/BsplineLinearBlend.h"
#include "../../Path/BsplineAngularBlend.h"
#include "../../Trajectory/Trajectory.h"
#include "../../Trajectory/TrajectoryBlending.h"
#include <fstream>
#include <QTime>

using namespace CNCLite;
using namespace Eigen;

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);
    QTime timer;
    /// Kinematic model.
    uint8_t dim = 5;
    VectorXd paraMachine(3);
    VectorXd lb(dim), ub(dim);
    paraMachine << 0.0, 0.0, 172.0;
    /// No limit on C axis.
    lb << 0, 0, 0, -30, std::numeric_limits<double>::lowest();
    ub << 600, 500, 250, 100, std::numeric_limits<double>::max();
    /// Kinematic constraints.
    uint8_t order = 3;
    VectorXd constLin(order), constOri(order);
    constLin << 50, 200, 2000;
    constOri << DEG2RAD(5), DEG2RAD(50), DEG2RAD(500);
    MatrixXd kineConstraint(dim, order);
    kineConstraint << 50, 200, 2000,
            50, 200, 2000,
            30, 150, 1500,
            DEG2RAD(5), DEG2RAD(20), DEG2RAD(200),
            DEG2RAD(10), DEG2RAD(50), DEG2RAD(500);
    /// Samping period.
    const double ts = 0.0044; /// sampling period, in s.
    double tolLin = 0.04; // Blending error for linear path, in mm.
    double tolOri = DEG2RAD(0.04); // Blending error for angular path, in radius.
    double cRatio = 0.25;
    double mRatio = 1 / 3.0;
    double ceLin = 0.01; // Chord error when traversing parametric curves.
    double ceOri = DEG2RAD(0.01); // Chord error for tool orientation, in radius.
    //    constLin << 50, 200, 2000;
    //    constOri << DEG2RAD(5), DEG2RAD(20), DEG2RAD(200); // in degree.
    VMC_C20 vmc(paraMachine, lb, ub, kineConstraint);

    /// Cutter location data.
    QString fn("../../../../NC Program/APT/Square_apt.apt");
    /// QString fn("../../../../NC Program/APT/Tulsyan_Altintas_IJMTM_2015.apt");
    /// QString fn("../../../../NC Program/APT/Iso_apt_S.apt");
    AptParser apt(fn);
    MatrixXd clWCS;
    uint8_t interval = 1;
    apt.readAll(clWCS, interval);
    MatrixXd clMCS(clWCS.rows(), dim);
    clMCS.row(0) = vmc.inverseKinematic(clWCS.row(0), VectorXd::Zero(5) );
    for (uint i=1; i<clWCS.rows(); i++)
    {
        clMCS.row(i) = vmc.inverseKinematic(clWCS.row(i), clMCS.row(i-1) );
    }

    std::cout << "Storing the WCS and MCS data...\n";
    std::ofstream outWCS("../../PathBlending5Axis/wcsData.txt");
    std::ofstream outMCS("../../PathBlending5Axis/mcsData.txt");
    if (outWCS.is_open() && outMCS.is_open() )
    {
        for (uint i=0; i<clWCS.rows(); i++)
        {
            outWCS << clWCS(i, 0) << " " << clWCS(i, 1) << " "
                   << clWCS(i, 2) << " " << clWCS(i, 3) << " "
                   << clWCS(i, 4) << " " << clWCS(i, 5) << std::endl;
            outMCS << clMCS(i, 0) << " " << clMCS(i, 1) << " "
                   << clMCS(i, 2) << " " << clMCS(i, 3) << " "
                   << clMCS(i, 4) << std::endl;
        }
    }
    outWCS.close();
    outMCS.close();

    uint32_t numberTraj = 2*clWCS.rows() - 3; /// Number of trajectories.
    std::vector<Path* > pathLin(numberTraj, nullptr);
    std::vector<Path* > pathOri(numberTraj, nullptr);
    /// Blend the corner forming by three points.
    VectorXd paramBlendLinear(3);
    VectorXd paramBlendAngular(3);
    paramBlendLinear << mRatio, cRatio, tolLin;
    paramBlendAngular << mRatio, cRatio, tolOri;
    std::vector<VectorXd> ctrlPtsLinear(4, VectorXd::Zero(dim) );
    std::vector<VectorXd> ctrlPtsAngular(4, VectorXd::Zero(dim) );
    ctrlPtsLinear[3] = paramBlendLinear;
    ctrlPtsAngular[3] = paramBlendAngular;

    std::cout << "Begin to blend the path..." << std::endl;
    timer.start();
    for (uint32_t i = 0; i < clWCS.rows()-2; i++)
    {
        /// Inserted Bspline.
        ctrlPtsLinear[0] = clWCS.row(i).head(3);
        ctrlPtsLinear[1] = clWCS.row(i+1).head(3);
        ctrlPtsLinear[2] = clWCS.row(i+2).head(3);
        ctrlPtsAngular[0] = clWCS.row(i).tail(3);
        ctrlPtsAngular[1] = clWCS.row(i+1).tail(3);
        ctrlPtsAngular[2] = clWCS.row(i+2).tail(3);
        pathLin[2*i + 1] = new BsplineLinearBlend(ctrlPtsLinear);
        pathOri[2*i + 1] = new BsplineAngularBlend(ctrlPtsAngular);
    }

    /// Then insert lines.
    /// The remaining lines and arcs have not been converted to Nurbs curves yet.
    VectorXd ps, pe;
    for (uint32_t i = 0; i <= clWCS.rows()-2; i++)
    {
        /// First line.
        if (0 == i)
        {
            ps = clWCS.row(0).head(3);
            pe = pathLin.at(1)->getBeginPoint();
            pathLin[0] = new Line(ps, pe);
            ps = clWCS.row(0).tail(3);
            pe = pathOri.at(1)->getBeginPoint();
            /// Start, end, and center points, respectively.
            pathOri[0] = new Arc(ps, pe, VectorXd::Zero(3) );
        }
        else if (i == clWCS.rows() - 2)
        {
            ps = pathLin.at(2*(clWCS.rows()-3) + 1)->getEndPoint();
            pe = clWCS.row(clWCS.rows() - 1).head(3); // Last point.
            pathLin[2*clWCS.rows() - 4] = new Line(ps, pe);
            ps = pathOri.at(2*(clWCS.rows()-3) + 1)->getEndPoint();
            pe = clWCS.row(clWCS.rows() - 1).tail(3); // Last point.
            pathOri[2*clWCS.rows() - 4] = new Arc(ps, pe, VectorXd::Zero(3) );
        }
        else
        {
            ps = pathLin.at(2*i - 1)->getEndPoint();
            pe = pathLin.at(2*i + 1)->getBeginPoint();
            pathLin[2*i] = new Line(ps, pe);
            ps = pathOri.at(2*i - 1)->getEndPoint();
            pe = pathOri.at(2*i + 1)->getBeginPoint();
            pathOri[2*i] = new Arc(ps, pe, VectorXd::Zero(3) );
        }
    }

    int tCost = timer.elapsed();
    std::cout << "Path blending finished in " << tCost << "ms. \n";

    std::deque<Trajectory* > trajLin(numberTraj, nullptr);
    std::deque<Trajectory* > trajOri(numberTraj, nullptr);
    std::vector< std::deque<Trajectory* > > trajList;

    std::cout << "Initializing the trajectories..." << std::endl;
    timer.restart();
    double vsLin, veLin, vsOri, veOri, vr;
    double maxCurv, uMid;
    Path* tp = nullptr; // temporary pointer.
    // tempLaw is used to check the feasibility of the boundary velocities.
    TimeLaw* tempLaw;
    switch (order)
    {
    case 2:
        tempLaw = new AccelerationBounded(); break;
    case 3:
        tempLaw = new JerkBounded(); break;
    case 4:
        tempLaw = new SnapBounded(); break;
    }

    double ks = 0.4;
    /// Trajectory for the inserted parametric curves.
    for (uint32_t i = 0; i < clWCS.rows()-2; i++)
    {
        /// Linear path
        tp = pathLin[2*i + 1];
        uMid = 0.5 * (tp->getBeginParameter() +
                      tp->getEndParameter() );
        /// Maximum curvature is attained at the middle of the parametric curve.
        maxCurv = tp->calculateCurvature(uMid);
        vsLin = Path::calculateTraverseVelocityR3(maxCurv, constLin(0), constLin(1),
                                                  constLin(2), ceLin, ts);
        /// Compatibility check.
        tempLaw->initialize(constLin, true);
        vr = tempLaw->lawVD(0.0, ks*tp->getLength() );
        vsLin = fmin(vsLin, vr);
        trajLin[2*i + 1] = new TrajectoryBlending(vsLin, vsLin, constLin, tp);

        /// Angular path
        tp = pathOri[2*i + 1];
        uMid = 0.5 * (tp->getBeginParameter() +
                      tp->getEndParameter() );
        /// Maximum curvature is attained at the middle of the parametric curve.
        maxCurv = tp->calculateCurvature(uMid);
        vsOri = Path::calculateTraverseVelocityS2(maxCurv, constOri(0), constOri(1),
                                                  constOri(2), ceOri, ts);
        /// Compatibility check.
        tempLaw->initialize(constOri, true);
        vr = tempLaw->lawVD(0.0, ks*tp->getLength() );
        vsOri = fmin(vsOri, vr);
        trajOri[2*i + 1] = new TrajectoryBlending(vsOri, vsOri, constOri, tp);
    }

    /// Trajectory for the remaining curves.
    for (uint32_t i = 0; i <= clWCS.rows()-2; i++)
    {
        if (0 == i)
        {
            vsLin = 0.0;
            vsOri = 0.0;
        }
        else
        {
            vsLin = trajLin.at(2*i - 1)->getVe();
            vsOri = trajOri.at(2*i - 1)->getVe();
        }
        if (i == clWCS.rows() - 2)
        {
            veLin = 0.0;
            veOri = 0.0;
        }
        else
        {
            veLin = trajLin.at(2*i + 1)->getVs();
            veOri = trajOri.at(2*i + 1)->getVs();
        }

        /// Compatibility check.
        /// Only the lower one of the two boundary velocities should be checked.
        tempLaw->initialize(constLin, true);
        vr = tempLaw->lawVD(0.0, 0.5*pathLin[2*i]->getLength() );
        if (veLin > vsLin)
        {
            if (vsLin > vr)
            {
                // Start velocity should be decreased.
                vsLin = vr;
                trajLin[2*i - 1]->setVe(vr);
            }
        }
        else
        {
            if (veLin > vr)
            {
                // End velocity should be decreased.
                veLin = vr;
                trajLin[2*i + 1]->setVs(vr);
            }
        }

        tempLaw->initialize(constOri, true);
        vr = tempLaw->lawVD(0.0, 0.5*pathOri[2*i]->getLength() );
        if (veOri > vsOri)
        {
            if (vsOri > vr)
            {
                // Start velocity should be decreased.
                vsOri = vr;
                trajOri[2*i - 1]->setVe(vr);
            }
        }
        else
        {
            if (veOri > vr)
            {
                // End velocity should be decreased.
                veOri = vr;
                trajOri[2*i + 1]->setVs(vr);
            }
        }

        trajLin[2*i] = new Trajectory(vsLin, veLin, constLin, pathLin[2*i] );
        trajOri[2*i] = new Trajectory(vsOri, veOri, constOri, pathOri[2*i] );
    }
    trajList.push_back(trajLin);
    trajList.push_back(trajOri);
    tCost = timer.elapsed();
    std::cout << "Trajectory initialization finished in " << tCost << "ms. \n";

    std::cout << "Planning the motions with SLATP....\n";
    timer.restart();
    Planner::SLATPEx(trajList);
    tCost = timer.elapsed();
    std::cout << "Planning finished in " << tCost << "ms. \n";

    std::cout << "Storing the trajectory data of SLATP...\n";
    std::ofstream out("../../PathBlending5Axis/traj5Axis.txt");
    double tAll = 0.0;
    double tLin = 0.0;
    double tOri = 0.0;
    uint idxLin = 0;
    uint idxOri = 0;
    double dLin, vLin, aLin, jLin, sLin;
    double dOri, vOri, aOri, jOri, sOri;
    double dAllLin = 0.0;
    double dAllOri = 0.0;
    VectorXd ptLin, ptOri; // current position in WCS.
    VectorXd coorWCS(6), coorMCS(5); // current coordinate in WCS and MCS.
    if (out.is_open() )
    {
        while (1)
        {
            while (idxLin < numberTraj &&
                   tLin > trajList.at(0).at(idxLin)->getDuration() )
            {
                dAllLin += trajList.at(0).at(idxLin)->getLength();
                tLin -= trajList.at(0).at(idxLin)->getDuration();
                /*ds = vLin*ts + aLin*ts*ts/2 + jLin*pow(ts, 3)/6;
                    ds = ds + dLin - trajList.at(0).at(idxLin)->getLength();
                    tLin = ds / trajList.at(0).at(idxLin)->getVe();
                    tLin = tLin > 0 ? tLin : 0.0;*/
                idxLin++;
                out << std::endl;
            }
            while (idxOri < numberTraj &&
                   tOri > trajList.at(1).at(idxOri)->getDuration() )
            {
                dAllOri += RAD2DEG(trajList.at(1).at(idxOri)->getLength() );
                tOri -= trajList.at(1).at(idxOri)->getDuration();
                /*ds = vOri*ts + aOri*ts*ts/2 + jOri*pow(ts, 3)/6;
                    ds = ds + dOri - trajList.at(1).at(idxOri)->getLength();
                    tOri = ds / trajList.at(1).at(idxOri)->getVe();
                    tOri = tOri > 0 ? tOri : 0.0;*/
                idxOri++;
                out << std::endl;
            }
            if (idxLin == numberTraj || idxOri == numberTraj)
                break;
            dLin = trajList.at(0).at(idxLin)->displacement(tLin);
            vLin = trajList.at(0).at(idxLin)->velocity(tLin);
            aLin = trajList.at(0).at(idxLin)->acceleration(tLin);
            jLin = trajList.at(0).at(idxLin)->jerk(tLin);
            sLin = trajList.at(0).at(idxLin)->snap(tLin);
            dOri = trajList.at(1).at(idxOri)->displacement(tOri);
            vOri = trajList.at(1).at(idxOri)->velocity(tOri);
            aOri = trajList.at(1).at(idxOri)->acceleration(tOri);
            jOri = trajList.at(1).at(idxOri)->jerk(tOri);
            sOri = trajList.at(1).at(idxOri)->snap(tOri);

            ptLin = trajList.at(0).at(idxLin)->axialPosition(tLin);
            ptOri = trajList.at(1).at(idxOri)->axialPosition(tOri);
            coorWCS.head(3) = ptLin;
            coorWCS.tail(3) = ptOri;
            if (0 == tLin)
                coorMCS = vmc.inverseKinematic(coorWCS,
                                               VectorXd::Zero(dim) );
            else
                coorMCS = vmc.inverseKinematic(coorWCS, coorMCS);
            out << tAll << " " << dLin+dAllLin << " " << vLin << " "
                << aLin << " " << jLin << " " << sLin << " "
                << RAD2DEG(dOri)+dAllOri << " " << RAD2DEG(vOri) << " "
                << RAD2DEG(aOri) << " " << RAD2DEG(jOri) << " "
                << RAD2DEG(sOri) << " "
                << coorWCS(0) << " " << coorWCS(1) << " " << coorWCS(2) << " "
                << coorWCS(3) << " " << coorWCS(4) << " " << coorWCS(5) << " "
                << coorMCS(0) << " " << coorMCS(1) << " " << coorMCS(2) << " "
                << coorMCS(3) << " " << coorMCS(4) << std::endl;
            tAll += ts;
            tLin += ts;
            tOri += ts;
        }
    }
    out.close();

    delete tempLaw;
    tempLaw = nullptr;
    /// Deallocate the path and trajectory memory.
    for (uint32_t i=0; i<numberTraj; i++)
    {
        if (pathLin.at(i) != nullptr)
        {
            delete pathLin[i];
            pathLin[i] = nullptr;
        }
        if (pathOri.at(i) != nullptr)
        {
            delete pathOri[i];
            pathOri[i] = nullptr;
        }
        delete trajLin[i];
        trajLin[i] = nullptr;
        delete trajOri[i];
        trajOri[i] = nullptr;
    }

    std::cout << "\nEnd of Program." << std::endl;
    return a.exec();
}
