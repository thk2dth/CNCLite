#include <QCoreApplication>
#include "../../Kinematic/Kinematic.h"
#include "../../Kinematic/VMC_Hust.h"
#include "../../Kinematic/AptParser.h"
#include "../../Planner/Planner.h"
#include "../../Path/Line.h"
#include "../../Path/Arc.h"
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
    VectorXd param(3);
    VectorXd lb(dim), ub(dim);
    // param << 0.0, 0.0, 172.0; // VMC C20
    param << 0.0, 0.0, 41.0; // VMC HUST
    /// No limit on C axis.
    lb << 0, 0, 0, -30, std::numeric_limits<double>::lowest();
    ub << 600, 500, 250, 100, std::numeric_limits<double>::max();
    /// Kinematic constraints.
    uint8_t order = 3;
    VectorXd constLin(order), constOri(order);
    MatrixXd kineConstraint(dim, order);
    kineConstraint << 50, 500, 5000,
            50, 500, 5000,
            30, 300, 3000,
            DEG2RAD(5), DEG2RAD(50), DEG2RAD(500),
            DEG2RAD(20), DEG2RAD(200), DEG2RAD(2000);
    /// Samping period.
    const double ts = 0.004; /// sampling period, in s.
    double perr = 0.05; // mm
    double oerr = DEG2RAD(0.05); // radius
    //    constLin << 50, 200, 2000;
    //    constOri << DEG2RAD(5), DEG2RAD(20), DEG2RAD(200); // in degree.
    // VMC_C20 vmc(param, lb, ub, kineConstraint); // VMC C20
    VMC_Hust vmc(param, lb, ub, kineConstraint);

    /// Cutter location data.
    QString fn("../../../../NC Program/APT/Square_apt.apt");
    /// QString fn("../../../../NC Program/APT/Tulsyan_Altintas_IJMTM_2015.apt");
    /// QString fn("../../../../NC Program/APT/lsopt_apt_S_Part_Test.apt");
    AptParser apt(fn);
    MatrixXd clWCS;
    uint8_t interval = 1;
    /// uint8_t interval = 20;
    apt.readAll(clWCS, interval);
    MatrixXd clMCS(clWCS.rows(), dim);
    clMCS.row(0) = vmc.inverseKinematic(clWCS.row(0), VectorXd::Zero(5) );
    for (uint i=1; i<clWCS.rows(); i++)
    {
        clMCS.row(i) = vmc.inverseKinematic(clWCS.row(i), clMCS.row(i-1) );
    }


    std::cout << "Storing the WCS and MCS data...\n";
    std::ofstream outWCS("../../Planner/wcsData.txt");
    std::ofstream outMCS("../../Planner/mcsData.txt");
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

    std::cout << "Initializing the trajectories...\n";
    timer.start();
    /// Initialize the trajectories.
    Path* pathLin;
    Path* pathOri;
    uint numberTraj = clWCS.rows()-1;
    std::deque<Trajectory* > trajLin(numberTraj);
    std::deque<Trajectory* > trajOri(numberTraj);
    std::vector< std::deque<Trajectory* > > trajList;
    std::deque<Trajectory* > trajLinOrg(numberTraj);
    std::deque<Trajectory* > trajOriOrg(numberTraj);
    double vsLin, veLin, vsOri, veOri;
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
    double vr;

    for (uint i=0; i<numberTraj; i++)
    {
        pathLin = new Line(clWCS.row(i).head(3), clWCS.row(i+1).head(3) );
        pathOri = new Arc(clWCS.row(i).tail(3), clWCS.row(i+1).tail(3), Vector3d::Zero() );
        /// Case 1: Determine the linear and angular kinematic constraints by axial constraints.
        /// vmc.kinematicConstraintMapping(clMCS.row(i), clMCS.row(i+1), constLin, constOri);
        /// Case 2: Use constant linear and angular kinematic constraints.
        constLin << 50, 200, 2000;
        constOri << DEG2RAD(5), DEG2RAD(50), DEG2RAD(500);
        if (0 == i)
        {
            vsLin = 0.0;
            vsOri = 0.0;
        }
        else
        {
            /// Note that the previous end velocity can exceed the maximum
            /// velocity constraint of the current block. Therefore, the
            /// end velocity of the previous block and the start velocity of
            /// the current block are restricted to the maximum velocity of
            /// the current block.
            if (trajLin.at(i-1)->getVe() > constLin(0) )
                vsLin = constLin(0);
            else
                vsLin = trajLin.at(i-1)->getVe();
            if (trajOri.at(i-1)->getVe() > constOri(0) )
                vsOri = constOri(0);
            else
                vsOri = trajOri.at(i-1)->getVe();
        }
        if (i == numberTraj-1 )
        {
            veLin = 0.0;
            veOri = 0.0;
        }
        else
        {

            veLin = Line::calculateCornerVelocity(clWCS.row(i).head(3),
                                                  clWCS.row(i+1).head(3),
                                                  clWCS.row(i+2).head(3),
                                                  constLin(0),
                                                  constLin(1), ts, perr );
            veOri = Arc::calculateCornerVelocityS2(clWCS.row(i).tail(3),
                                                   clWCS.row(i+1).tail(3),
                                                   clWCS.row(i+2).tail(3),
                                                   constOri(0),
                                                   constOri(1), ts, oerr );
        }
        /// In traditional scanning, the boundary velocities are not limited by vr.
        trajLinOrg[i] = new Trajectory(vsLin, veLin, constLin, pathLin);
        if (i > 0 && vsLin < trajLinOrg.at(i-1)->getVe() )
            trajLinOrg[i-1]->setVe(vsLin); // Update previous profile.
        /// Feasibility test.
        tempLaw->initialize(constLin, true);
        vr = tempLaw->lawVD(0.0, pathLin->getLength()*0.5);
        vsLin = vsLin < vr ? vsLin : vr;
        if (i > 0 && vsLin < trajLin.at(i-1)->getVe() )
            trajLin[i-1]->setVe(vsLin); // Update previous profile.
        trajLin[i] = new Trajectory(vsLin, fmin(veLin, vr),
                                    constLin, pathLin);
        /// In traditional scanning, the boundary velocities are not limited by vr.
        trajOriOrg[i] = new Trajectory(vsOri, veOri, constOri, pathOri);
        if (i > 0 && vsOri < trajOriOrg.at(i-1)->getVe() )
            trajOriOrg[i-1]->setVe(vsOri); // Update previous profile.
        tempLaw->initialize(constOri, true);
        vr = tempLaw->lawVD(0.0, pathOri->getLength()*0.5);
        vsOri = vsOri < vr ? vsOri : vr;
        if (i > 0 && vsOri < trajOri.at(i-1)->getVe() )
            trajOri[i-1]->setVe(vsOri); // Update previous profile.
        trajOri[i] = new Trajectory(vsOri, fmin(veOri, vr),
                                    constOri, pathOri);

    }
    trajList.push_back(trajLin);
    trajList.push_back(trajOri);
    int tCost = timer.elapsed();
    std::cout << "Initialization finished in " << tCost << "ms. \n";

    std::cout << "Planning the motions with SLATP....\n";
    timer.restart();
    for (uint16_t itr = 0; itr < 500; itr++)
        Planner::SLATP(trajList);
    tCost = timer.elapsed();
    std::cout << "Planning finished in " << tCost*0.002 << "ms. \n";

    std::cout << "Planning the motions individually....\n";
    timer.restart();
    Planner::BDS(trajLinOrg);
    Planner::BDS(trajOriOrg);
    tCost = timer.elapsed();
    std::cout << "Planning finished in " << tCost << "ms. \n";

    std::cout << "Storing the trajectory data of SLATP...\n";
    std::ofstream out("../../Planner/trajSLATP.txt");
    double tAll = 0.0;
    double tLin = 0.0;
    double tOri = 0.0;
    uint idxLin = 0;
    uint idxOri = 0;
    double dLin, vLin, aLin, jLin, sLin;
    double dOri, vOri, aOri, jOri, sOri;
    double dAllLin = 0.0;
    double dAllOri = 0.0;
    double ds = 0.0;
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
    
    std::cout << "Storing the trajectory data of traditional method...\n";
    std::ofstream outLin("../../Planner/trajLinear.txt");
    /// Reset the parameters in interpolation.
    tAll = 0.0;
    tLin = 0.0;
    tOri = 0.0;
    idxLin = 0;
    idxOri = 0;
    dAllLin = 0.0;
    dAllOri = 0.0;
    double uc = 0.0; // current parameter.
    if (outLin.is_open() )
    {
        while (1)
        {
            /// Orientation by parameter synchronization.
            while (idxLin < numberTraj &&
                   tLin > trajLinOrg.at(idxLin)->getDuration() )
            {
                dAllLin += trajLinOrg.at(idxLin)->getLength();
                tLin -= trajLinOrg.at(idxLin)->getDuration();
                /*ds = vLin*ts + aLin*ts*ts/2 + jLin*pow(ts, 3)/6;
                ds = ds + dLin - trajLinOrg.at(idxLin)->getLength();
                tLin = ds / trajLinOrg.at(idxLin)->getVe();
                tLin = tLin > 0 ? tLin : 0.0;*/
                idxLin++;
            }
            if (idxLin == numberTraj)
                break;
            dLin = trajLinOrg.at(idxLin)->displacement(tLin);
            vLin = trajLinOrg.at(idxLin)->velocity(tLin);
            aLin = trajLinOrg.at(idxLin)->acceleration(tLin);
            jLin = trajLinOrg.at(idxLin)->jerk(tLin);
            sLin = trajLinOrg.at(idxLin)->snap(tLin);
            uc = trajLinOrg.at(idxLin)->getUc();
            /// Trajectory information of the orientation is not available directly.
            dOri = 0.0;
            vOri = 0.0;
            aOri = 0.0;
            jOri = 0.0;
            sOri = 0.0;
            ptLin = trajLinOrg.at(idxLin)->axialPosition(tLin);
            /// Orientation by parameter synchronization.
            ptOri = trajOriOrg.at(idxLin)->getPath()->calculatePoint(uc);
            coorWCS.head(3) = ptLin;
            coorWCS.tail(3) = ptOri;
            if (0 == tLin)
                coorMCS = vmc.inverseKinematic(coorWCS,
                                               VectorXd::Zero(dim) );
            else
                coorMCS = vmc.inverseKinematic(coorWCS, coorMCS);
            outLin << tAll << " " << dLin+dAllLin << " " << vLin << " "
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
        }
    }
    outLin.close();

    std::ofstream outOri("../../Planner/trajAngular.txt");
    /// Reset the parameters in interpolation.
    tAll = 0.0;
    tLin = 0.0;
    tOri = 0.0;
    idxLin = 0;
    idxOri = 0;
    dAllLin = 0.0;
    dAllOri = 0.0;
    uc = 0.0; // current parameter.
    if (outOri.is_open() )
    {
        while (1)
        {
            /// Position by parameter synchronization.
            while (idxOri < numberTraj &&
                   tOri > trajOriOrg.at(idxOri)->getDuration() )
            {
                dAllOri += RAD2DEG(trajOriOrg.at(idxOri)->getLength() );
                tOri -= trajOriOrg.at(idxOri)->getDuration();
                /*ds = vOri*ts + aOri*ts*ts/2 + jOri*pow(ts, 3)/6;
                ds = ds + dOri - trajOriOrg.at(idxOri)->getLength();
                tOri = ds / trajOriOrg.at(idxOri)->getVe();
                tOri = tOri > 0 ? tOri : 0.0;*/
                idxOri++;
            }
            if (idxOri == numberTraj)
                break;

            /// Trajectory information of the position is not available directly.
            dLin = 0.0;
            vLin = 0.0;
            aLin = 0.0;
            jLin = 0.0;
            sLin = 0.0;

            dOri = trajOriOrg.at(idxOri)->displacement(tOri);
            vOri = trajOriOrg.at(idxOri)->velocity(tOri);
            aOri = trajOriOrg.at(idxOri)->acceleration(tOri);
            jOri = trajOriOrg.at(idxOri)->jerk(tOri);
            sOri = trajOriOrg.at(idxOri)->snap(tOri);
            uc = trajOriOrg.at(idxOri)->getUc();
            ptOri = trajOriOrg.at(idxOri)->axialPosition(tOri);
            /// Position by parameter synchronization.
            ptLin = trajLinOrg.at(idxOri)->getPath()->calculatePoint(uc);
            coorWCS.head(3) = ptLin;
            coorWCS.tail(3) = ptOri;
            if (0 == tLin)
                coorMCS = vmc.inverseKinematic(coorWCS,
                                               VectorXd::Zero(dim) );
            else
                coorMCS = vmc.inverseKinematic(coorWCS, coorMCS);
            outOri << tAll << " " << dLin+dAllLin << " " << vLin << " "
                   << aLin << " " << jLin << " " << sLin << " "
                   << RAD2DEG(dOri)+dAllOri << " " << RAD2DEG(vOri) << " "
                   << RAD2DEG(aOri) << " " << RAD2DEG(jOri) << " "
                   << RAD2DEG(sOri) << " "
                   << coorWCS(0) << " " << coorWCS(1) << " " << coorWCS(2) << " "
                   << coorWCS(3) << " " << coorWCS(4) << " " << coorWCS(5) << " "
                   << coorMCS(0) << " " << coorMCS(1) << " " << coorMCS(2) << " "
                   << coorMCS(3) << " " << coorMCS(4) << std::endl;
            tAll += ts;
            tOri += ts;
        }
    }
    outOri.close();
    

    delete tempLaw;
    /// Deallocate the path memory.
    Path* tp; // Temp path.
    for (uint i=0; i<trajLinOrg.size(); i++)
    {
        tp = trajLinOrg.at(i)->getPath();
        if (tp != nullptr)
        {
            delete tp;
            tp = nullptr;
        }
    }
    /// Deallocate the trajectory memory.
    for (uint i=0; i<trajLin.size(); i++)
    {
        delete trajLin[i];
        trajLin[i] = nullptr;
    }
    for (uint i=0; i<trajOri.size(); i++)
    {
        delete trajOri[i];
        trajOri[i] = nullptr;
    }
    
    for (uint i=0; i<trajLinOrg.size(); i++)
    {
        delete trajLinOrg[i];
        trajLinOrg[i] = nullptr;
    }
    for (uint i=0; i<trajOriOrg.size(); i++)
    {
        delete trajOriOrg[i];
        trajOriOrg[i] = nullptr;
    }

    std::cout << "\nEnd of Program." << std::endl;
    return a.exec();
}
