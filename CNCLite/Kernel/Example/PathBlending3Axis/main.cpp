/**********************************************************************
 * Example PathBlending3Axis demonstrates the path blending for 3-axis
 * linear toolpath.
 * The path is firstly blended with specified tolerance.
 * Then, the trajectory of the blended path is planned. The inserted
 * parametric curves are traversed with constant velocities.
 * HJ, 20180414.
**********************************************************************/
#include <QCoreApplication>
#include <QTime>
#include <QByteArray>
#include <QList>
#include <QString>
#include <QFile>
#include <vector>
#include <deque>
#include <fstream>
#include "../../Path/Nurbs.h"
#include "../../Path/BsplineLinearBlend.h"
#include "../../Trajectory/Trajectory.h"
#include "../../Planner/Planner.h"

using namespace CNCLite;
using namespace Eigen;

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);
    QTime timer;
    /// Load linear points from a file.
    uint8_t dim = 3; // dimension
    /// QFile fd("../../../../NC Program/Test.dat");
    QFile fd("../../../../NC Program/OriginalSpiralDataTest.dat");
    if (!fd.open(QFile::ReadOnly) )
    {
        std::cerr << "Failed to open the file storing the linear toolpath.\n";
        return -1;
    }
    QByteArray str = fd.readLine().simplified();
    QList<QByteArray> strList;
    VectorXd temp = VectorXd::Zero(dim);
    std::vector<VectorXd> pts;
    uint32_t numberPoints = 0;
    uint8_t interval = 1;
    uint32_t numLine = 0;

    std::cout << "Begin to load data..." << std::endl;

    while (str != "\0")
    {
        strList = str.split(' '); // Suppose data is delimited by backspace.
        assert(strList.size() <= dim);
        for (uint8_t i=0; i<strList.size(); i++)
            temp(i) = strList.at(i).toDouble();
        if (numLine % interval == 0)
        {
            pts.push_back(temp);
            numberPoints++;
        }
        str = fd.readLine().simplified();
        numLine++;
    }
    fd.close();

    /// Blending parameters.
    double tolerance = 0.05; // blending error tolerance, in mm.
    double cRatio = 0.25;
    double mRatio = 1 / 3.0;

    VectorXd param(dim);
    param << mRatio, cRatio, tolerance;

    /// Blend the line segments by inserting Bsplines.
    /// After blending, number of path segments is 2*numberPoints-3.
    std::vector<Path* > smoothPath(2*numberPoints - 3, nullptr);
    Path* blb = nullptr; // Inserted Bspline.
    std::vector<VectorXd> ctrlPts(4, VectorXd::Zero(dim) ); // Blend the corner forming by three points.
    ctrlPts[3] = param;

    std::cout << "Begin to blend the path..." << std::endl;
    timer.start();
    for (uint32_t i=0; i<numberPoints-2; i++)
    {
        /// Inserted Bspline.
        ctrlPts[0] = pts.at(i);
        ctrlPts[1] = pts.at(i+1);
        ctrlPts[2] = pts.at(i+2);
        blb = new BsplineLinearBlend(ctrlPts);
        smoothPath[2*i + 1] = blb;
    }

    /// Then insert lines.
    /// Finally, convert the remaining line segments to cubic Bsplines by inserting
    /// two intermediate points.
    std::vector<VectorXd> linePts(4, VectorXd::Zero(dim)); // Four points for remaining line segments.
    Vector3d lineSlope; // slope of the remaining line.
    double ne, nb; // norms of begin and end tangents.

    for (uint32_t i=0; i<=numberPoints-2; i++)
    {
        /// First line.
        if (0 == i)
        {
            linePts[0] = pts.at(0);
            linePts[3] = smoothPath.at(1)->getBeginPoint();
            lineSlope = linePts.at(3) - linePts.at(0);
            ne = smoothPath.at(1)->getBeginTangent().norm();
            linePts[2] = linePts.at(3) - lineSlope * (ne / (lineSlope.norm() * 3.0) );
            linePts[1] = (linePts.at(0) + linePts.at(2) ) * 0.5; // V1 is not constrained.
            /// By default, the knot vector is [0,0,0,0,1,1,1,1].
            smoothPath[0] = new Nurbs(linePts);
        }
        else if (i == numberPoints-2)
        {
            linePts[0] = smoothPath.at(2 * (numberPoints - 3) + 1)->getEndPoint();
            linePts[3] = pts.at(numberPoints - 1);
            lineSlope = linePts.at(3) - linePts.at(0);
            nb = smoothPath.at(2 * (numberPoints - 3) )->getEndTangent().norm();
            linePts[1] = linePts.at(0) + lineSlope * (nb / (lineSlope.norm() * 3.0) );
            linePts[2] = (linePts.at(1) + linePts.at(3) ) * 0.5;
            smoothPath[2*numberPoints-4] = new Nurbs(linePts);
        }
        else
        {
            linePts[0] = smoothPath.at(2*i - 1)->getEndPoint();
            linePts[3] = smoothPath.at(2*i + 1)->getBeginPoint();
            lineSlope = linePts.at(3) - linePts.at(0);
            nb = smoothPath.at(2*i - 1)->getEndTangent().norm();
            ne = smoothPath.at(2*i + 1)->getBeginTangent().norm();
            linePts[1] = linePts.at(0) + lineSlope * (nb / (lineSlope.norm() * 3.0) );
            linePts[2] = linePts.at(3) - lineSlope * (ne / (lineSlope.norm() * 3.0) );
            smoothPath[2*i] = new Nurbs(linePts);
        }
    }

    /// After path blending, plan the trajectories.
    /// Kinematic constraints.
    uint8_t order = 4; // order of trajectories.
    VectorXd con(order);
    con << 30.0, 1000.0, 1e4, 1e5;
    double ts = 1e-2; // sampling period.
    double ce = 0.001; // chord error for trajectory planning.
    std::deque<Trajectory* > traj(2*numberPoints - 3, nullptr);
    double vs, ve, maxCurv, uMid;

    std::cout << "Begin to initialize the trajectories..." << std::endl;

    for (uint32_t i=0; i<smoothPath.size(); i++)
    {
        if (smoothPath.at(i)->getCurveType() == Path::CURVE_TYPE::NURBS)
        {
            if (0 == i)
            {
                vs = 0.0;
                uMid = 0.5 * (smoothPath.at(i+1)->getBeginParameter() +
                              smoothPath.at(i+1)->getEndParameter() );
                maxCurv = smoothPath.at(i+1)->calculateCurvature(uMid);
                ve = Path::calculateTraverseVelocityR3(maxCurv, con(0), con(1),
                                                       con(2), ce, ts);
                traj[i] = new Trajectory(vs, ve, con, smoothPath[i] );
            }
            else if (i == smoothPath.size() - 1)
            {
                ve = 0.0;
                uMid = 0.5 * (smoothPath.at(i-1)->getBeginParameter() +
                              smoothPath.at(i-1)->getEndParameter() );
                maxCurv = smoothPath.at(i-1)->calculateCurvature(uMid);
                vs = Path::calculateTraverseVelocityR3(maxCurv, con(0), con(1),
                                                       con(2), ce, ts);
                traj[i] = new Trajectory(vs, ve, con, smoothPath[i] );
            }
            else
            {
                uMid = 0.5 * (smoothPath.at(i-1)->getBeginParameter() +
                              smoothPath.at(i-1)->getEndParameter() );
                maxCurv = smoothPath.at(i-1)->calculateCurvature(uMid);
                vs = Path::calculateTraverseVelocityR3(maxCurv, con(0), con(1),
                                                       con(2), ce, ts);
                uMid = 0.5 * (smoothPath.at(i+1)->getBeginParameter() +
                              smoothPath.at(i+1)->getEndParameter() );
                maxCurv = smoothPath.at(i+1)->calculateCurvature(uMid);
                ve = Path::calculateTraverseVelocityR3(maxCurv, con(0), con(1),
                                                       con(2), ce, ts);
                traj[i] = new Trajectory(vs, ve, con, smoothPath[i] );
            }
        }
        /// Use a trajectory with constant velocity for the inserted Bspline.
        else
        {
            uMid = 0.5 * (smoothPath.at(i)->getBeginParameter() +
                          smoothPath.at(i)->getEndParameter() );
            maxCurv = smoothPath.at(i)->calculateCurvature(uMid);
            VectorXd tempCon = con;
            vs = Path::calculateTraverseVelocityR3(maxCurv, con(0), con(1),
                                                   con(2), ce, ts);
            tempCon(0) = vs; // velocity less than vs.
            traj[i] = new Trajectory(vs, vs, tempCon, smoothPath[i] );
        }
    }
    std::cout << "Complete the trajectory initialization in " << timer.elapsed() << "ms." << std::endl;

    std::cout << "Begin to plan the trajectories..." << std::endl;
    timer.restart();
    for (uint16_t itrNum = 0; itrNum < 500; itrNum++)
        Planner::BDS(traj);
    std::cout << "Complete the feedrate scheduling in " << timer.elapsed()*0.002 << "ms." << std::endl;

    /// Store the trajectory data
    std::ofstream out("../../PathBlending3Axis/traj3Axis.txt");
    std::ofstream gcode("../../PathBlending3Axis/gCodeTransition.MPF");
    double tAll = 0.0;
    double tLin = 0.0;
    uint32_t idx = 0;
    double dLin, vLin, aLin, jLin, sLin;
    double dAllLin = 0.0;
    VectorXd ptLin, vtLin;
    double uc = 0.0;
    double cur = 0.0;
    std::cout << "Begin to store the trajectory data..." << std::endl;

    if (out.is_open() && gcode.is_open() )
    {
        while (1)
        {
            while (idx < traj.size() &&
                   tLin > traj.at(idx)->getDuration() )
            {
                dAllLin += traj.at(idx)->getLength();
                /*ds = vLin*ts + aLin*ts*ts/2 + jLin*pow(ts, 3) / 6
                        + sLin*pow(ts, 4) / 24;
                ds = ds + dLin -  traj.at(idx)->getLength();
                tLin = ds / traj.at(idx)->getVe();
                tLin = tLin > 0 ? tLin : 0.0;*/
                tLin -= traj.at(idx)->getDuration();
                idx++;
                out << std::endl;
            }
            if (idx == traj.size() )
                break;
            dLin = traj.at(idx)->displacement(tLin);
            vLin = traj.at(idx)->velocity(tLin);
            aLin = traj.at(idx)->acceleration(tLin);
            jLin = traj.at(idx)->jerk(tLin);
            sLin = traj.at(idx)->snap(tLin); // snap
            ptLin = traj.at(idx)->axialPosition(tLin); // axial position
            vtLin = traj.at(idx)->axialVelocity(tLin); // axial velocity
            uc = traj.at(idx)->getUc();
            cur = traj.at(idx)->getPath()->calculateCurvature(uc);
            out << tAll << " " << dLin + dAllLin << " " << vLin << " "
                << aLin << " " << jLin << " " << sLin << " " << ptLin(0)
                << " " << ptLin(1) << " " << vtLin(0) << " " << vtLin(1)
                << " " << uc << " " << cur << std::endl;
            gcode << "G01 " << "X" << ptLin(0) << " Y" << ptLin(1) << " Z" <<
                     ptLin(2) << " F" << vLin*60 <<std::endl;
            tAll += ts;
            tLin += ts;
        }
    }
    out.close();
    gcode.close();

    /// Deallocate the path memory.
    Path* tp;
    for (uint32_t i=0; i<traj.size(); i++)
    {
        tp = traj.at(i)->getPath();
        if (tp != nullptr)
        {
            delete tp;
            tp = nullptr;
        }
    }

    /// Deallocate the trajectory memory.
    for (uint32_t i=0; i<traj.size(); i++)
    {
        if (traj.at(i) != nullptr)
        {
            delete traj[i];
            traj[i] = nullptr;
        }
    }


    std::cout << "Elapsed time for this program is: " << timer.elapsed() * 0.001 << std::endl;
    std::cout << "\nEnd of Program." << std::endl;
    return a.exec();
}
