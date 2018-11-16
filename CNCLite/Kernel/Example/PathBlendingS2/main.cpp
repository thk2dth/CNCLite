/*****************************************************************************
 * Path blending for tool orientation.
 ***************************************************************************/
#include <QCoreApplication>
#include "../../../Kinematic/Kinematic.h"
#include "../../../Kinematic/AptParser.h"
#include "../../../Planner/Planner.h"
#include "../../Path/Line.h"
#include "../../Path/Arc.h"
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
    timer.start();

    /// Load linear points from a file.
    uint8_t dim = 3; // dimension
    QFile fd("../../../../NC Program/toolOrientation.txt");
    if (!fd.open(QFile::ReadOnly) )
    {
        std::cerr << "Failed to open the file storing the tool orientation.\n";
        return -1;
    }
    QByteArray str = fd.readLine().simplified();
    QList<QByteArray> strList;
    VectorXd temp = VectorXd::Zero(dim);
    std::vector<VectorXd> pts;
    uint32_t numberPoints = 0;

    std::cout << "Begin to load data..." << std::endl;

    while (str != "\0")
    {
        strList = str.split(' '); // Suppose data is delimited by backspace.
        assert(strList.size() <= dim);
        for (uint8_t i=0; i<strList.size(); i++)
            temp(i) = strList.at(i).toDouble();
        pts.push_back(temp);
        numberPoints++;
        str = fd.readLine().simplified();
    }
    fd.close();

    /// Blending parameters.
    double tolerance = DEG2RAD(10); // blending error tolerance, in mm.
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

    for (uint32_t i=0; i<numberPoints-2; i++)
    {
        /// Inserted Bspline.
        ctrlPts[0] = pts.at(i);
        ctrlPts[1] = pts.at(i+1);
        ctrlPts[2] = pts.at(i+2);
        blb = new BsplineAngularBlend(ctrlPts);
        smoothPath[2*i + 1] = blb;
    }

    /// Then insert lines.
    /// The remaining lines and arcs have not been converted to Nurbs curves yet.
    VectorXd ps, pe;
    for (uint32_t i = 0; i <=numberPoints-2; i++)
    {
        /// First line.
        if (0 == i)
        {
            ps = pts.at(0);
            pe = smoothPath.at(1)->getBeginPoint();
            /// Start, end, and center points, respectively.
            smoothPath[0] = new Arc(ps, pe, VectorXd::Zero(3) );
        }
        else if (i == numberPoints - 2)
        {
            ps = smoothPath.at(2*(numberPoints-3) + 1)->getEndPoint();
            pe = pts.at(numberPoints - 1); // Last point.
            smoothPath[2*numberPoints - 4] = new Arc(ps, pe, VectorXd::Zero(3) );
        }
        else
        {
            ps = smoothPath.at(2*i - 1)->getEndPoint();
            pe = smoothPath.at(2*i + 1)->getBeginPoint();
            smoothPath[2*i] = new Arc(ps, pe, VectorXd::Zero(3) );
        }
    }

    /// After path blending, plan the trajectories.
    /// Kinematic constraints.
    uint8_t order = 4; // order of trajectories.
    VectorXd con(order);
    con << DEG2RAD(20.0), DEG2RAD(100.0), DEG2RAD(1000.0), DEG2RAD(10000.0);
    double ts = 4e-3; // sampling period.
    double ce = DEG2RAD(0.01); // chord error for trajectory planning.
    std::deque<Trajectory* > traj(2*numberPoints - 3, nullptr);
    double vs, ve, maxCurv, uMid;


    std::cout << "Begin to initialize the trajectories..." << std::endl;

    for (uint32_t i=0; i<smoothPath.size(); i++)
    {
        if (smoothPath.at(i)->getCurveType() == Path::CURVE_TYPE::ARC)
        {
            if (0 == i)
            {
                vs = 0.0;
                uMid = 0.5 * (smoothPath.at(i+1)->getBeginParameter() +
                              smoothPath.at(i+1)->getEndParameter() );
                maxCurv = smoothPath.at(i+1)->calculateCurvature(uMid);
                ve = Path::calculateTraverseVelocityS2(maxCurv, con(0), con(1),
                                                                   con(2), ce, ts);
                traj[i] = new Trajectory(vs, ve, con, smoothPath[i] );
            }
            else if (i == smoothPath.size() - 1)
            {
                ve = 0.0;
                uMid = 0.5 * (smoothPath.at(i-1)->getBeginParameter() +
                              smoothPath.at(i-1)->getEndParameter() );
                maxCurv = smoothPath.at(i-1)->calculateCurvature(uMid);
                vs = Path::calculateTraverseVelocityS2(maxCurv, con(0), con(1),
                                                                   con(2), ce, ts);
                traj[i] = new Trajectory(vs, ve, con, smoothPath[i] );
            }
            else
            {
                uMid = 0.5 * (smoothPath.at(i-1)->getBeginParameter() +
                              smoothPath.at(i-1)->getEndParameter() );
                maxCurv = smoothPath.at(i-1)->calculateCurvature(uMid);
                vs = Path::calculateTraverseVelocityS2(maxCurv, con(0), con(1),
                                                                   con(2), ce, ts);
                uMid = 0.5 * (smoothPath.at(i+1)->getBeginParameter() +
                              smoothPath.at(i+1)->getEndParameter() );
                maxCurv = smoothPath.at(i+1)->calculateCurvature(uMid);
                ve = Path::calculateTraverseVelocityS2(maxCurv, con(0), con(1),
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
            VectorXd tempCon(order);
            vs = Path::calculateTraverseVelocityS2(maxCurv, con(0), con(1),
                                                               con(2), ce, ts);
            tempCon << vs, 0.0, 0.0, 0.0; // Constant velocity.
            traj[i] = new Trajectory(vs, vs, tempCon, smoothPath[i] );
        }
    }

    std::cout << "Begin to plan the trajectories..." << std::endl;
    Planner::BDS(traj);

    /// Store the trajectory data
    std::ofstream out("../../PathBlendingS2/trajS2.txt");
    double tAll = 0.0;
    double tOri = 0.0;
    uint32_t idx = 0;
    double dOri, vOri, aOri, jOri, sOri;
    double vOriPre, aOriPre, jOriPre, sOriPre;
    double dAllOri = 0.0;
    VectorXd ptOri;
    double ds = 0.0;
    double dt = ts;
    double uc = 0.0;
    double cur = 0.0;

    std::cout << "Begin to store the trajectory data..." << std::endl;

    if (out.is_open() )
    {
        while (1)
        {
            while (idx < traj.size() &&
                   tOri > traj.at(idx)->getDuration() )
            {
                dAllOri += traj.at(idx)->getLength();
                tOri -= traj.at(idx)->getDuration();
                dt = tOri;
                vOriPre = traj.at(idx)->getVe();
                aOriPre = traj.at(idx)->acceleration(traj.at(idx)->getDuration() );
                jOriPre = traj.at(idx)->jerk(traj.at(idx)->getDuration() );
                sOriPre = traj.at(idx)->snap(traj.at(idx)->getDuration() );
                uc = traj.at(idx)->getUe();
                dOri = 0.0;
                idx++;
                out << std::endl;
            }
            if (idx == traj.size() )
                break;
            dOri = traj.at(idx)->displacement(tOri);
            vOri = traj.at(idx)->velocity(tOri);
            aOri = traj.at(idx)->acceleration(tOri);
            jOri = traj.at(idx)->jerk(tOri);
            sOri = traj.at(idx)->snap(tOri); // snap
            ds = vOriPre*ts + aOriPre*ts*ts/2 + jOriPre*ts*ts*ts/6;
            uc = traj.at(idx)->getPath()->calculateNextParameter(ds, uc, Path::SECONDRUKUCOM);
            cur = traj.at(idx)->getPath()->calculateCurvature(uc);
            ptOri = traj.at(idx)->axialPosition(tOri); // axial position
            out << tAll << " " << RAD2DEG(dOri + dAllOri) << " " << RAD2DEG(vOri) << " "
                << RAD2DEG(aOri) << " " << RAD2DEG(jOri) << " " << RAD2DEG(sOri) << " " << ptOri(0)
                << " " << ptOri(1) << " " << ptOri(2) << " " << cur << std::endl;
            tAll += ts;
            tOri += ts;
            vOriPre = vOri;
            aOriPre = aOri;
            jOriPre = jOri;
            sOriPre = sOri;
            dt = ts;
        }
    }
    out.close();

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

    std::cout << "\nEnd of Program." << std::endl;
    return a.exec();
}
