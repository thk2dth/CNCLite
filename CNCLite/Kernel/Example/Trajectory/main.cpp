#include <QCoreApplication>
#include "../../Trajectory/Trajectory.h"
#include "../../Path/Line.h"
#include <fstream>

using namespace Eigen;
using namespace CNCLite;
using namespace std;

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);
    VectorXd p0(3), p1(3);
    p0 << 6.0, 0.0, 0.0;
    p1 << 0.0, 0.0, 0.0;
    Path *path = new Line(p0, p1);
    VectorXd kineConst(3);
    kineConst << 1.0e2, 5.0e2, 5.0e3;
    Trajectory traj;
    double vs = 40.0;
    double ve = 75.0;
    double Tg = 5.80;
    traj.initialize(vs, ve, kineConst, path);
    traj.plan();
    TimeLaw* law = new JerkBounded();
    law->initialize(kineConst, true);
    /// Reachable velocity.
    double vr = law->lawVD(0.0, 0.5*path->getLength());
    traj.synchronize(Tg );
    double duration = traj.getDuration();
    int num = 1001;
    VectorXd t = VectorXd::LinSpaced(num, 0.0, duration);
    ofstream out("../../Trajectory/traj.txt");
    double u = 0.0;
    double pd = 0.0;
    double dd = 0.0;
    if (out.is_open() )
    {
        for (int i=0; i<num; i++)
        {
            double d = traj.displacement(t(i) );
            double v = traj.velocity(t(i) );
            double a = traj.acceleration(t(i) );
            double j = traj.jerk(t(i) );
            double s = traj.snap(t(i) );
            dd = d - pd; // displacement increase
            pd = d;
            u = traj.getPath()->calculateNextParameter(dd, u, Path::SECONDTE);
            Vector3d pt = traj.getPath()->calculatePoint(u);
            Vector3d vec = traj.getPath()->calculateDer1(u).normalized();
            Vector3d vv = vec * v;
            out << t(i) << " " << d << " " << v << " " << a << " " << j
                << " " << s << " "<< pt(0) << " " << pt(1) << " " << pt(2)
                << " " << vv(0) << " " << vv(1) << " " << vv(2) <<endl;
        }
    }
    out.close();

    delete path;
    path = nullptr;
    delete law;
    law = nullptr;
    cout << "End of program." << endl;
    return a.exec();
}
