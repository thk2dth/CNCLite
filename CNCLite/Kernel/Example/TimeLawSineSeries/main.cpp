/************************************************************
 * Debug the time law with sine series representation.
***********************************************************/
#include <QCoreApplication>
#include "../../TimeLaw/JerkSineSeries.h"
#include <fstream>

using namespace Eigen;
using namespace CNCLite;
using namespace std;

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);
    double vs = 10.0;
    double ve = 100.0;
    unsigned int order = 4;
    VectorXd kineConst(order);
    kineConst << 1e2, 1e3, 1e4, 1e5;
    JerkSineSeries traj;
    traj.initialize(kineConst);
    traj.lawVV(vs, ve);
    double duration = traj.getDuration();
    int num = 1001;
    VectorXd t = VectorXd::LinSpaced(num, 0.0, duration);
    ofstream out("../../TimeLawSineSeries/traj.txt");
    if ( out.is_open() )
    {
        for (int i=0; i<num; i++)
        {
            double s = traj.snap( t(i) );
            double j = traj.jerk( t(i) );
            double a = traj.acceleration( t(i) );
            double v = traj.velocity( t(i) );
            double d = 0.0;
            out << t(i) <<" " << d <<" " << v <<" " << a <<" " << j <<" " << s <<endl;
        }
        out.close();
    }
    cout << "End of program." << endl;
    return a.exec();
}
