#include <QCoreApplication>
#include "../../Path/Arc.h"
#include "../../Path/Line.h"
#include "../../Path/Nurbs.h"
#include "../../Path/BsplineAngularBlend.h"
#include <fstream>

using namespace Eigen;
using namespace std;
using namespace CNCLite;

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);
    int num = 1000;
    VectorXd u = VectorXd::LinSpaced(num, 0.0, 1.0);
    VectorXd bgnPt(3), endPt(3), ctPt(3);
    bgnPt << 1.0, 0.0, 0.0;
    endPt << 0.0, 1.0, 0.0;
    ctPt << 0.0, 0.0, 0.0;
    Arc arc(bgnPt, endPt, ctPt);
    Line line(bgnPt, endPt);
    VectorXd pt(3);
    cout << "Arc: " << endl;
    for (int i=0; i<num; i++)
    {
        pt = arc.calculatePoint( u(i) );
        cout << "The " << i << "th point is: "<< endl << pt <<";" << endl;
    }
    cout << "Line: " << endl;
    for (int i=0; i<num; i++)
    {
        pt = line.calculatePoint( u(i) );
        cout << "The " << i << "th point is: "<< endl << pt <<";" << endl;
    }

    vector<VectorXd> ctrlPts;
    ctrlPts.push_back(Vector4d(0.0, 0.0, 0.0, 1.0) );
    ctrlPts.push_back(Vector4d(0.0, 1.0, 0.0, 1.0) );
    ctrlPts.push_back(Vector4d(1.0, 1.0, 0.0, 1.0) );
    ctrlPts.push_back(Vector4d(1.0, 0.0, 0.0, 1.0) );
    vector<double> knotVec(8);
    unsigned int degree = 3;
    for (int i=0; i<=3; i++)
    {
        knotVec[i] = 0.0;
        knotVec[i+degree+1] = 1.0;
    }
    Nurbs nrb(ctrlPts, knotVec, true);
    for (int i=0; i<num; i++)
    {
        pt = nrb.calculatePoint( u(i) );
        cout << "The " << i << "th point is: "<< endl << pt <<";" << endl;
    }

    Nurbs nrbJson;
    nrbJson.initializeFromJson("../../../../NC Program/NurbsData/Fan.json");
    ofstream contour("../../../../Results/Kernel/Path/contour.txt");
    if (contour.is_open() )
    {
        for (int i=0; i<num; i++)
        {
            contour << u(i) << ' ';
            pt = nrbJson.calculatePoint( u(i) );
            for (int j=0; j<pt.size(); j++)
                contour << pt(j) << ' ';
            contour << endl;
        }
    }
    contour.close();

    BsplineAngularBlend bab;
    ctrlPts.clear();
    ctrlPts.push_back(Vector3d(1.0, 0.0, 0.0) );
    ctrlPts.push_back(Vector3d(0.0, 1.0, 0.0) );
    ctrlPts.push_back(Vector3d(0.0, 0.0, 1.0) );
    ctrlPts.push_back(Vector3d(1.0/3, 0.25, 5) ); // blending error, 5 degrees.
    bab.initializeFromPoints(ctrlPts);

    ofstream spc("../../../../Results/Kernel/Path/sphericalCurve.txt");
    double cur;
    if (spc.is_open() )
    {
        for (int i=0; i<num; i++)
        {
            spc << u(i) << ' ';
            pt = bab.calculatePoint( u(i) );
            cur = bab.calculateCurvature( u(i) );
            for (int j=0; j<pt.size(); j++)
                spc << pt(j) << ' ';
            spc << cur << endl;
        }
    }
    spc.close();


    cout<< "End of program!" <<endl;

    return a.exec();
}
