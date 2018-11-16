#include "AptParser.h"

using namespace CNCLite;
using namespace Eigen;

AptParser::AptParser(const QString& fn)
{
    aptFile.setFileName(fn);
    isFileValid = aptFile.open(QIODevice::ReadOnly);
    /// Valid data in apt file has the following pattern, eg, -123.456
    reFormat = QRegExp("(-?\\d*\\.?\\d+)");
}

AptParser::~AptParser()
{
    if (isFileValid)
    {
        aptFile.close();
        isFileValid = false;
    }
}

VectorXd AptParser::readLine()
{
    if (!isFileValid || aptFile.atEnd() )
        return VectorXd::Zero(5);
    std::vector<double> v;
    int idx = 0;
    QString line = aptFile.readLine();
    while ((idx = reFormat.indexIn(line, idx)) != -1)
    {
        /// There are remaining data in the line.
        v.push_back(reFormat.cap(1).toDouble() );
        idx += reFormat.matchedLength();
    }
    VectorXd temp(v.size() );
    for (uint8_t i=0; i<v.size(); i++)
        temp(i) = v.at(i);
    Vector3d Ori = temp.tail(3);
    Ori.normalize(); // normalize the orientation vector.
    temp.tail(3) = Ori;
    return temp;
}

void AptParser::readAll(MatrixXd& data, uint interval)
{
    if (!isFileValid || aptFile.atEnd() )
        return;
    uint i = 0;
    std::vector<VectorXd > tempData;
    VectorXd t;
    while (!aptFile.atEnd() )
    {
        t = readLine();
        if (i % interval == 0)
        {
           tempData.push_back(t);
        }
        i++;
    }
    uint num = tempData.size();
    uint dim = t.size();
    /// Resize the matrix according to the number and dimension of the data.
    data.resize(num, dim);
    for (i=0; i<num; i++)
        data.row(i) = tempData.at(i);
}
