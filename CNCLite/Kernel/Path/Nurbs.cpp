#include "Nurbs.h"

/// Use QJson to parse json data.
#include <QString>
#include <QFile>
#include <QTextStream>
#include <QByteArray>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QList>



using namespace CNCLite;
using namespace Eigen;

Nurbs::Nurbs()
{
    curveType = NURBS;
    isInitialized = false;
    isHomogeneous = false;
}

Nurbs::~Nurbs()
{

}

Nurbs::Nurbs(const std::vector<VectorXd> &ctrlPts,
             double bgnPara, double endPara)
{
    initializeFromPoints(ctrlPts, bgnPara, endPara);
}

Nurbs::Nurbs(const std::vector<VectorXd> &ctrlPts,
             const std::vector<double> &knotVec, bool isHomo)
{
    initializeFromPoints(ctrlPts, knotVec, isHomo);
}

void Nurbs::initializeFromPoints(const std::vector<VectorXd> &ctrlPts,
                                 double bgnPara, double endPara)
{
    degree = ctrlPts.size() - 1;
    knotVector.resize(degree*2 + 2);
    beginParameter = bgnPara;
    endParameter = endPara;
    for (unsigned int i=0; i<=degree; i++)
    {
        knotVector[i] = beginParameter;
        knotVector[i+degree+1] = endParameter;
    }
    isHomogeneous = false;
    controlPoints = ctrlPts;
    numberControlPoints = ctrlPts.size();
    dimension = controlPoints.front().size();
    curveType = NURBS;
    isInitialized = true;
    initialize();
}

void Nurbs::initializeFromPoints(const std::vector<VectorXd> &ctrlPts,
                                 const std::vector<double> &knotVec, bool isHomo)
{
    knotVector = knotVec;
    controlPoints = ctrlPts;
    numberControlPoints = ctrlPts.size();
    beginParameter = knotVector.front();
    endParameter = knotVector.back();
    isHomogeneous = isHomo;
    degree = knotVector.size() - ctrlPts.size() - 1;
    dimension = ctrlPts.front().size();
    if (isHomogeneous)
        dimension -= 1; /// If homogeneous coordinate, the last element is weight.
    curveType = NURBS;
    isInitialized = true;
    initialize();
}

void Nurbs::initializeFromFiles(const std::string &ctrlPtsFileName,
                                const std::string &knotVecFileName, bool isHomo, double r)
{
    controlPoints.clear();
    knotVector.clear();
    numberControlPoints = 0;
    QFile fc( ctrlPtsFileName.data() );
    QFile fk( knotVecFileName.data() );
    if (!fc.open(QFile::ReadOnly) )
    {
        std::cerr << "Failed to open the file that stores the control points." << std::endl;
        return;
    }
    if (!fk.open(QFile::ReadOnly) )
    {
        std::cerr << "Failed to open the file that stores the knot vector." << std::endl;
        return;
    }
    QByteArray str = fc.readLine().simplified();
    QList<QByteArray> strList = str.split(' '); // Suppose data is delimited by backspace.
    if (isHomo)
        dimension = strList.size() - 1;
    else
        dimension = strList.size();
    VectorXd temp(strList.size() );
    double w = 1.0; // weight.
    while (str != "\0")
    {
        if (isHomo)
        {
            w = strList.back().toDouble();
            for (int i=0; i<dimension; i++)
            {
                /// Homogeneous point should formatted as [w*x, w*y, w*z, w].
                temp(i) = strList.at(i).toDouble() * w * r;
            }
            temp(dimension) = w;
        }
        else
        {
            for (int i=0; i<dimension; i++)
                temp(i) = strList.at(i).toDouble() * r;
        }
        controlPoints.push_back(temp);
        numberControlPoints++;
        str = fc.readLine().simplified();
    }
    fc.close();

    str = fk.readLine().simplified();
    while (str != "\0")
    {
        knotVector.push_back(str.toDouble() );
        str = fk.readLine().simplified();
    }
    fk.close();
    beginParameter = knotVector.front();
    endParameter = knotVector.back();
    isHomogeneous = isHomo;
    degree = knotVector.size() - controlPoints.size() - 1;
    curveType = NURBS;
    isInitialized = true;
    initialize();
}

void Nurbs::initializeFromJson(const std::string &fileName, double r, bool isValidate)
{
    QFile fj(fileName.data() );
    if (!fj.open(QFile::ReadOnly) )
    {
        std::cerr << "Failed to open json file." << std::endl;
        return;
    }
    else
    {
        QByteArray data = fj.readAll();
        fj.close();
        QJsonDocument jsonDoc(QJsonDocument::fromJson(data) );
        QJsonObject jsonObj = jsonDoc.object();
        degree = jsonObj["Degree"].toInt();
        dimension = jsonObj["Dimension"].toInt();
        numberControlPoints = jsonObj["NumberControlPoints"].toInt();
        isHomogeneous = jsonObj["IsHomogeneous"].toBool();
        unsigned int size = (isHomogeneous ? dimension+1 : dimension);
        VectorXd temp(size);
        /// json object stores control points and knot vector.
        /// All the data is in an array and organized as [x, y, z, w]
        QJsonArray pts = jsonObj["ControlPoints"].toArray();
        QJsonArray knotVec = jsonObj["KnotVector"].toArray();
        /// Data validation
        if (isValidate)
        {
            if (!((knotVec.size() == numberControlPoints+degree+1)
                  && (pts.size() == size*numberControlPoints)) )
            {
                std::cerr << "Json data validation failed." << std::endl;
                return;
            }
        }
        for (int i=0; i<numberControlPoints; i++)
        {
            if (isHomogeneous)
            {
                double w = pts.at( (size-1) * numberControlPoints + i).toDouble();
                for (int j=0; j<size-1; j++)
                {
                    temp(j) = pts.at(j*numberControlPoints + i).toDouble() * w * r;
                }
                temp(size-1) = w;
            }
            else
            {
                for (int j=0; j<size; j++)
                {
                    temp(j) = pts.at(j*numberControlPoints + i).toDouble() * r;
                }
            }
            controlPoints.push_back(temp);
        }
        for (int i=0; i<numberControlPoints+degree+1; i++)
            knotVector.push_back(knotVec.at(i).toDouble() );
        beginParameter = knotVector.front();
        endParameter = knotVector.back();
        curveType = NURBS;
        isInitialized = true;
        initialize();
    }
}

bool Nurbs::toJsonFile(const std::string &fileName) const
{
    return false;
}

VectorXd Nurbs::calculatePoint(double u) const
{
    if (fabs(u-endParameter) <= EPS_NUM)
        u = endParameter;
    if (fabs(u-beginParameter) <= EPS_NUM)
        u = beginParameter;
    if (u - beginParameter < -EPS_NUM)
    {
        std::cout << "The given parameter less than the begin paramter. "
                     "The begin parameter is used."<<std::endl;
        u = beginParameter;
    }
    if (u - endParameter > EPS_NUM)
    {
        std::cout<<"The given parameter is great than the end parameter. "
                   "The end parameter is used."<<std::endl;
        u = endParameter;
    }
    int uIndex, i, j;
    std::vector<std::vector<VectorXd> > iterativeMat(degree + 1);
    for (i=0; i<=degree; i++)
        iterativeMat[i].resize(degree+1);
    uIndex = findSpan(u);
    MatrixXd alphaMat = alphaMatrix(u, uIndex);
    for (i=degree; i>=0; i--)
    {
        for (j=0; j<=i; j++)
        {
            if (i == degree)
                iterativeMat[j][i] = tempIterativeN(uIndex-i+j, 0);
            else
                iterativeMat[j][i] = iterativeMat[j][i+1]*(1-alphaMat(j,i) ) +
                        iterativeMat[j+1][i+1]*alphaMat(j,i);
        }
    }
    if (isHomogeneous)
        return affineHomogenous(iterativeMat[0][0]);
    else
        return iterativeMat[0][0];
}

VectorXd Nurbs::calculateDer1(double u) const
{
    std::vector<VectorXd> rst(2);
    calculateDerN(u, 1, rst);
    return rst.at(1);
}

VectorXd Nurbs::calculateDer2(double u) const
{
    std::vector<VectorXd> rst(3);
    calculateDerN(u, 2, rst);
    return rst.at(2);
}

void Nurbs::calculateDerN(double u, unsigned int order, std::vector<VectorXd> &rst) const
{
    int uIndex, i, j, k;
    std::vector<std::vector<VectorXd> > iterativeMat(degree+1);
    std::vector<VectorXd> rstTemp(order+1);
    for (i=0; i<=degree; i++)
        iterativeMat[i].resize(degree+1);
    uIndex = findSpan(u);
    MatrixXd alphaMat = alphaMatrix(u, uIndex);
    VectorXd vec;
    if (isHomogeneous)
    {
        vec = VectorXd::Zero(dimension+1);
        vec(dimension) = 1.0; // A zero homogenous coordinate.
    }
    else
        vec = VectorXd::Zero(dimension);
    for (k=0; k<=order; k++)
    {
        if (k < degree+1)
        {
            for (i=degree-k; i>=0; i--)
                for (j=0; j<=i; j++)
                {
                    if (i == degree-k)
                        iterativeMat[j][i] = tempIterativeN(uIndex-i+j, k);
                    else
                        iterativeMat[j][i] = iterativeMat[j][i+1] *(1-alphaMat(j, i) ) +
                                iterativeMat[j+1][i+1]*alphaMat(j, i);
                }
            rstTemp[k] = iterativeMat[0][0];
        }
        else
            rstTemp[k] = vec;
    }
    /// C=A/W, C'=(A'-W'*C)/W, etc. The NURBS BOOK Eq. (4.6).
    if (isHomogeneous)
    {
        std::vector<VectorXd> tempC(order+1);
        rst[0] = affineHomogenous( rstTemp[0] );
        tempC[0] = rst[0];
        if (order >= 1)
        {
            rst[1] = (rstTemp[1].head(dimension) - tempC[0]*rstTemp[1](dimension) ) / rstTemp[0](dimension);
            tempC[1] = rst[1];
        }
        if (order >= 2)
        {
            rst[2] = (rstTemp[2].head(dimension) -tempC[1]*(2*rstTemp[1](dimension) ) -
                    tempC[0]*rstTemp[2](dimension) ) / rstTemp[0](dimension);
            tempC[2] = rst[2];
        }
        if (order >= 3)
        {
            rst[3] = (rstTemp[3].head(dimension) - tempC[1]*(3*rstTemp[2](dimension)) -
                    tempC[2]*(3*rstTemp[1](dimension)) - tempC[0]*rstTemp[3](dimension) ) / rstTemp[0](dimension);
            tempC[3] = rst[3];
        }
        if (order >= 4)
        {
            rst[4] = (rstTemp[4].head(dimension) - tempC[1]*(4*rstTemp[3](dimension)) -
                    tempC[2]*(6*rstTemp[2](dimension)) - tempC[3]*(4*rstTemp[1](dimension)) -
                    tempC[0]*rstTemp[4](dimension) ) / rstTemp[0](dimension);
            tempC[4] = rst[4];
        }
        if (order >= 5)
        {
            rst[5] = (rstTemp[5].head(dimension) - tempC[1]*(5*rstTemp[4](dimension)) -
                    tempC[2]*(10*rstTemp[3](dimension)) - tempC[3]*(10*rstTemp[2](dimension)) -
                    tempC[4]*(5*rstTemp[1](dimension) ) - tempC[0]*rstTemp[5](dimension) ) / rstTemp[0](dimension);
            tempC[5] = rst[5];
        }
        if (order >= 6)
        {
            std::cout << "The Sixth and higher orders are not available yet!\n"<<std::endl;
            rst[6] = VectorXd::Zero(dimension);
        }
    }
    else
    {
        /// For B-spline, the derivatives of W are always 0.
        rst = rstTemp;
    }
}

double Nurbs::calculateCurvature(double u) const
{
    return calculateCurvatureByDerivatives(u);
}

double Nurbs::calculateLengthInterval(double us, double ue) const
{
    double len = 0.0;
    unsigned int count = 0;
    calculateLengthIntervalSimpson(us, ue, &len, &count);
    return len;
}

double Nurbs::calculateNextParameter(double ds, double us, const INTERPOLATION_METHOD method) const
{
    double ue = beginParameter;
    switch (method) {
    case FIRSTTE:
        ue = calculateNextParameterFirstTaylor(ds, us);
        break;
    case SECONDTE:
        ue = calculateNextParameterSecondTaylor(ds, us);
        break;
    case NEWTON:
        ue = calculateNextParameterNewton(ds, us);
        break;
    case SECONDRUKUCOM:
        ue = calculateNextParameterSecondRuKuCompensation(ds, us);
        break;
    default:
        std::cerr << "Undefined interpolation method: " << method << std::endl;
        break;
    }
    return ue;
}

void Nurbs::scaleByParemeter(double r)
{
    double k = r / (endParameter - beginParameter);
    for (unsigned int i=1; i<=numberControlPoints+degree; i++)
        knotVector[i] = beginParameter + k * (knotVector[i]-beginParameter);
    beginTangent /= r;
    endTangent /= r;
}

void Nurbs::calculateBeginPoint()
{
    if (!isHomogeneous)
        beginPoint = controlPoints.front();
    else
        /// [w*x, w*y, w*z] / w
        beginPoint = affineHomogenous(controlPoints.front() );
}

void Nurbs::calculateEndPoint()
{
    if (!isHomogeneous)
        endPoint = controlPoints.back();
    else
        endPoint = affineHomogenous(controlPoints.back() );
}

/// The NURBS BOOK Eqs. (4.9, 4.10)
void Nurbs::calculateBeginTangent()
{
    if (!isHomogeneous)
        beginTangent = (controlPoints.at(1) - controlPoints.at(0) ) *
                (degree / (knotVector.at(degree+1)-knotVector.at(degree) ) );
    else
        beginTangent = (affineHomogenous(controlPoints.at(1) ) - affineHomogenous(controlPoints.at(0) ) ) *
                (degree*controlPoints.at(1)(dimension) / ( (knotVector.at(degree+1)-knotVector.at(degree))
                                                           * controlPoints.at(0)(dimension) ) );
}

void Nurbs::calculateEndTangent()
{
    if (!isHomogeneous)
        endTangent = (controlPoints.back() - controlPoints.at(numberControlPoints-2) ) *
                (degree / (knotVector.at(numberControlPoints)-knotVector.at(numberControlPoints-1) ) );
    else
        endTangent = (affineHomogenous(controlPoints.back()) - affineHomogenous(controlPoints.at(numberControlPoints-2) ) ) *
                (degree*controlPoints.at(numberControlPoints-2)(dimension) / ( (knotVector.at(numberControlPoints)-knotVector.at(numberControlPoints-1) )
                                                                               * controlPoints.back()(dimension) ) );
}

void Nurbs::calculateTotalLength()
{
    length = calculateLengthInterval(beginParameter, endParameter);
}

void Nurbs::calculateBeginCurvature()
{
    beginCurvature = calculateCurvature(beginParameter);
}

void Nurbs::calculateEndCurvature()
{
    endCurvature = calculateCurvature(endParameter);
}

unsigned int Nurbs::getDegree() const
{
    return degree;
}

std::vector<double> Nurbs::getKnotVector() const
{
    return knotVector;
}

bool Nurbs::getIsHomogeneous() const
{
    return isHomogeneous;
}

unsigned int Nurbs::findSpan(double u) const
{
    if (u == knotVector.at(numberControlPoints) )
        return numberControlPoints - 1;
    unsigned int low = degree;
    unsigned int high = numberControlPoints;
    unsigned int mid = 0.5 * (low + high);
    while (u<knotVector.at(mid) || u>=knotVector.at(mid+1) )
    {
        if (u < knotVector.at(mid) )
            high = mid;
        else
            low = mid;
        mid = 0.5 * (low + high);
    }
    return mid;
}

MatrixXd Nurbs::alphaMatrix(double u, unsigned int uIndex) const
{
    MatrixXd mat(degree, degree);
    unsigned int i, j;
    for (i=0; i<degree; i++)
    {
        for (j=0; j<=i; j++)
        {
            mat(j, i) = (u-knotVector.at(uIndex-i+j) ) /
                    (knotVector.at(uIndex+1+j) - knotVector.at(uIndex-i+j) );
        }
    }
    return mat;
}

VectorXd Nurbs::tempIterative(unsigned int index) const
{
    return (controlPoints.at(index)-controlPoints.at(index-1) ) *
            (degree/(knotVector.at(index+degree)-knotVector.at(index)) );
}

VectorXd Nurbs::tempIterativeN(unsigned int index, unsigned int order) const
{
    if (0 == order)
        return controlPoints.at(index);
    else
    {
        return (tempIterativeN(index, order-1) - tempIterativeN(index-1, order-1) ) *
                (degree-order+1) / (knotVector.at(index+degree-order+1) - knotVector.at(index) );
    }
}


VectorXd Nurbs::affineHomogenous(const VectorXd &homoCoor)
{
    unsigned int dim = homoCoor.size() - 1;
    return homoCoor.head(dim) / homoCoor(dim);
}

unsigned int Nurbs::findSpan(double u, unsigned int deg, const std::vector<double> &knotVec)
{
    unsigned int numberBasisFun = knotVec.size() - deg - 2;
    if (u >= knotVec.at(numberBasisFun+1) ) // safety check.
        return numberBasisFun;
    /// Find the index which u belongs to by binary search.
    unsigned int low, high, mid;
    low = deg;
    high = numberBasisFun + 1;
    mid = 0.5 * (low + high);
    while (u < knotVec.at(mid) || u >= knotVec.at(mid+1) )
    {
        if (u < knotVec.at(mid) )
            high = mid;
        else
            low = mid;
        mid = 0.5 * (low + high);
    }
    return mid;
}

unsigned int Nurbs::findMultiplicity(double u, const std::vector<double> &knotVec)
{
    unsigned int m = 0;
    for (unsigned int i=0; i<knotVec.size(); i++)
    {
        /// Account for the numerical issue.
        if (fabs(knotVec.at(i)-u) < EPS_NUM )
            m++;
    }
    return m;
}

/// Let index is the index of u. The p+1 non-zero basis function will be N_{index-p,p}..N_{index,p}
VectorXd Nurbs::basisFunction(double u, unsigned int index, unsigned int deg,
                              const std::vector<double> &knotVec)
{
    VectorXd left = VectorXd::Zero(deg+1);
    VectorXd right = VectorXd::Zero(deg+1);
    VectorXd basisVecNonZero = VectorXd::Zero(deg+1);
    VectorXd basisVecAll = VectorXd::Zero(knotVec.size()-deg-1);
    double saved, temp;
    basisVecNonZero(0) = 1.0;
    for (unsigned int j=0; j<=deg; j++)
    {
        left(j) =u - knotVec.at(index+1-j);
        right(j) = knotVec.at(index+j) - u;
        saved = 0.0;
        for (unsigned int r=0; r<j; r++)
        {
            temp = basisVecNonZero(r) / (right(r+1) + left(j-r) );
            basisVecNonZero(r) = saved + right(r+1) * temp;
            saved = left(j-r)*temp;
        }
        basisVecNonZero(j) = saved;
    }
    for (unsigned int j=0; j<=deg; j++)
        basisVecAll(index-deg+j) = basisVecNonZero(j);
    return basisVecAll;
}

/// Basis function, by size numberParameter * numberControlPoints
MatrixXd Nurbs::basisMatrix(const VectorXd &u, unsigned int deg,
                            const std::vector<double> &knotVec)
{
    unsigned int num = u.size();
    MatrixXd mat(num, knotVec.size()-deg-1);
    for (unsigned int i=0; i<num; i++)
    {
        unsigned int index = findSpan(u(i), deg, knotVec);
        mat.row(i) = basisFunction(u(i), index, deg, knotVec);
    }
    return mat;
}
