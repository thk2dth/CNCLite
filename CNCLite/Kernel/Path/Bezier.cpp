#include "Bezier.h"

using namespace CNCLite;
using namespace  Eigen;

Bezier::Bezier()
{
    curveType = BEZIER;
    isInitialized = false;
}

Bezier::Bezier(const std::vector<VectorXd> &ctrlPts, double bgnPara, double endPara)
{
    initializeFromPoints(ctrlPts, bgnPara, endPara);
}

Bezier::~Bezier()
{

}

void Bezier::initializeFromPoints(const std::vector<VectorXd> &ctrlPts, double bgnPara, double endPara)
{
    controlPoints = ctrlPts;
    beginParameter = bgnPara;
    endParameter = endPara;
    dimension = ctrlPts.front().size();
    numberControlPoints = ctrlPts.size();
    degree = numberControlPoints - 1;
    curveType = BEZIER;
    isInitialized = true;
    initialize();
}

VectorXd Bezier::calculatePoint(double u) const
{
    VectorXd temp = VectorXd::Zero(dimension);
    for (unsigned int i=0; i<=degree; i++)
    {
        temp += controlPoints.at(i) * bernsteinPolynomial(degree, i, u);
    }
    return temp;
}

VectorXd Bezier::calculateDer1(double u) const
{
    VectorXd temp = VectorXd::Zero(dimension);
    for (unsigned int i=0; i<=degree-1; i++)
    {
        temp += ( controlPoints.at(i+1) - controlPoints.at(i)) * bernsteinPolynomial(degree-1, i, u);
    }
    return temp*degree;
}

VectorXd Bezier::calculateDer2(double u) const
{
    VectorXd temp = VectorXd::Zero(dimension);
    for (unsigned int i=0; i<=degree-2; i++)
    {
        temp += ( controlPoints.at(i+2) + controlPoints.at(i) - controlPoints.at(i+1) * 2.0 ) *
                bernsteinPolynomial(degree-2, i, u);
    }
    return temp*degree;
}

void Bezier::calculateDerN(double u, unsigned int order, std::vector<VectorXd> &rst) const
{
    rst[0] = calculatePoint(u);
    if (order >= 1)
        rst[1] = calculateDer1(u);
    if (order >= 2)
        rst[2] = calculateDer2(u);
    if (order > 2)
    {
        std::cerr << "Undefined derivatives yet.";
    }
}

double Bezier::calculateLengthInterval(double us, double ue) const
{
    double len = 0.0;
    unsigned int count = 0;
    calculateLengthIntervalSimpson(us, ue, &len, &count);
    return len;
}

void Bezier::scaleByParemeter(double r)
{
    endParameter = beginParameter + (endParameter - beginParameter) * r;
    beginTangent /= r;
    endTangent /= r;
}

double Bezier::calculateNextParameter(double ds, double us, const INTERPOLATION_METHOD method) const
{
    double ue = beginParameter;
    switch (method) {
    case FIRSTTE:
        ue = calculateNextParameterNewton(ds, us);
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

double Bezier::calculateCurvature(double u) const
{
    return calculateCurvatureByDerivatives(u);
}

void Bezier::calculateBeginPoint()
{
    beginPoint = controlPoints.front();
}

void Bezier::calculateEndPoint()
{
    endPoint = controlPoints.back();
}

void Bezier::calculateBeginTangent()
{
    beginTangent = ( controlPoints.at(1) - controlPoints.at(0) ) * degree;
}

void Bezier::calculateEndTangent()
{
    endTangent = ( controlPoints.at(numberControlPoints-1) -
                   controlPoints.at(numberControlPoints-2) ) * degree;
}

void Bezier::calculateBeginCurvature()
{
    beginCurvature = calculateCurvature(beginParameter);
}

void Bezier::calculateEndCurvature()
{
    endCurvature = calculateCurvature(endParameter);
}

void Bezier::calculateTotalLength()
{
    length = calculateLengthInterval(beginParameter, endParameter);
}

/// Additional functions.
unsigned int Bezier::getDegree() const
{
    return degree;
}

double Bezier::bernsteinPolynomial(unsigned int degree, unsigned int index, double u) const
{
    return pow(u, index) * pow( (1-u), degree-index ) * factorialInterval(degree, index+1)
            / factorialInterval(degree-index, 1);
}

unsigned int Bezier::factorialInterval(unsigned int end, unsigned int start) const
{
    if (start > end)
        return 1;
    else if (end == 0)
        return 1;
    else
        return end*factorialInterval(end-1, start);
}
