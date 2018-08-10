#include "BsplineApproximation.h"
using namespace CNCLite;
using namespace Eigen;

BsplineApproximation::BsplineApproximation()
{
    fittingError = 1e-4;
    chordError = 0.1;
    chordErrorRatio = 0.5;
    degree = 3;
}

BsplineApproximation::BsplineApproximation(const std::vector<Vector3d> &pts, uint8_t deg,
                                           double ce, double fer, double mer)
{
    approximateModification(pts, deg, ce, fer, mer);
}

BsplineApproximation::~BsplineApproximation()
{

}

VectorXd BsplineApproximation::approximatePIA(const std::vector<Eigen::Vector3d> &pts,
                                                         uint8_t deg, double fe)
{
    degree = deg;
    fittingError = fe;
    numberControlPoints = pts.size();
    uint32_t numberKnotVector = numberControlPoints + degree + 1;
    VectorXd paraVec = VectorXd::Zero(numberControlPoints); // 对应控制点的曲线参数
    knotVector.assign(numberKnotVector, 0.0);
    chordLengthParameterization(pts, deg, knotVector, paraVec);
    MatrixXd basisMat = basisMatrix(paraVec, deg, knotVector);
    // 迭代过程的误差量，初始化为零
    std::vector<VectorXd> deltaPts(numberControlPoints, Vector3d::Zero() );
    std::vector<VectorXd> ctrlPts(numberControlPoints, Vector3d::Zero() ); // 迭代过程的控制点
    std::vector<VectorXd> ptsTemp(numberControlPoints, Vector3d::Zero() ); // 迭代过程的中间点
    for (uint32_t i=0; i<numberControlPoints; i++)
    {
        ctrlPts[i] = pts.at(i);
        ptsTemp[i] = pts.at(i);
    }
    double maxDis = std::numeric_limits<double>::max();
    uint32_t itrNum = 0;
    do
    {
        itrNum++;
        for (uint32_t i = 0; i < numberControlPoints; i++)
        {
            ctrlPts[i] += deltaPts[i];
            ptsTemp[i] = ctrlPts[0] * basisMat(i, 0); // 如此是为了方便初始化ptsTemp[i]
            for (uint32_t j = 1; j < numberControlPoints; j++)
            {
                ptsTemp[i] += ctrlPts[j] * basisMat(i, j);
            }
            deltaPts[i] = pts[i] - ptsTemp[i];
        }
        maxDis = maximumNormPoints(deltaPts);

    }while(maxDis >= fe);
    std::cout << "Iteration number for PIA is: " << itrNum << "." << std::endl;
    initializeFromPoints(ctrlPts, knotVector, false);
    return paraVec;
}

void BsplineApproximation::approximateModification(const std::vector<Vector3d> &pts,
                                                   uint8_t deg, double ce, double fer, double mer)
{
    chordError = ce;
    chordErrorRatio = 1 - fer;
    uint32_t numPara = pts.size();
    VectorXd paraVec = approximatePIA(pts, deg, fer*ce);
    // 插入弦长参数化时得到的曲线参数。曲线参数一个个插入。
    // 对于每一个参数，先计算在已有节点矢量中的重复度r，然后插入degree-r次。
    uint32_t i, j;
    std::vector<uint8_t> mulVec(numPara, 0); // 曲线参数的重复度
    // 中间四点到首末两点构成直线的距离标志
    //   flag, determine the points to be altered.
    //           1 (0001), the first point; 2 (0010), the second one;
    //          4 (0100), the third one; 8 (1000), the fourth one.
    //           0 (0000), none.
    //           e.g., flag = 6 implies the second and third points are out of the constraints.
    std::vector<uint8_t> flagVec(numPara-1, 0);
    // 节点插入，同时保存所有的节点重复度。
    for (i=0; i<numPara-1; i++)
    {
        mulVec[i] = findMultiplicity(paraVec[i], knotVector);
        knotVectorInsert(paraVec[i], degree-mulVec[i]);
    }
    uint32_t startIndex, endIndex, numIndex; // 曲线参数u_i, u_(i+1)对应的节点区间。
    endIndex = knotVector.size();
    startIndex = 0;
    numIndex = endIndex - startIndex;
    // 补全其他节点矢量，保证曲线参数u_i, u_(i+1)对应的节点区间共有8个控制点
    double uMid; // 中间节点矢量，用于判断节点的插入位置和插入的节点值
    for (i=0; i<numPara-1; i++)
    {
        for (j=startIndex; j<knotVector.size(); j++)
            if (knotVector.at(j) == paraVec(i+1) )
            {
                endIndex = j;
                break;
            }
        for (j=endIndex; j>0; j--)
            if (knotVector.at(j) == paraVec(i) )
            {
                startIndex = j;
                break;
            }
        uMid = 0.5 * (paraVec(i) + paraVec(i+1) );
        numIndex = endIndex - startIndex;
        if (1 == numIndex)
        {
            // 需插入四个节点，((5-j)*u_i + j*u_(i+1))/5, j=1,2,3,4
            for(j=1; j<5; j++)
                knotVectorInsert(((5-j)*paraVec(i) + j*paraVec(i+1) ) * 0.2);
            endIndex += 4; // 插入4个节点后，最后的节点索引也要+4
        }
        else if(2 == numIndex)
        {
            // 插入三个节点。分两种情况。
            if(knotVector[startIndex+1] < uMid)
            {
                for(j=1; j<4; j++)
                {
                    knotVectorInsert(((4-j)*knotVector[startIndex+1] + j*paraVec(i+1) ) * 0.25);
                }
            }
            else
            {
                for(j=1; j<4; j++)
                {
                    knotVectorInsert(((4-j)*paraVec(i) + j*knotVector[startIndex+1]) * 0.25);
                }
            }
            endIndex += 3; // 插入3个节点后，最后的节点索引也要+3
        }
        else if(3 == numIndex)
        {
            // 插入两个节点。分四种情况。
            if(knotVector[startIndex+2] <= uMid)
            {
                knotVectorInsert((2*knotVector[startIndex+2] + paraVec(i+1) ) / 3.0);
                knotVectorInsert((knotVector[startIndex+2] + 2*paraVec(i+1) ) / 3.0);
            }
            else if(knotVector[startIndex+1] >= uMid)
            {
                knotVectorInsert((2*paraVec(i) + knotVector[startIndex+1]) / 3.0);
                knotVectorInsert((paraVec(i) + 2*knotVector[startIndex+1]) / 3.0);
            }
            else
            {
                // 分两种情况
                if(knotVector[startIndex+2]-knotVector[startIndex+1] >=
                        (paraVec(i+1)-paraVec(i))/3.0)
                {
                    knotVectorInsert((2*knotVector[startIndex+1] + knotVector[startIndex+2]) / 3.0);
                    knotVectorInsert((knotVector[startIndex+1] + 2*knotVector[startIndex+2]) / 3.0);
                }
                else
                {
                    knotVectorInsert((paraVec(i) + knotVector[startIndex+1]) * 0.5);
                    knotVectorInsert((knotVector[startIndex+2] + paraVec(i+1) ) * 0.5);
                }
            }
            endIndex += 2; // 插入2个节点后，最后的索引也要+2
        }
        else
            printf("In knot insertion, the number of knot parameters to "
                   "be inserted is greater than 2 between %f and %f.\n",
                   knotVector[startIndex], knotVector[endIndex] );
    }
    // 判断8个控制点中的中间6个控制点到线段的弓高误差
    startIndex = 0;
    double dis;
    // 生成controlPoint和knotVector的副本
    std::vector<Vector3d> ctrlPtsTemp(8, Vector3d::Zero() );
    std::vector<double> knotVecTemp(6, 0.0);
    Vector3d lineVec = Vector3d::Zero();
    for(i=0; i<numPara-1; i++)
    {
        for(j=0; j<8; j++)
            // 定义曲线的8个控制控制点
            ctrlPtsTemp[j] = controlPoints[7*i+j];
        for(j=0; j<6; j++)
            // 定义曲线的6个控制点。u(i)和u(i+1)均只记一次。
            knotVecTemp[j] = knotVector[7*i+3+j];
        // 注意判断的是控制点到**原始线段**的距离
        lineVec = (pts[i+1]-pts[i]).normalized();
        for(j=1; j<7; j++)
        {
            // 注意判断的是控制点到原始线段的距离
            // The line used in distanceToLine function is defined
            // by a point and a unit vector direction.
            /// dis = ctrlPtsTemp[j].distanceToLine(pts[i], lineVec);
            Vector3d vecDis = ctrlPtsTemp[j] - pts[i];
            dis = vecDis.cross(lineVec).norm();
            if(dis > ce)
            {
                flagVec[i] += 1<<(j-1); // 左移位代替幂乘
            }
        }
        // 有控制点不满足误差要求，则调整该区间的控制点
        if(0 != flagVec[i])
        {
            localModification(ctrlPtsTemp, knotVecTemp, ce*(1-fer), mer); /// 保证误差总量为弓高误差
            // 更新局部调整前的控制点和节点矢量
            for(j=1; j<7; j++)
            {
                controlPoints[7*i+j] = ctrlPtsTemp[j];
            }
            for(j=1; j<5; j++)
                knotVector[7*i+3+j] =  knotVecTemp[j];
        }
    }
    // 调用初始化函数
    initializeFromPoints(controlPoints, knotVector, false);
}

// 调整控制点。
// 注意这个函数中的ce为总弓高误差
void BsplineApproximation::localModification(std::vector<Vector3d> &ctrlPtsAfterMod,
                                             std::vector<double> &knotVecAfterMod, double me,
                                             double mer)
{
    std::cout << "The curve is locally modified based on the strong convex hull property." << std::endl;
    if (ctrlPtsAfterMod.size() !=8 || knotVecAfterMod.size() != 6)
    {
        std::cout << "The number of control points and of knot vector is not 8 or 6." << std::endl;
        return;
    }
    double d0, d1, d2, d3;
    Vector3d startDer1, startDer2, endDer1, endDer2;
    double u1 = knotVecAfterMod[0];
    double v1 = knotVecAfterMod[1];
    double v2 = knotVecAfterMod[2];
    double v3 = knotVecAfterMod[3];
    double v4 = knotVecAfterMod[4];
    double u2 = knotVecAfterMod[5];
    Vector3d R0 = ctrlPtsAfterMod[0];
    Vector3d R1 = ctrlPtsAfterMod[1];
    Vector3d R2 = ctrlPtsAfterMod[2];
    // Vector3d R3 = ctrlPtsAfterMod[3];
    // Vector3d R4 = ctrlPtsAfterMod[4];
    Vector3d R5 = ctrlPtsAfterMod[5];
    Vector3d R6 = ctrlPtsAfterMod[6];
    Vector3d R7 = ctrlPtsAfterMod[7];
    Vector3d lineVec = R7 - R0;

    // double cf = me*mer; // 为保证调整后的节点参数尽可能均分分布，第一个点处的误差需乘一个系数
    startDer1 = (R1-R0) * (3.0/(v1-u1));
    startDer2 = (R2-R1)*(6.0/((v1-u1)*(v2-u1))) - (R1-R0)*(6.0/((v1-u1)*(v1-u1)));
    endDer1 = (R7-R6) * (3.0/(u2-v4));
    endDer2 = (R7-R6)*(6.0/((u2-v4)*(u2-v4))) - (R6-R5)*(6.0/((u2-v4)*(u2-v3)));
    double du = 0.2*(u2-u1);
    Vector3d n11, n12, n21, n22;
    double c1, c2, c3, c4, c5, c6, c7, c8;
    // 可以先计算n12
    n12 = startDer1.cross(lineVec).normalized();
    n11 = n12.cross(lineVec).normalized();
    // 可以先计算n22
    n22 = endDer1.cross(lineVec).normalized();
    n21 = n22.cross(lineVec).normalized();
    c1 = n11.dot(startDer1);
    c2 = n12.dot(startDer1); // c2 = 0
    c3 = n11.dot(startDer2);
    c4 = n12.dot(startDer2); // 对于平面情况，c4 = 0
    c5 = n21.dot(endDer1);
    c6 = n22.dot(endDer1); // c6 = 0
    c7 = n21.dot(endDer2);
    c8 = n22.dot(endDer2);// 对于平面情况，c8 = 0
    calculateLambda(c1, c3, c4, du, me, d0, d1);
    calculateLambda(-c5, c7, c8, du, me, d3, d2);
    // 更新之后的节点矢量和控制点
    // 插入的节点矢量记为w0..3
    knotVecAfterMod[1] = u1 + d0; // w0
    knotVecAfterMod[4] = u2 - d3; // w3
    knotVecAfterMod[2] = u1 + d1; // w1
    knotVecAfterMod[3] = u2 - d2; // w2
    // 插入的控制点记为S1..6
    ctrlPtsAfterMod[1] = R0 + startDer1 * (d0/3.0); // S1
    ctrlPtsAfterMod[6] = R7 - endDer1 * (d3/3.0); // S6
    ctrlPtsAfterMod[2] = R0 + startDer1*((d0+d1)/3.0) + startDer2*((d0*d1)/6.0); // S2
    ctrlPtsAfterMod[5] = R7 - endDer1*((d2+d3)/3.0) + endDer2*(d2*d3/6.0); // S5
    ctrlPtsAfterMod[3] = ctrlPtsAfterMod[2]*(2.0/3) + ctrlPtsAfterMod[5]*(1.0/3);
    ctrlPtsAfterMod[4] = ctrlPtsAfterMod[2]*(1.0/3) + ctrlPtsAfterMod[5]*(2.0/3);
}

void BsplineApproximation::chordLengthParameterization(const std::vector<Vector3d> &pts,
                                                       uint8_t deg, std::vector<double> &knotVec,
                                                       VectorXd &paraVec)
{
    uint32_t num = pts.size();
    paraVec = VectorXd::Zero(num);
    knotVec.resize(num+deg+1);
    paraVec[0] = 0.0;
    uint32_t i, j;
    for (i=1; i<num; i++) paraVec[i] = paraVec[i-1] + (pts[i]-pts[i-1]).norm();
    for (i=1; i<num; i++) paraVec[i] = paraVec[i] / paraVec[num-1];
    // 节点矢量
    for (i=0; i<=deg; i++) knotVec[i] = paraVec[0];
    for (i=num; i<=num+deg; i++) knotVec[i] = paraVec[num-1];
    for (j=1; j<=num-deg-1; j++)
    {
        knotVec[j+deg] = 0.0;
        for (i=j; i<=j+deg-1; i++)
            knotVec[j+deg] += paraVec[i];
        knotVec[j+deg] /= deg;
    }
}

double BsplineApproximation::maximumNormPoints(const std::vector<VectorXd> &pts) const
{
    double dis, temp;
    dis = 0.0;
    uint32_t i = 0;
    while (i<pts.size() )
    {
        temp = pts.at(i).norm();
        if (temp > dis)
            dis = temp;
        i++;
    }
    return dis;
}

void BsplineApproximation::calculateLambda(double c1, double c3, double c4, double du,
                                           double em, double &lambda0, double &lambda1) const
{
    double x1, x2;
    if (c3 == 0)
        x1 = 1.5 * em / std::fabs(c1);
    else if (c3 < 0)
        x1 = du;
    else // c3 > 0
    {
        if (c1 > 0)
            x1 = (-2*c1 + sqrt(4*c1*c1+6*c3*em) ) / c3;
        else
        {
            if (c1*c1 > 1.5*c3*em)
                x1 = (-2*c1 - sqrt(4*c1*c1-6*c3*em) ) / c3;
            else
                x1 = (-2*c1 + sqrt(4*c1*c1+6*c3*em) ) / c3;
        }
    }
    lambda0 = std::fmin(std::fmin(du, 1.5*em/fabs(c1) ), x1);
    double k1 = c1*lambda0/3.0;
    double k2 = c1/3.0 + c3*lambda0/6.0;
    double k3 = c4*lambda0/6.0;
    x2 = (-k1*k2+sqrt(k2*k2*em*em+k3*k3*(em*em-k1*k1))) / (k2*k2+k3*k3);
    lambda1 = std::fmin(2*du, x2);
}
