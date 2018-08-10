#ifndef BSPLINEAPPROXIMATION_H
#define BSPLINEAPPROXIMATION_H
#include "Nurbs.h"

namespace CNCLite{
/*************************************************************************
 * 全局光顺的B样条类，bsplineApproximationCurve。默认采用3次B样条。
 * 根据给定数据点，使用：（1）PIA（progressive iteration approxiamtion）方法;
 * （2）插值方法；（3）WPIA方法（techsoft提供的Matrix计算特征值不准，目前未添加）；
 * （4）前面方法+控制点调节方法，来生成B样条。
 * PIA参考论文Totally positive bases and progressive iteration approximation
 * WPIA参考论文Weighted progressive iteration approximation and convergence analysis
 * **************************************************************************/
class BsplineApproximation : public Nurbs
{
public:
    BsplineApproximation();
    BsplineApproximation(const std::vector<Eigen::Vector3d> &pts, uint8_t deg = 3, double ce = 0.1,
                         double fer = 0.5, double mer = 0.5);
    ~BsplineApproximation();

public:
    // [W]PIA及直接插值求逆的方法需要返回弦长参数化得到的曲线参数。弦长参数化的参数在后续的计算中可能需要。
    Eigen::VectorXd approximatePIA(const std::vector<Eigen::Vector3d> &pts, uint8_t deg = 3,
                                       double fe = 1e-4);
    void approximateModification(const std::vector<Eigen::Vector3d> &pts, uint8_t deg = 3, double ce = 0.1,
                                 double fer = 0.5, double mer = 0.5);

    // 在节点插入中，保证每个[u(i), u(i+1))区间总有4个节点矢量参数，
    // 因此修改时，不需要另外增加控制点和节点矢量
private:
    // 注意这个函数中的ce为总弓高误差
    // 修正后的控制点和节点矢量直接存储于ctrlPtsAfterMod与knotVecAfterMod中
    void localModification(std::vector<Eigen::Vector3d> &ctrlPtsAfterMod, std::vector<double> &knotVecAfterMod,
                           double me, double mer = 0.5);
    // Calculate lambda_0 and lambda_1
    void calculateLambda(double c1, double c3, double c4, double du, double em,
                         double &lambda0, double &lambda1) const;

    // 给定数据点，采用弦长参数化方法计算节点矢量及数据点对应的曲线参数。
    void chordLengthParameterization(const std::vector<Eigen::Vector3d> &pts, uint8_t deg,
                                     std::vector<double> &knotVec, Eigen::VectorXd &paraVec);
    // 找出点集中最大的2-范数
    double maximumNormPoints(const std::vector<Eigen::VectorXd> &pts) const;

private:
    double fittingError; // 拟合误差。用于指定PIA方法的拟合误差。
    double chordError; // 弓高误差。用于指定调整控制点中的弓高误差。
    double chordErrorRatio; // 控制点调整中误差的分配系数。
};

} // End of CNCLite

#endif // BSPLINEAPPROXIMATION_H
