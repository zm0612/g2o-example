#include <Eigen/Core>
#include <iostream>

#include "g2o/stuff/sampler.h"
#include "g2o/core/sparse_optimizer.h"
#include "g2o/core/block_solver.h"
#include "g2o/core/optimization_algorithm_levenberg.h"
#include "g2o/core/base_vertex.h"
#include "g2o/core/base_unary_edge.h"
#include "g2o/solvers/dense/linear_solver_dense.h"
#include "g2o/core/robust_kernel_impl.h"

using namespace std;

/*!
 * 继承BaseVertex类，构造顶点
 */
class VertexParams : public g2o::BaseVertex<3, Eigen::Vector3d> {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

    VertexParams() = default;

    bool read(std::istream & /*is*/) override {
        cerr << __PRETTY_FUNCTION__ << " not implemented yet" << endl;
        return false;
    }

    bool write(std::ostream & /*os*/) const override {
        cerr << __PRETTY_FUNCTION__ << " not implemented yet" << endl;
        return false;
    }

    //该函数作用是更新顶点的估计值
    void setToOriginImpl() override {
        cerr << __PRETTY_FUNCTION__ << " not implemented yet" << endl;
    }

    //更新优化之后的顶点
    void oplusImpl(const double *update) override {
        Eigen::Vector3d::ConstMapType v(update);
        _estimate += v;
    }
};

/*!
 * 从BaseUnaryEdge继承得到一元边
 */
class EdgePointOnCurve : public g2o::BaseUnaryEdge<1, Eigen::Vector2d, VertexParams> {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    EdgePointOnCurve() = default;

    bool read(std::istream & /*is*/) override {
        cerr << __PRETTY_FUNCTION__ << " not implemented yet" << endl;
        return false;
    }

    bool write(std::ostream & /*os*/) const override {
        cerr << __PRETTY_FUNCTION__ << " not implemented yet" << endl;
        return false;
    }

    //边的误差计算
    void computeError() override {
        const VertexParams *params = dynamic_cast<const VertexParams *>(vertex(0));//顶点
        const double &a = params->estimate()(0);
        const double &b = params->estimate()(1);
        const double &lambda = params->estimate()(2);
        double fval = a * exp(-lambda * measurement()(0)) + b;
        _error(0) = std::abs(fval - measurement()(1));
    }
};

int main(int argc, char **argv) {
    int numPoints = 50;
    int maxIterations = 50;
    bool verbose = true;

    double a = 2.;
    double b = 0.4;
    double lambda = 0.2;
    Eigen::Vector2d *points = new Eigen::Vector2d[numPoints];
    ofstream points_file("../points.txt", ios::out);

    //准备用于拟合的数据
    for (int i = 0; i < numPoints; ++i) {
        double x = g2o::Sampler::uniformRand(0, 10);
        double y = a * exp(-lambda * x) + b;
        y += g2o::Sampler::gaussRand(0, 0.02);

        if (i == 20) {
            x = 8;
            y = 2.5;
        }

        points[i].x() = x;
        points[i].y() = y;
        points_file << x << " " << y << endl;
    }
    points_file.close();

    typedef g2o::BlockSolver<g2o::BlockSolverTraits<Eigen::Dynamic, Eigen::Dynamic> > MyBlockSolver;
    typedef g2o::LinearSolverDense<MyBlockSolver::PoseMatrixType> MyLinearSolver;

    g2o::SparseOptimizer optimizer;

    g2o::OptimizationAlgorithmLevenberg *solver = new g2o::OptimizationAlgorithmLevenberg(
            g2o::make_unique<MyBlockSolver>(g2o::make_unique<MyLinearSolver>()));

    optimizer.setAlgorithm(solver);

    VertexParams *params = new VertexParams();
    params->setId(0);
    params->setEstimate(Eigen::Vector3d(1, 1, 1));//初始化顶点的估计值
    optimizer.addVertex(params);

    for (int i = 0; i < numPoints; ++i) {
        EdgePointOnCurve *e = new EdgePointOnCurve;
        e->setInformation(Eigen::Matrix<double, 1, 1>::Identity());

        if (i == 20) {
            e->setInformation(Eigen::Matrix<double, 1, 1>::Identity() * 10);
        }

        e->setVertex(0, params);
        e->setMeasurement(points[i]);

        g2o::RobustKernelHuber *robust_kernel_huber = new g2o::RobustKernelHuber;
        robust_kernel_huber->setDelta(0.3);
        e->setRobustKernel(robust_kernel_huber);

        optimizer.addEdge(e);
    }

    optimizer.initializeOptimization();
    optimizer.setVerbose(verbose);
    optimizer.optimize(maxIterations);

    ofstream result_file("../result.txt");
    result_file << params->estimate()[0] << " "
                << params->estimate()[1] << " "
                << params->estimate()[2];
    result_file.close();


    cout << endl << "a, b, lambda: "
         << params->estimate()[0] << ", "
         << params->estimate()[1] << ", "
         << params->estimate()[2] << endl;

    delete[] points;

    return 0;
}
