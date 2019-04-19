#include <cfloat>
#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  // Verify that a measurement was taken.
  assert(estimations.size() > 0);
  // Verify that there are as many measurements as ground truth values.
  assert(estimations.size() == ground_truth.size());

  // Compute cumulative RMSE.
  VectorXd rmse = VectorXd::Zero(4);

  VectorXd sum = VectorXd::Zero(4);
  for (int i = 0; i < estimations.size(); ++i) {
    VectorXd error = estimations[i] - ground_truth[i];
    VectorXd square = error.array() * error.array();
    sum += square;
  }
  VectorXd mean = sum / estimations.size();

  rmse << sqrt(mean.array());
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
  // Converting from state space into radar measurement space.
  MatrixXd Hj(3, 4);
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  // Check for division by 0.
  if (fabs(px) < FLT_EPSILON && fabs(py) < FLT_EPSILON) {
    Hj = MatrixXd::Zero(3, 4);
    return Hj;
  }

  float norm_2 = px * px + py * py;
  float norm = sqrt(norm_2);
  float norm_3_2 = norm_2 * norm;

  Hj(0, 0) = px / norm;
  Hj(0, 1) = py / norm;
  Hj(0, 2) = 0;
  Hj(0, 3) = 0;

  Hj(1, 0) = -py / norm_2;
  Hj(1, 1) = px / norm_2;
  Hj(1, 2) = 0;
  Hj(1, 3) = 0;

  Hj(2, 0) = py * (vx * py - vy * px) / norm_3_2;
  Hj(2, 1) = px * (vy * px - vx * py) / norm_3_2;
  Hj(2, 2) = Hj(0, 0);
  Hj(2, 3) = Hj(0, 1);

  return Hj;
}
