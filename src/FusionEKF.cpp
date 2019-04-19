#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  // Hj_ will be reset after each prediction step in ProcessMeasurement.
  Hj_ = MatrixXd::Zero(3, 4);
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    // Based on the measurement covariance matrices, lidar seems to give quite
    // accurate position data but has no idea about velocity, whereas radar seems
    // to be pretty accurate about both. Therefore, initialize the state
    // covariance matrix to have low position errors and high(er) velocity errors.
    MatrixXd P(4, 4);
    P << 1, 0, 0, 0,
         0, 1, 0, 0,
         0, 0, 3, 0,
         0, 0, 0, 3;

    MatrixXd F = MatrixXd::Zero(4, 4);
    MatrixXd Q = MatrixXd::Zero(4, 4);

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float rho = measurement_pack.raw_measurements_(0);
      float phi = measurement_pack.raw_measurements_(1);
      float rho_dot = measurement_pack.raw_measurements_(2);

      VectorXd x(4);
      x(0) = rho * cos(phi);
      x(1) = rho * sin(phi);
      x(2) = rho_dot * cos(phi);
      x(3) = rho_dot * sin(phi);

      MatrixXd H = Hj_;
      MatrixXd R = R_radar_;

      ekf_.Init(x, P, F, H, R, Q);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      float px = measurement_pack.raw_measurements_(0);
      float py = measurement_pack.raw_measurements_(1);

      VectorXd x(4);
      x(0) = px;
      x(1) = py;
      x(2) = 0;
      x(3) = 0;

      MatrixXd H = H_laser_;
      MatrixXd R = R_laser_;

      ekf_.Init(x, P, F, H, R, Q);
    }

    // done initializing, no need to predict or update
    previous_timestamp_ = measurement_pack.timestamp_;
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  float noise_ax = 9.0;
  float noise_ay = 9.0;
  float delta_t = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;

  MatrixXd F(4, 4);
  F << 1, 0, delta_t, 0,
        0, 1, 0, delta_t,
        0, 0, 1, 0,
        0, 0, 0, 1;
  ekf_.F_ = F;

  float delta_t_2 = delta_t * delta_t;
  float delta_t_3 = delta_t_2 * delta_t;
  float delta_t_4 = delta_t_3 * delta_t;
  MatrixXd Q(4, 4);
  Q(0, 0) = delta_t_4 / 4 * noise_ax;
  Q(0, 1) = 0;
  Q(0, 2) = delta_t_3 / 2 * noise_ax;
  Q(0, 3) = 0;

  Q(1, 0) = 0;
  Q(1, 1) = delta_t_4 / 4 * noise_ay;
  Q(1, 2) = 0;
  Q(1, 3) = delta_t_3 / 2 * noise_ay;

  Q(2, 0) = delta_t_3 / 2 * noise_ax;
  Q(2, 1) = 0;
  Q(2, 2) = delta_t_2 * noise_ax;
  Q(2, 3) = 0;

  Q(3, 0) = 0;
  Q(3, 1) = delta_t_3 / 2 * noise_ay;
  Q(3, 2) = 0;
  Q(3, 3) = delta_t_2 * noise_ay;
  ekf_.Q_ = Q;

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */
  VectorXd z = measurement_pack.raw_measurements_;

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    Hj_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.H_ = Hj_;
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(z);
  } else {
    // Laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(z);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
