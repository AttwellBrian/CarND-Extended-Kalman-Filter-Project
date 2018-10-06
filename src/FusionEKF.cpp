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
  Hj_ = MatrixXd(3, 4); // matrix with jacobian for ignoring radar

  // measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  // measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  // Project 4d state to 2d observation space.
  // Just use this to ignore the velocity field.s
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  // I straight up don't understand this ...
  Hj_ << 1, 1, 0, 0,
         1, 1, 0, 0,
         1, 1, 1, 1; 

  // The transition matrix that takes pose and converts it to a predicted pose.
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, 1, 0, // position_x_new = position_x + time*velocity_x to predict
             0, 1, 0, 1, // position_y_new = position_y + time*velocity_y to predict
             0, 0, 1, 0, // velocity_x = velocity_x
             0, 0, 0, 1; // velocity_y = velocity_y

  // Initial state covariance matrix. Ie, overall state covariance matrix.
  // Before we measure anything we have huge uncertainy about position and
  // velocity.
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ << 0xFFFFFF, 0, 0, 0,  
             0, 0xFFFFFF, 0, 0,
             0, 0, 0xFFFFFF, 0,
             0, 0, 0, 0xFFFFFF;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {

  std::cout << "measurement = " << measurement_pack.raw_measurements_ << std::endl;
  const auto& raw_measurement = measurement_pack.raw_measurements_;

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // Convert radar from polar to cartesian coordinates and initialize state.
      // I could have just ignored this for simplicity.
      float radius     = raw_measurement(0); // aka rho
      float angle    = raw_measurement(1); // aka phi
      float radial_velocity = raw_measurement(2); // aka rho_dot
      ekf_.x_(0) = radius     * cos(angle);
      ekf_.x_(1) = radius     * sin(angle);      
      ekf_.x_(2) = radial_velocity * cos(angle);
      ekf_.x_(3) = radial_velocity * sin(angle);
      
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      ekf_.x_ << 0, 0, 0, 0;  
      ekf_.x_(0) = raw_measurement.x();
      ekf_.x_(1) = raw_measurement.y();
    }

    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  // elapsed time from last measurement, raised to several powers for convenience
  uint32_t delta_time_us = (measurement_pack.timestamp_ - previous_timestamp_);
  float delta_time = delta_time_us / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;
  float delta_time_2 = delta_time   * delta_time;
  float delta_time_3 = delta_time_2 * delta_time;
  float delta_time_4 = delta_time_3 * delta_time;

  // Update the f matrix so that our predictions are for the current time.
  ekf_.F_(0, 2) = delta_time;
  ekf_.F_(1, 3) = delta_time;

  // set the acceleration noise components. 
  // known as sigma_ax in our notes.
  // It wasn't clear to me if the spec was saying noise=9 or noise squared=9.
  // So I tried both and found that the difference didn't really matter.
  // MSE = [0.08, 0.09, 0.44, 0.46] with 81
  // MSE = [0.09, 0.085, 0.42, 0.46] with 9
  float noise_ax = 9;
  float noise_ay = 9;

  // Set the process covariance matrix Q
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ << delta_time_4 / 4 * noise_ax, 0, delta_time_3 / 2 * noise_ax, 0,
             0, delta_time_4 / 4 * noise_ay, 0, delta_time_3 / 2 * noise_ay,
             delta_time_3 / 2 * noise_ax, 0, delta_time_2*noise_ax, 0,
             0, delta_time_3 / 2 * noise_ay, 0, delta_time_2*noise_ay;

  ekf_.Predict(); 

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  // this makes things much worse...
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    Tools tools;
    Hj_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.H_ = Hj_;
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(raw_measurement);
  } else {
    // Laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(raw_measurement);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
