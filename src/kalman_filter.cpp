#include "kalman_filter.h"
#include <iostream>

#define _USE_MATH_DEFINES
#include <cmath>

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Predict() {
  // F applies the timestep to x which contains both the 
  // current position and current velocity.
  // "u" is the external motion/noise. This is gaussian with mean 0.
  // So adding it is a no-op. Thus I only show it in comments.
  x_ = F_ * x_; /* + u; */

  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

/// update the state by using Kalman Filter equations
void KalmanFilter::Update(const VectorXd &z) {
  // error calculation (y). H is just so we ignore one of the dimensions in x. 
  VectorXd y = z - H_ * x_;

  UpdateCommon(y);
}

float normalize_radians(float radians) {
  if (radians > M_PI) {
    radians = radians - M_PI;
  } else if (radians < -M_PI) {
    radians = radians + M_PI;
  } else {
    return radians;
  }
  // Recurse because we could be more than one pi off.
  return normalize_radians(radians);
}

// update the state by using Extended Kalman Filter equations
void KalmanFilter::UpdateEKF(const VectorXd &z) {
  float radius = sqrt(x_(0)*x_(0) + x_(1)*x_(1));
  float angle = atan2(x_(1), x_(0));

  float radius_derivative;

  if (fabs(radius) < 0.0001) { // value treated as 0 inside tools.h
    std::cout << " less than zero " << std::endl;
    radius_derivative = 0;
  } else {
    radius_derivative = (x_(0)*x_(2) + x_(1)*x_(3))/radius;
  }

  VectorXd z_pred(3);
  z_pred << radius, angle, radius_derivative;
  VectorXd y = z - z_pred;  
  y(1) = normalize_radians(y(1));

  UpdateCommon(y);
}

// Using the prediction_error, perform an update.
void KalmanFilter::UpdateCommon(const Eigen::VectorXd &prediction_error) {
   // Identity matrix.
  MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
   
  MatrixXd Ht = H_.transpose();
  // Sensor
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  // Kalman gain.
  MatrixXd K =  P_ * Ht * Si;

  //new state
  x_ = x_ + (K * prediction_error);
  // new covariance P
  P_ = (I - K * H_) * P_;
}