#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  // external motion. aka noise. this is a gaussian with mean 0.
  // so we can just ignore it. no-op placeholder is put here as a reminder.
  VectorXd u = VectorXd(2);
	u << 0, 0;

  // F applies the timestep to x which contains both the 
  // current position and current velocity.
  x_ = F_ * x_ + u;

  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

/// update the state by using Kalman Filter equations
void KalmanFilter::Update(const VectorXd &z) {

  // Identity matrix.
  MatrixXd I = MatrixXd::Identity(2, 2);

  // error calculation (y). H is just so we ignore one of the dimensionsn in x. 
  VectorXd y = z - H_ * x_;
  MatrixXd Ht = H_.transpose();
  // Sensor
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  // Kalman gain.
  MatrixXd K =  P_ * Ht * Si;

  //new state
  x_ = x_ + (K * y);
  // new covariance P
  P_ = (I - K * H_) * P_;
}

// TODO: let's start with this??? 
void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
}
