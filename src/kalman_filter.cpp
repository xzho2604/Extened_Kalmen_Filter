#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

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
  /**
   * TODO: predict the state
   */
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * TODO: update the state by using Kalman Filter equations
   */
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  int x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   */

  //given x state compute the h(x) to polar coordinate 
  // recover state parameters
  double px = x_(0);
  double py = x_(1);
  double vx = x_(2);
  double vy = x_(3);
  //precomute some terms
  double c1 = px*px+py*py;
  double c2 = sqrt(c1);
  double c3 = (c1*c2);
  
  // check division by zero
  if (fabs(c1) < 0.0001) {
      std::cout << "CalculateJacobian () - Error - Division by Zero" << std::endl;
    return;
  }

  //std::cout << "Compute Innoavation"<<std::endl;
  //compute the h(x)
  VectorXd hx(3);
  hx(0)  = c2;
  hx(1) = atan2(py,px);
  hx(2) = (px*vx + py*vy)/c2;
  //std::cout << "Compute Innoavation Passed"<<std::endl;

  //compute the innnovation
  VectorXd y = z - hx;

  //normalised the angle to -pi to pi
  while(y[1] < -M_PI)
    y[1] += 2 * M_PI;
  while (y[1] > M_PI)
    y[1] -= 2 * M_PI;

  //Mesaurement updates
  MatrixXd Ht = H_.transpose();
  MatrixXd PHt = P_ * Ht;
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K = PHt * Si;

  //New estimate
  x_ = x_ + (K * y);
  int x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;

}
