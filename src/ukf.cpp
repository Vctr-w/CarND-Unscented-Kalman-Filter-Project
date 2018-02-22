#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

#define EPS 0.0001

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  //std_a_ = 30;
  std_a_ = 1;

  // Process noise standard deviation yaw acceleration in rad/s^2
  //std_yawdd_ = 30;
  std_yawdd_ = 0.2;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  is_initialized_ = false;

  //set state dimension
  n_x_ = 5;

  //set augmented dimension
  n_aug_ = 7;

  //define spreading parameter
  lambda_ = 3 - n_aug_;

  previous_timestamp_ = 0;

  weights_ = VectorXd(2 * n_aug_ + 1);

  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

  if (!use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER) {
    return;
  }

  if (!use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    return;
  }

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float rho = meas_package.raw_measurements_(0);
      float phi = meas_package.raw_measurements_(1);
      float rho_dot = meas_package.raw_measurements_(2);

      x_(0) = rho * cos(phi);
      x_(1) = rho * sin(phi);
      x_(2) = 0;
      x_(3) = 0;
      x_(4) = 0;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */

      x_(0) = meas_package.raw_measurements_(0);
      x_(1) = meas_package.raw_measurements_(1);
      x_(2) = 0;
      x_(3) = 0;
      x_(4) = 0;
    }

    P_ = MatrixXd::Identity(n_x_, n_x_);
    /*
    P_ << 1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1000, 0,
        0, 0, 0, 1000;
    */

    previous_timestamp_ = meas_package.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  double delta_t = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0; //dt - expressed in seconds
  previous_timestamp_ = meas_package.timestamp_;

  if (delta_t > 0) {
    Prediction(delta_t);

    /*****************************************************************************
     *  Update
     ****************************************************************************/

    UpdateGeneric(meas_package, meas_package.sensor_type_);
    /*
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      UpdateRadar(meas_package);
    } else {
      UpdateLidar(meas_package);
    }*/
  }

  return;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  AugmentedSigmaPoints(&Xsig_aug);

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
      Xsig_pred_.col(i) = Xsig_aug.col(i).head(5);
      
      float v = Xsig_aug(2, i);
      float psi = Xsig_aug(3, i);
      float psi_dot = Xsig_aug(4, i);
      float v_a = Xsig_aug(5, i);
      float v_psi_dotdot = Xsig_aug(6, i);
      
      Xsig_pred_(3, i) += psi_dot * delta_t;

      if (Xsig_aug(4, i) != 0) {
          Xsig_pred_(0, i) += (v / psi_dot) * (sin(psi + psi_dot * delta_t) - sin(psi));
          Xsig_pred_(1, i) += (v / psi_dot) * (-cos(psi + psi_dot * delta_t) + cos(psi));
      }
      else {
          Xsig_pred_(0, i) += v * cos(psi) * delta_t;
          Xsig_pred_(1, i) += v * sin(psi) * delta_t;
      }
      
      Xsig_pred_(0, i) += 0.5 * delta_t * delta_t * cos(psi) * v_a;
      Xsig_pred_(1, i) += 0.5 * delta_t * delta_t * sin(psi) * v_a;
      Xsig_pred_(2, i) += delta_t * v_a;
      Xsig_pred_(3, i) += 0.5 * delta_t * delta_t * v_psi_dotdot;
      Xsig_pred_(4, i) += delta_t * v_psi_dotdot;
  }
    
  //set weights
  
  weights_.fill(0.5 / (lambda_ + n_aug_));
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  
  //predict state mean
  
  x_ = Xsig_pred_ * weights_;
  
  //predict state covariance matrix
  VectorXd ones = VectorXd(2 * n_aug_ + 1);
  ones.fill(1);
  
  P_ = (Xsig_pred_ - x_ * ones.transpose()) * 
    (weights_.asDiagonal()) * 
    (Xsig_pred_ - x_ * ones.transpose()).transpose();
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
}

void UKF::UpdateGeneric(MeasurementPackage meas_package, MeasurementPackage::SensorType sensor_type) {
  int n_z;

  if (sensor_type == MeasurementPackage::RADAR) {
    n_z = 3;
  }
  else if (sensor_type == MeasurementPackage::LASER) {
    n_z = 2;
  }

  VectorXd z = meas_package.raw_measurements_;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    float p_x = Xsig_pred_(0, i);
    float p_y = Xsig_pred_(1, i);
    float v = Xsig_pred_(2, i);
    float psi = Xsig_pred_(3, i);
    float psi_dot = Xsig_pred_(4, i);
    
    double v1 = cos(psi) * v;
    double v2 = sin(psi) * v;
    
    if (sensor_type == MeasurementPackage::RADAR) {
      Zsig(0, i) = sqrt(p_x * p_x + p_y * p_y);
      Zsig(1, i) = atan2(p_y, p_x);
      Zsig(2, i) = (p_x * v1 + p_y * v2 ) / sqrt(p_x * p_x + p_y * p_y);
    }
    else if (sensor_type == MeasurementPackage::LASER) {
      Zsig(0, i) = p_x;
      Zsig(1, i) = p_y;
    }
  }
  
  //calculate mean predicted measurement
  
  z_pred = Zsig * weights_;
  
  //calculate innovation covariance matrix S
  VectorXd ones = VectorXd(2 * n_aug_ + 1);
  ones.fill(1);
  
  MatrixXd R = MatrixXd(n_z, n_z);

  if (sensor_type == MeasurementPackage::RADAR) {
    R <<    std_radr_ * std_radr_, 0, 0,
          0, std_radphi_ * std_radphi_, 0,
          0, 0, std_radrd_ * std_radrd_;
  }
  else if (sensor_type == MeasurementPackage::LASER) {
    R <<    std_laspx_ * std_laspx_, 0,
          0, std_laspy_ * std_laspy_;
  }

  S = (Zsig - z_pred * ones.transpose()) * 
    (weights_.asDiagonal()) *
    (Zsig - z_pred * ones.transpose()).transpose() + R;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  Tc = (Xsig_pred_ - x_ * ones.transpose()) *
    weights_.asDiagonal() *
    (Zsig - z * ones.transpose()).transpose();
  
  //calculate Kalman gain K;
  
  MatrixXd K = Tc * S.inverse();
  
  //update state mean and covariance matrix
  
   //residual
  VectorXd z_diff = z - z_pred;
  
  //angle normalization
  while (z_diff(1) > M_PI) {
    z_diff(1) -= 2. * M_PI;
  }

  while (z_diff(1) < -M_PI) {
    z_diff(1) += 2. * M_PI;
  }

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();
}

void UKF::AugmentedSigmaPoints(MatrixXd* Xsig_out) {
  //create augmented mean vector
  VectorXd x_aug = VectorXd(7);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(7, 7);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
 
  //create augmented mean state
  
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;
  
  //create augmented covariance matrix
  
  P_aug.fill(0);
  P_aug.topLeftCorner(5, 5) = P_;
  P_aug(5,5) = std_a_ * std_a_;
  P_aug(6,6) = std_yawdd_ * std_yawdd_;

  //create square root matrix
  
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  
  Xsig_aug.col(0) = x_aug;
  
  for (int i = 0; i < n_aug_; i++) {
      Xsig_aug.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
      Xsig_aug.col(n_aug_ + i + 1) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
  }

  *Xsig_out = Xsig_aug;
}
