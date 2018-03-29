#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

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
  std_a_ = 0.3;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.3;
  
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
//  cout<< "ProcessMeasurement starts" << endl;
  if(!is_initialized_){
	  x_ <<  0,0,0,0,0;
		 
		P_ << 	   1,0,0,0,0,
					0,1,0,0,0,
					0,0,1,0,0,
					0,0,0,1,0,
					0,0,0,0,1;
		n_x_ = 5;
		n_aug_ = 7;
	
	if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
	  //cout<< "initializing radar"<< endl;
	  float ro = meas_package.raw_measurements_(0);
	  float phi = meas_package.raw_measurements_(1);
	  float ro_dot = meas_package.raw_measurements_(2);
	  x_(0) = ro*cos(phi);
	  x_(1) = ro*sin(phi);
	  //cout<< "initializing radar done"<< endl;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
	  // cout<< "initializing lidar"<< endl;
	  x_(0) = meas_package.raw_measurements_(0);
	  x_(1) = meas_package.raw_measurements_(1);
	 // cout<< "initializing lidar done"<< endl;
    }
    time_us_ = meas_package.timestamp_;
    // done initializing, no need to predict or update
    is_initialized_ = true;
	//cout << "initialization done" << endl; 
    return;

	 
  }
  //cout<<"Before prediction"<< endl;
   double delta_t = (meas_package.timestamp_ - time_us_)/1000000.0;
   time_us_ = meas_package.timestamp_;
   Prediction(delta_t);
   //cout<<"After prediction"<< endl;
   //cout<< meas_package.sensor_type_ <<endl;
   if(meas_package.sensor_type_ == MeasurementPackage::RADAR){
	   UpdateRadar(meas_package);
   }
   else{
	   UpdateLidar(meas_package);
   }
   //cout<< "ProcessMeasurement end" << endl;
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
  //cout << "Prediction starts"<< endl;
  VectorXd x_aug_ = VectorXd(7);
  
  MatrixXd P_aug_ = MatrixXd(7,7);
  
  MatrixXd Xsig_aug_ = MatrixXd(n_aug_, 2*n_aug_ + 1);
  
  x_aug_.head(5) = x_;
  x_aug_(5) = 0;
  x_aug_(6) = 0;
  P_aug_.fill(0.0);
  P_aug_.topLeftCorner(5,5) = P_;
  P_aug_(5,5) = std_a_*std_a_;
  P_aug_(6,6) = std_yawdd_*std_yawdd_;
  
  lambda_ = 3-n_aug_;
  
  MatrixXd L = P_aug_.llt().matrixL();
  Xsig_aug_.col(0) = x_aug_;
  
  for(int i=0; i<n_aug_; i++){
	  Xsig_aug_.col(i+1) = x_aug_ + sqrt(lambda_ + n_aug_)*L.col(i);
	  Xsig_aug_.col(i+n_aug_+1) = x_aug_ - sqrt(lambda_ + n_aug_)*L.col(i);
  }
  //cout<< "flag1"<< endl;
  //create matrix with predicted sigma points as columns
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  
  for(int i = 0;i< 2*n_aug_+1; i++){
      double px = Xsig_aug_(0,i);
      double py = Xsig_aug_(1,i);
      double v = Xsig_aug_(2,i);
      double yaw = Xsig_aug_(3,i);
      double yawr = Xsig_aug_(4,i);
      double nu_a = Xsig_aug_(5,i);
      double nu_yaw = Xsig_aug_(6,i);
      
      if(fabs(yawr)> 0.001){
          px += (v/yawr)*(sin(yaw +  yawr*delta_t)-sin(yaw)) + 0.5*delta_t*delta_t*cos(yaw)*nu_a;
          py += (v/yawr)*(-cos(yaw +  yawr*delta_t)+cos(yaw)) + 0.5*delta_t*delta_t*sin(yaw)*nu_a;
          v += delta_t*nu_a;
          yaw += yawr*delta_t + 0.5*delta_t*delta_t*nu_yaw;
          yawr += delta_t*nu_yaw;
      }
      else{
          px += (v)*cos(yaw)*delta_t + 0.5*delta_t*delta_t*cos(yaw)*nu_a;
          py += (v)*sin(yaw)*delta_t + 0.5*delta_t*delta_t*sin(yaw)*nu_a;
          v += delta_t*nu_a;
          yaw += 0.5*delta_t*delta_t*nu_yaw;
          yawr += delta_t*nu_yaw;
      }
      Xsig_pred_.col(i) << px,py,v,yaw,yawr;
  }
  
  weights_ = VectorXd(2*n_aug_ +1);
  
  weights_(0) = lambda_/(lambda_ + n_aug_); 
  
  for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights_
    weights_(i) = 0.5/(n_aug_+lambda_);
  }
  
  //predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    x_ = x_+ weights_(i) * Xsig_pred_.col(i);
  }

  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
  }
 // cout << "Prediction done" << endl;
  //cout << Xsig_pred_.size()<<endl;
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
  //cout << "UpdateLidar starts"<< endl;
   int n_z_ = 2;
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z_, 2 * n_aug_ + 1);

/*******************************************************************************
 * Student part begin
 ******************************************************************************/

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    //double v  = Xsig_pred_(2,i);
    //double yaw = Xsig_pred_(3,i);

    //double v1 = cos(yaw)*v;
    //double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0,i) = p_x;                        //r
    Zsig(1,i) = p_y;                                 //phi
    //Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //innovation covariance matrix S
  MatrixXd S = MatrixXd(n_z_,n_z_);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z_,n_z_);
  R <<    std_laspx_*std_laspx_, 0,
          0, std_laspy_*std_laspy_;
          
  S = S + R;
  
  VectorXd z = meas_package.raw_measurements_;
  
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_);

/*******************************************************************************
 * Student part begin
 ******************************************************************************/

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();

  
  //cout << "UpdateLidar ends"<< endl;
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
  //cout << "UpdateRadar starts"<< endl;
    int n_z_ = 3;
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z_, 2 * n_aug_ + 1);

/*******************************************************************************
 * Student part begin
 ******************************************************************************/

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
//cout<< Xsig_pred_.size()<< endl;
    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);
//cout<<"here2"<< endl;
    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig(1,i) = atan2(p_y,p_x);                                 //phi
    Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }
//cout<< "flag 2"<< endl;
  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //innovation covariance matrix S
  MatrixXd S = MatrixXd(n_z_,n_z_);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z_,n_z_);
  R <<    std_radr_*std_radr_, 0, 0,
          0, std_radphi_*std_radphi_, 0,
          0, 0,std_radrd_*std_radrd_;
  S = S + R;
  
  VectorXd z = meas_package.raw_measurements_;
  
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_);

/*******************************************************************************
 * Student part begin
 ******************************************************************************/

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();
  
  
 // cout << "UpdateRadar ends"<< endl;
}
