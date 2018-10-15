#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include "tools.h"

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
  P_ << 1,0,0,0,0,
  0,1,0,0,0,
  0,0,1,0,0,
  0,0,0,1,0,
  0,0,0,0,1;
  

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 3;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 3;
  
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
  
  n_x_=5;
  
  n_aug_=7;
  Xsig_pred_ = MatrixXd(n_x_,2*n_aug_+1);

  
  lambda_=3-n_x_;
  
  //set measurement dimension, radar can measure r, phi, and r_dot
  n_z = 3;
  
  
  // initializing matrices
  R_Laser = MatrixXd(2, 2);
  H_ = MatrixXd(2,5);

  //measurement covariance matrix - laser
  R_Laser << std_laspx_*std_laspx_	,	 0	,
						0			, std_laspy_ * std_laspy_;

  H_ << 1, 0, 0, 0,0,
		0, 1, 0, 0,0;
  // Create vector for weights
  weights_ = VectorXd(2 * n_aug_ + 1);
  
  //add measurement noise covariance matrix
  R_Radar = MatrixXd(n_z,n_z);
  R_Radar <<    std_radr_*std_radr_		,			0				,			0			,
						0				, std_radphi_*std_radphi_	,			0			,
						0				,			0				,std_radrd_*std_radrd_	;
  ///* the current NIS for radar
  NIS_radar_ = 0.0;

  ///* the current NIS for laser
  NIS_laser_ = 0.0;
  
}

UKF::~UKF() {}


/**********************************************************************************
							UKF::ProcessMeasurement
 **********************************************************************************/

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
  
  double 	previous_timestamp_ = 0.0;
  Tools tools;


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
    cout << "UKF: " << endl;
    x_ << 1, 1, 0,0,0;

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */

	  double rho = meas_package.raw_measurements_[0]; // range
	  double phi = meas_package.raw_measurements_[1]; // bearing
	  double rho_dot = meas_package.raw_measurements_[2]; // velocity of rho

	  double px = rho * cos(phi); // position x
	  double py = rho * sin(phi); // position y
	  double vx = rho_dot * cos(phi); // velocity x
	  double vy = rho_dot * sin(phi); // velocity y
	  double v = sqrt(vx * vx + vy * vy); // velocity
	
	  x_ << px, py, v, 0, 0;
		 
	  
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
	   x_(0)= meas_package.raw_measurements_(0) ;
	   x_(1)= meas_package.raw_measurements_(1);
    }
	
	previous_timestamp_ = meas_package.timestamp_;
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }
  
  //compute the time elapsed between the current and previous measurements
	float dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
	previous_timestamp_ = meas_package.timestamp_;
	// Call prediction
	cout << "Prediction started." << endl;
	Prediction(dt);
	cout << "Prediction completed." << endl;
	
	
	if(meas_package.sensor_type_==MeasurementPackage::LASER){
		cout << "UpdateLidar started." << endl;
		UpdateLidar(meas_package);
		cout << "UpdateLidar Finished." << endl;
	}
	else if(meas_package.sensor_type_==MeasurementPackage::RADAR){
		cout << "UpdateRadar started." << endl;
		UpdateRadar(meas_package);
		cout << "UpdateRadar Finished." << endl;
	}
	
  
}


/**********************************************************************************
***********************************************************************************
							UKF::Prediction
 **********************************************************************************
 **********************************************************************************/

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
  
  
  /*********************Creating Sigma Points *************************************/
  
  
  //create augmented mean vector
  VectorXd x_aug = VectorXd(7);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(7, 7);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, (2 * n_aug_ + 1));
  
  //create augmented mean state
  x_aug<<x_,0,0;
  //create augmented covariance matrix
  MatrixXd Q = MatrixXd(2,2);
  Q<<std_a_,0,
  0,std_yawdd_;
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;
  //create square root matrix
  MatrixXd A = MatrixXd(7,7);
  A = P_aug.llt().matrixL();
  //std::cout<<"A is :"<<A<<std::endl;
  //create augmented sigma points
  Xsig_aug.col(0)=x_aug;
  double lambda_constant = sqrt(lambda_+n_aug_);
  //std::cout<<"The value of Lambda constant is :"<<lambda_constant<<std::endl;
  
  for(int i=1;i<=n_aug_;i++)
  {
      Xsig_aug.col(i)=x_aug + (lambda_constant*A.col(i-1));
      Xsig_aug.col(i+n_aug_)=x_aug - (lambda_constant*A.col(i-1));
  }
  
  /******************************************************************************
						PREDICT SIGMA POINTS 	
	*****************************************************************************/
  for (int i = 0; i< 2*n_aug_+1; i++)
  {
    //extract values for better readability
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }
  
  /******************************************************************************
						PREDICT MEAN AND COVARIANCE
	*****************************************************************************/
  // set weights
  double weight_0 = lambda_/(lambda_+n_aug_);
  weights_(0) = weight_0;
  double weight=0.0;
  for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
    weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
  }

  //predicted state mean
  //x.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    x_ = x_+ weights_(i) * Xsig_pred_.col(i);
  }

  //predicted state covariance matrix
  //P.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = VectorXd(5);
	Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
  }
}


/**********************************************************************************
 **********************************************************************************
							UKF::UpdateLidar
 **********************************************************************************
 **********************************************************************************/

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
  
		
	//Creating Z matrix
	VectorXd z =VectorXd(2);
	z<< meas_package.raw_measurements_(0), 
	meas_package.raw_measurements_(1);
  /* KF Measurement update step*/
		MatrixXd I = MatrixXd::Identity(5, 5);
		VectorXd Y = VectorXd(2);
		Y = z - H_ * x_;
		MatrixXd Ht = MatrixXd(5,2);
		Ht=H_.transpose();
		MatrixXd S = MatrixXd(2,2);
		S = H_ * P_ * Ht + R_Laser;
		MatrixXd K =  MatrixXd(5,5);
		K = P_ * Ht * S.inverse();

		//new state
		x_ = x_ + (K * Y);
		P_ = (I - K * H_) * P_;
		
		// Calculate NIS update
	NIS_laser_ = Y.transpose() * S.inverse() * Y;
  
}


 
 /**********************************************************************************
  **********************************************************************************
							UKF::UpdateRadar
  **********************************************************************************
  **********************************************************************************/
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
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  
  
  //create matrix for sigma points in measurement space
  VectorXd z = VectorXd(n_z);
  
  //innovation covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  
  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  
  
  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig(1,i) = atan2(p_y,p_x);                                 //phi
    Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }
  cout<<"Sigma points are transformed"<<endl;

  //calculate mean predicted measurement
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //calculate innovation covariance matrix S
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = VectorXd(3);
	z_diff=Zsig.col(i) - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  // Add measurement noise to covariance matrix
  S = S + R_Radar;
  
  cout<<"Measurement noise added"<<endl;
  
  
  
  
  /******************************************************************************
						UPDATE RADAR 	
	*****************************************************************************/
  
  z = meas_package.raw_measurements_;
  cout<<"Z : "<<z<<endl;
  
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff1 = VectorXd(3);
	z_diff1=Zsig.col(i) - z_pred;
	cout<<"Z_diff : "<<z_diff1<<endl;
    //angle normalization
    while (z_diff1(1)> M_PI) z_diff1(1)-=2.*M_PI;
    while (z_diff1(1)<-M_PI) z_diff1(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = VectorXd(5);
	x_diff=Xsig_pred_.col(i) - x_;
	cout<<"X_diff : "<<x_diff<<endl;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff1.transpose();
  }
  
  cout<<"TC : "<<Tc<<endl;

  //Kalman gain K;
  MatrixXd K = MatrixXd(n_x_,n_z);
  K=Tc * S.inverse();

  //residual
  VectorXd z_res = VectorXd(3);
  z_res=z - z_pred;

  //angle normalization
  while (z_res(1)> M_PI) z_res(1)-=2.*M_PI;
  while (z_res(1)<-M_PI) z_res(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K * z_res;
  P_ = P_ - K*S*K.transpose();
  
  // Calculate NIS update
  NIS_radar_ = z_res.transpose() * S.inverse() * z_res;
  
  
}
