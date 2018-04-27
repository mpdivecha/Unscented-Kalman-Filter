#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

#define PI 3.14159265

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF()
{
    // if this is false, laser measurements will be ignored (except during init)
    use_laser_ = true;

    // if this is false, radar measurements will be ignored (except during init)
    use_radar_ = true;

    // initial state vector
    x_ = VectorXd(5);

    // initial covariance matrix
    P_ = MatrixXd(5, 5);
    P_.setZero();

    // Process noise standard deviation longitudinal acceleration in m/s^2
    std_a_ = 2;

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
    n_x_ = 5;   // State vector dimension
    n_aug_ = 7; // Augmented state dimension
    lambda_ = 3 - n_aug_;
    Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
    Xsig_pred_.setZero();
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package)
{
    /**
    TODO:

    Complete this function! Make sure you switch between lidar and radar
    measurements.
    */
    if (!is_initialized_)
    {
        /**
         TODO:
        *  Initialize the state x_ with the first measurement
        *  Create the covariance matrix
        *  Remember: you'll need to convert radar from polar to cartesian coordinates
        */
        // First measurement
        cout << "EKF: " << endl;

        previous_timestamp_ = meas_package.timestamp_;

        if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
        {
            /**
            * Convert radar from polar to cartesian coordinates and initialize state
            */
            float rho = meas_package.raw_measurements_[0];
            float phi = meas_package.raw_measurements_[1];
            //float rho_dot = meas_package.raw_measurements_[2];

            float px = rho * cos(phi);
            float py = rho * sin(phi);
            x_ << px, py, 0, 0, 0;
        }
        else if (meas_package.sensor_type_ == MeasurementPackage::LASER)
        {
            /*
            Initialize state
            */
            x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
        }
        //create example covariance matrix
        // This example is taken from the exercise code
        MatrixXd P = MatrixXd(n_x_, n_x_);
        P << 0.0043, -0.0013, 0.0030, -0.0022, -0.0020,
            -0.0013, 0.0077, 0.0011, 0.0071, 0.0060,
            0.0030, 0.0011, 0.0054, 0.0007, 0.0008,
            -0.0022, 0.0071, 0.0007, 0.0098, 0.0100,
            -0.0020, 0.0060, 0.0008, 0.0100, 0.0123;
        
        //P_ = P;
        P_.diagonal() = VectorXd::Ones(5);
        P_(0, 0) = 1.5;
        P_(1, 1) = 1.5;

        cout << "Intialized" << endl;
        // done initializing, no need to predict or update
        is_initialized_ = true;
        return;
    }

    /*****************************************************************************
    * Prediction
    *****************************************************************************/
    cout << "Going to predict" << endl;

    //compute the time elapsed between the current and previous measurements
    double dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0; //dt - expressed in seconds
    previous_timestamp_ = meas_package.timestamp_;

    Prediction(dt);

    /*****************************************************************************
    * Update
    *****************************************************************************/
    if (use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {
        cout << "Going to update Radar" << endl;
        UpdateRadar(meas_package);
    }
    else if (use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER)
    {
        cout << "Going to update Laser" << endl;
        UpdateLidar(meas_package);
    }
}

/* Generate sigma points
 * @param {MatrixXd*} Xsig_out   The generated sigma points
 */
void UKF::AugmentedSigmaPoints(MatrixXd *Xsig_out)
{
    // Create augmented covariance matrix
    MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
    P_aug.setZero();

    // Define the spreading parameter
    double lambda = 3 - n_aug_;

    // Create sigma point matrix
    //MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

    /*
    TODO:
    * Create augmented mean state
    * Create augmented covariance matrix
    * Create square root matrix
    * Create augmented sigma points
    */
    //cout << "Inside AugmentedSigmaPoints, generating sigma pts" << endl;
    int nExtraDims = n_aug_ - n_x_;

    // Create augmented mean vector
    VectorXd x_aug_(n_aug_);
    x_aug_.setZero();
    //cout << "Going to set the x_aug_ vector";
    //cout << "x_aug_.size: " << x_aug_.size() << " x_.size: " << x_.size() << endl;
    x_aug_.head(n_x_) = x_;

    // The Process noise covariance matrix Q
    //cout << "Getting matrix Q" << endl;
    MatrixXd Q(nExtraDims, nExtraDims);
    Q << std_a_ * std_a_, 0,
        0, std_yawdd_ * std_yawdd_;

    //cout << "Filling matrix P_aug" << endl;
    P_aug.topLeftCorner(n_x_, n_x_) = P_;
    P_aug.bottomRightCorner(nExtraDims, nExtraDims) = Q;

    // The square root of P_aug
    //cout << "Getting sqrt matrix A" << endl;
    MatrixXd A = P_aug.llt().matrixL();
    MatrixXd times = sqrt(lambda + n_aug_) * A;

    //cout << "Filling out Xsig_out" << endl;
    Xsig_out->col(0) = x_aug_;
    Xsig_out->block(0, 1, n_aug_, n_aug_) = times.colwise() + x_aug_;
    Xsig_out->block(0, n_aug_ + 1, n_aug_, n_aug_) = ((-(times.colwise() - x_aug_).array()+1.0)-1.0);   // -(times.colwise() - x_aug);

    //cout << "\nXsig_out: " << *Xsig_out << endl << endl;
    //cout << "x_: " << x_ << endl << endl << "x_aug_: " << x_aug_ << endl;
    //*Xsig_out = Xsig_aug;
}

/* Predicts sigma points
 * @param {MatrixXd} Xsig_aug The input augmented sigma points
 * @param {double} delta_t Time elapsed between previous measurement, in seconds.
 * @param {MatrixXd*} Xsig_out The predicted sigma points matrix
 */
void UKF::SigmaPointPrediction(MatrixXd Xsig_aug, double delta_t, MatrixXd *Xsig_out)
{
    //create matrix with predicted sigma points as columns
    //MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);

    // Predict sigma points
    for (int i = 0; i < (2 * n_aug_ + 1); i++)
    {
        // sv is the state_vec. Using abbreviation for brevity
        VectorXd sv = Xsig_aug.col(i).head(n_x_);
        // The remaining augmented vector
        VectorXd ag = Xsig_aug.col(i).tail(n_aug_ - n_x_);

        // Some short-cuts
        double sin_psi = sin(sv(3));             // sin of psi
        double cos_psi = cos(sv(3));             // cos of psi
        double psidot_dt = sv(4) * delta_t;      // psi_dot times delta_t
        double nu_psi = sv(2) / sv(4);           // nu divided by psi_dot
        double dtsq_2 = (delta_t * delta_t) / 2; // delta_t squared divided by 2

        if (sv(4) == 0) // If psi_dot is zero
        {
            Xsig_pred_(0, i) = sv(0) + sv(2) * cos_psi * delta_t + dtsq_2 * cos_psi * ag(0);
            Xsig_pred_(1, i) = sv(1) + sv(2) * sin_psi * delta_t + dtsq_2 * sin_psi * ag(0);
        }
        else
        {
            Xsig_pred_(0, i) = sv(0) + nu_psi * (sin(sv(3) + psidot_dt) - sin_psi) + dtsq_2 * cos_psi * ag(0);
            Xsig_pred_(1, i) = sv(1) + nu_psi * (-cos(sv(3) + psidot_dt) + cos_psi) + dtsq_2 * sin_psi * ag(0);
        }
        Xsig_pred_(2, i) = sv(2) + 0 + delta_t * ag(0);
        Xsig_pred_(3, i) = sv(3) + psidot_dt + dtsq_2 * ag(1);
        Xsig_pred_(4, i) = sv(4) + 0 + delta_t * ag(1);
    }
    //*Xsig_out = Xsig_pred;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t)
{
    /**
   TODO:

   Complete this function! Estimate the object's location. Modify the state
   vector, x_. Predict sigma points, the state, and the state covariance matrix.
   */
    // The augmented sigma points
    //cout << "Going to generate  sigma points" << endl;
    MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
    AugmentedSigmaPoints(&Xsig_aug);

    //cout << "\nXsig_aug: " << Xsig_aug << endl << endl;
    //cout << "x_: " << x_ << endl;

    // The predicted sigma points after the non-linear transformation
    //MatrixXd Xsig_pred;
    //cout << "Going to predict sigma points" << endl;
    SigmaPointPrediction(Xsig_aug, delta_t, &Xsig_pred_);

    // Compute the mean state and covariance

    // Set the weights
    double vec_len = 2 * n_aug_ + 1;
    VectorXd weights(vec_len);
    //cout << "lambda_: " << lambda_ << endl;
    weights(0) = lambda_ / (lambda_ + n_aug_);
    for (int i = 1; i < vec_len; i++)
    {
        weights(i) = 1 / (2 * (lambda_ + n_aug_));
    }

    // Predicted mean state from the sigma points
    x_ = Xsig_pred_ * weights;

    // Predicted covariance matrix from the sigma points
    MatrixXd residual = Xsig_pred_.array().colwise() - x_.array();
    MatrixXd wResidual = residual.array().rowwise() * weights.transpose().array();
    P_ = wResidual * residual.transpose();
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package)
{
    /**
    TODO:

    Complete this function! Use lidar data to update the belief about the object's
    position. Modify the state vector, x_, and covariance, P_.

    You'll also need to calculate the lidar NIS.
    */
   // Predict the measurement

   int n_z = 2;     // Measurement dimension for Laser

   VectorXd z = meas_package.raw_measurements_;

   // Set the weights
   int vec_len = 2 * n_aug_ + 1;
   VectorXd weights(vec_len);
   weights(0) = lambda_ / (lambda_ + n_aug_);
   for (int i = 1; i < vec_len; i++)
   {
       weights(i) = 1 / (2 * (lambda_ + n_aug_));
   }

   // Create matrix for sigma points in measurement space
   MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

   // Mean predicted measurement
   VectorXd z_pred(n_z);

   // Measurement covariance matrix S
   MatrixXd S(n_z, n_z);

   // Transform the sigma points into measurement space
   int n_pts = 2 * n_aug_ + 1;

   for (int i = 0; i < n_pts; i++)
   {
       VectorXd sig_pt = Xsig_pred_.col(i);
       double denom = sig_pt(0)*sig_pt(0) + sig_pt(1)*sig_pt(1);
       Zsig(0, i) = sig_pt(0); //sqrt(denom);
       Zsig(1, i) = sig_pt(1); //atan2(sig_pt(1), sig_pt(0));
   }

   //cout << "Going to compute z_pred" << endl;
   z_pred = Zsig * weights;

  // cout << "Going to compute R" << endl;
   MatrixXd R(2, 2);
   R.setZero();
   R(0, 0) = std_laspx_ * std_laspx_;
   R(1, 1) = std_laspy_ * std_laspy_;

   //cout << "Going to compute residual" << endl;
   MatrixXd residual = Zsig.array().colwise() - z_pred.array();
   MatrixXd wResidual = residual.array().rowwise() * weights.transpose().array();
   //cout << "Going to compute S" << endl;
   S = wResidual * residual.transpose() + R;
   MatrixXd Sinv = S.inverse();

   //cout << "Going to compute Tc" << endl;
   MatrixXd Tc = MatrixXd(n_x_, n_z);

    // The cross-correlation matrix
   Tc = (Xsig_pred_.array().colwise() - x_.array());
   Tc = Tc.array().rowwise() * weights.transpose().array();
   Tc = Tc * MatrixXd((Zsig.array().colwise() - z.array()).transpose());

   //cout << "Going to compute K" << endl;
   // The Kalman gain matrix
   MatrixXd K = Tc * Sinv;

   x_ = x_ + K*(z - z_pred);

   P_ = P_ - K*S*K.transpose();

   VectorXd zdiff = (z - z_pred);
    //cout << "zdiff: " << zdiff << " S.inv: " << Sinv << endl;
    double nis = zdiff.transpose() * Sinv * zdiff;
    cout << "Calculated nis is " << nis << endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package)
{
    /**
    TODO:

    Complete this function! Use radar data to update the belief about the object's
    position. Modify the state vector, x_, and covariance, P_.

    You'll also need to calculate the radar NIS.
    */
    // Predict the measurement

    int n_z = 3; // Measurement dimension for Radar

    VectorXd z = meas_package.raw_measurements_;

    // Set the weights
    double vec_len = 2 * n_aug_ + 1;
    VectorXd weights(vec_len);
    weights(0) = lambda_ / (lambda_ + n_aug_);
    for (int i = 1; i < vec_len; i++)
    {
        weights(i) = 1 / (2 * (lambda_ + n_aug_));
    }

    // Create matrix for sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

    // Mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);

    // Measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z, n_z);

    // Transform the sigma points into measurement space
    int n_pts = 2 * n_aug_ + 1;

    //cout << "\nXsig_pred_: " << Xsig_pred_ << endl;

    for (int i = 0; i < n_pts; i++)
    {
        VectorXd sig_pt = Xsig_pred_.col(i);
        double px = sig_pt(0);
        double py = sig_pt(1);
        double v =  sig_pt(2);
        double psi = sig_pt(3);
        double denom = sqrt(px*px + py*py); //sig_pt(0) * sig_pt(0) + sig_pt(1) * sig_pt(1);
        Zsig(0, i) = denom; //sqrt(denom);
        Zsig(1, i) = atan2(py, px); //atan2(sig_pt(1), sig_pt(0));
        Zsig(2, i) = (px*cos(psi)*v + py*sin(psi)*v)/denom; //(sig_pt(0) * cos(sig_pt(3)) * sig_pt(2) + sig_pt(1) * sin(sig_pt(3)) * sig_pt(2)) / sqrt(denom);
    }

    //cout << "\nweights: " << weights.transpose() << endl;
    //cout << "\nZsig: " << Zsig << endl;

    //cout << "Going to compute z_pred" << endl;
    z_pred = Zsig * weights;
    if (z_pred(1) <= -PI)
    {
        z_pred(1) = z_pred(1) + 2 * PI;
    }
    else if (z_pred(1) >= PI)
    {
        z_pred(1) = z_pred(1) - 2 * PI;
    }
    //cout << "\nz_pred: " << z_pred.transpose() << " z: " << z.transpose() << endl << endl;


    //cout << "Going to compute R" << endl;
    MatrixXd R(3, 3);
    R.setZero();
    R(0, 0) = std_radr_ * std_radr_;
    R(1, 1) = std_radphi_ * std_radphi_;
    R(2, 2) = std_radrd_ * std_radrd_;

    //cout << "Going to compute residual" << endl;
    MatrixXd residual = Zsig.array().colwise() - z_pred.array(); //z - z_pred;
    //cout << "Going to compute weighted residual with r.size: " << residual.size() << " ,w.size(): " << weights.size() << endl;
    MatrixXd wResidual = residual.array().rowwise() * weights.transpose().array();  // residual * weights;
    //cout << "Going to compute S" << endl;
    S = wResidual * residual.transpose() + R;
    //cout << "S: " << S << endl;
    MatrixXd Sinv = S.inverse();

    MatrixXd Tc = MatrixXd(n_x_, n_z);

    //cout << "Going to compute Tc" << endl;
    Tc = (Xsig_pred_.array().colwise() - x_.array());
    Tc = Tc.array().rowwise() * weights.transpose().array();
    Tc = Tc * MatrixXd((Zsig.array().colwise() - z.array()).transpose());

    // The Kalman gain matrix
    //cout << "Going to compute K" << endl;
    MatrixXd K = Tc * Sinv;

    x_ = x_ + K * (z - z_pred);

    P_ = P_ - K * S * K.transpose();

    VectorXd zdiff = (z - z_pred);
    //cout << "zdiff: " << zdiff << " S.inv: " << Sinv << endl;
    double nis = zdiff.transpose() * Sinv * zdiff;
    cout << "Calculated nis is " << nis << endl;
}
