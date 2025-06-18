/* ----------------------------------------------------------------------------
 * Copyright 2018, Ross Hartley <m.ross.hartley@gmail.com>
 * All Rights Reserved
 * See LICENSE for the license information
 * -------------------------------------------------------------------------- */

/**
 *  @file   kinematics.cpp
 *  @author Ross Hartley
 *  @brief  Example of invariant filtering for contact-aided inertial navigation
 *  @date   September 25, 2018
 **/


#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <Eigen/Dense>
#include <boost/algorithm/string.hpp>
#include <vector>
#include "InEKF.h"

#define DT_MIN 1e-6
#define DT_MAX 1

using namespace std;
using namespace inekf;

double stod98(const std::string &s) {
    return atof(s.c_str());
}

int stoi98(const std::string &s) {
    return atoi(s.c_str());
}


void SavePoseAsTUM(const std::string& filename, Eigen::Matrix3d orient, Eigen::Vector3d tran, double t) {
  static bool first_call = true;
  std::ofstream save_points;

  if (first_call) {
    save_points.open(filename, std::ios::trunc);
    first_call = false;
  } else {
    save_points.open(filename, std::ios::app);
  }

  if (!save_points) {
    std::cerr << "Error opening file: " << filename << std::endl;
    return;
  }

  save_points.setf(std::ios::fixed, std::ios::floatfield);
  Eigen::Quaternion<double> q(orient);

  save_points.precision(9);
  save_points << t << " ";
  save_points.precision(10);
  save_points << tran(0) << " "
              << tran(1) << " "
              << tran(2) << " "
              << q.x() << " "
              << q.y() << " "
              << q.z() << " "
              << q.w() << std::endl;
}

int main() {
    //  ---- Initialize invariant extended Kalman filter ----- //
    RobotState initial_state; 

    // Initialize state mean
    Eigen::Matrix3d R0;
    Eigen::Vector3d v0, p0, bg0, ba0;
    // R0 << 1, 0, 0, // initial orientation
    //       0, 1, 0, // IMU frame is rotated 90deg about the x-axis
    //       0, 0, 1;
    R0 << 
    0.9999185963458515, -2.304557019609009e-05,    0.01275931631809694,
    0.01275933713025751,   0.001806025993930204,    -0.9999169653468294;
                     0,    -0.9999983688681978,  -0.001806173023014302,
    v0 << 0,0,0; // initial velocity
    p0 << 0,0,0; // initial position
    // bg0 << 0,0,0; // initial gyroscope bias
    // ba0 << 0,0,0; // initial accelerometer bias
    bg0 << -1.81176e-05, 5.99727e-05, -0.000257213;
    ba0 << 0.000218875, 2.95973e-05, -0.0169909;
    initial_state.setRotation(R0);
    initial_state.setVelocity(v0);
    initial_state.setPosition(p0);
    initial_state.setGyroscopeBias(bg0);
    initial_state.setAccelerometerBias(ba0);

    // Initialize state covariance
    NoiseParams noise_params;
    noise_params.setGyroscopeNoise(0.001);
    noise_params.setAccelerometerNoise(0.01);
    noise_params.setGyroscopeBiasNoise(0.0001);
    noise_params.setAccelerometerBiasNoise(0.0001);
    noise_params.setContactNoise(0.01);

    // Initialize filter
    InEKF filter(initial_state, noise_params);
    cout << "Noise parameters are initialized to: \n";
    cout << filter.getNoiseParams() << endl;
    cout << "Robot's state is initialized to: \n";
    cout << filter.getState() << endl;

    // Open data file
    ifstream infile("../src/data/imu_kinematic_measurements_cear.txt");
    string line;
    Eigen::Matrix<double,6,1> imu_measurement = Eigen::Matrix<double,6,1>::Zero();
    Eigen::Matrix<double,6,1> imu_measurement_prev = Eigen::Matrix<double,6,1>::Zero();
    double t = 0;
    double t_prev = 0;

    // int NUM_LINES = 638;
    // int count = 0;

    // ---- Loop through data file and read in measurements line by line ---- //
    while (getline(infile, line)){
        vector<string> measurement;
        boost::split(measurement,line,boost::is_any_of(" "));
        // // Handle measurements
        if (measurement[0].compare("IMU")==0){
            // cout << "Received IMU Data, propagating state\n";
            assert((measurement.size()-2) == 6);
            t = stod98(measurement[1]); 
            // Read in IMU data
            imu_measurement << stod98(measurement[2]), 
                               stod98(measurement[3]), 
                               stod98(measurement[4]),
                               stod98(measurement[5]),
                               stod98(measurement[6]),
                               stod98(measurement[7]);

            // Propagate using IMU data
            double dt = t - t_prev;
            // if (dt > DT_MIN && dt < DT_MAX) {
                filter.Propagate(imu_measurement_prev, dt);
            // }
            // std::cout << "imu_measurement: " << imu_measurement.transpose() << "\n";

            // Store previous timestamp
            t_prev = t;
            imu_measurement_prev = imu_measurement;
        }
        else if (measurement[0].compare("CONTACT")==0){
            // cout << "Received CONTACT Data, setting filter's contact state\n";
            assert((measurement.size()-2)%2 == 0);
            vector<pair<int,bool> > contacts;
            int id;
            bool indicator;
            // t = stod98(measurement[1]); 
            // Read in contact data
            for (int i=2; i<measurement.size(); i+=2) {
                id = stoi98(measurement[i]);
                indicator = bool(stod98(measurement[i+1]));
                contacts.push_back(pair<int,bool> (id, indicator));
                // std::cout << "contact id: " << id << ", indicator: " << indicator << "\n";
            }
            // Set filter's contact state
            filter.setContacts(contacts);
        }
        else if (measurement[0].compare("KINEMATIC")==0){
            // cout << "Received KINEMATIC observation, correcting state\n";  
            assert((measurement.size()-2)%4 == 0);
            int id;
            Eigen::Quaternion<double> q;
            Eigen::Vector3d p;
            Eigen::Matrix4d pose = Eigen::Matrix4d::Identity();
            Eigen::Matrix<double,6,6> covariance;
            vectorKinematics measured_kinematics;
            double t_meas = stod98(measurement[1]); 
            // Read in kinematic data
            for (int i=2; i<measurement.size(); i+=4) {
                id = stoi98(measurement[i]);
                q = Eigen::Quaternion<double>(1, 0, 0, 0);
                q.normalize();
                p << stod98(measurement[i+1]),stod98(measurement[i+2]),stod98(measurement[i+3]);
                pose.block<3,3>(0,0) = q.toRotationMatrix();
                pose.block<3,1>(0,3) = p;
                covariance = Eigen::Matrix<double,6,6>::Identity() * 1e-2;
                Kinematics frame(id, pose, covariance);
                measured_kinematics.push_back(frame);
                // std::cout << "pose, id: " << id << "\n" << pose << "\n";
            }
            // cout << "measured_kinematics.size(): " << measured_kinematics.size() << "\n";
            // Correct state using kinematic measurements
            filter.CorrectKinematics(measured_kinematics);
            static std::string pose_path = "/home/s/data/cear/indoor/mocap1_well-lit_trot/algo_pose/InEKF.txt";
            SavePoseAsTUM(pose_path, filter.getState().getRotation(), filter.getState().getPosition(), t_meas);
        }

        // count++;
        // if (count > NUM_LINES)
        //     break;
    }

    // Print final state
    cout.precision(17);
    // cout << "imu:\n " << imu_measurement << endl;
    cout << "t:\n " << t << endl;

    cout << filter.getState() << endl;
    cout << "Covariance: \n" << filter.getState().getP() << endl;
    return 0;
}
