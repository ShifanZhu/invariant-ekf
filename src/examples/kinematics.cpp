/* ----------------------------------------------------------------------------
 * Copyright 2018, Ross Hartley
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
#include <Eigen/Dense>
#include <boost/algorithm/string.hpp>
#include <vector>
#include <tuple>
#include "InEKF.h"

#define DT_MIN 1e-6
#define DT_MAX 1

using namespace std;
using namespace inekf;

int main() {
    //  ---- Initialize invariant extended Kalman filter ----- //
    RobotState initial_state; 

    // Initialize state mean
    Eigen::Matrix3d R0;
    Eigen::Vector3d v0, p0, bg0, ba0;
    R0 << 1, 0, 0, // initial orientation
          0, -1, 0, // IMU frame is rotated 90deg about the x-axis
          0, 0, -1;
    v0 << 0,0,0; // initial velocity
    p0 << 0,0,0; // initial position
    bg0 << 0,0,0; // initial gyroscope bias
    ba0 << 0,0,0; // initial accelerometer bias
    initial_state.setRotation(R0);
    initial_state.setVelocity(v0);
    initial_state.setPosition(p0);
    initial_state.setGyroscopeBias(bg0);
    initial_state.setAccelerometerBias(ba0);

    // Initialize state covariance
    NoiseParams noise_params;
    noise_params.setGyroscopeNoise(0.01);
    noise_params.setAccelerometerNoise(0.1);
    noise_params.setGyroscopeBiasNoise(0.00001);
    noise_params.setAccelerometerBiasNoise(0.0001);
    noise_params.setContactNoise(0.01);

    // Initialize filter
    InEKF filter(initial_state, noise_params);
    cout << "Noise parameters are initialized to: \n";
    cout << filter.getNoiseParams() << endl;
    cout << "Robot's state is initialized to: \n";
    cout << filter.getState() << endl;

    // Open data file
    ifstream infile("../src/data/imu_kinematic_measurements.txt");
    string line;
    Eigen::Matrix<double,6,1> imu_measurement = Eigen::Matrix<double,6,1>::Zero();
    Eigen::Matrix<double,6,1> imu_measurement_prev = Eigen::Matrix<double,6,1>::Zero();
    double t = 0;
    double t_prev = 0;

    // ---- Loop through data file and read in measurements line by line ---- //
    while (getline(infile, line)){
        vector<string> measurement;
        boost::split(measurement,line,boost::is_any_of(" "));
        // // Handle measurements
        if (measurement[0].compare("IMU")==0){
            cout << "Received IMU Data, propagating state\n";
            assert((measurement.size()-2) == 6);
            t = stod(measurement[1]); 
            // Read in IMU data
            imu_measurement << stod(measurement[2]), 
                               stod(measurement[3]), 
                               stod(measurement[4]),
                               stod(measurement[5]),
                               stod(measurement[6]),
                               stod(measurement[7]);

            // Propagate using IMU data
            double dt = t - t_prev;
            if (dt > DT_MIN && dt < DT_MAX) {
                filter.Propagate(imu_measurement_prev, dt);
            }

        }
        else if (measurement[0].compare("CONTACT")==0){
            cout << "Received CONTACT Data, setting filter's contact state\n";
            assert((measurement.size()-2)%2 == 0);
            vector<pair<int,bool>> contacts;
            int id;
            bool indicator;
            t = stod(measurement[1]); 
            // Read in contact data
            for (int i=2; i<measurement.size(); i+=2) {
                id = stoi(measurement[i]);
                indicator = bool(stod(measurement[i+1]));
                contacts.push_back(pair<int,bool> (id, indicator));
            }       
            // Set filter's contact state
            filter.setContacts(contacts);
        }
        else if (measurement[0].compare("KINEMATIC")==0){
            cout << "Received KINEMATIC observation, correcting state\n";  
            assert((measurement.size()-2)%44 == 0);
            int id;
            Eigen::Quaternion<double> q;
            Eigen::Vector3d p;
            Eigen::Matrix4d H = Eigen::Matrix4d::Identity();
            Eigen::Matrix<double,6,6> Cov;
            vectorTupleIntMatrix4dMatrix6d measured_kinematics;
            t = stod(measurement[1]); 
            // Read in kinematic data
            for (int i=2; i<measurement.size(); i+=44) {
                id = stoi(measurement[i]); 
                q = Eigen::Quaternion<double> (stod(measurement[i+1]),stod(measurement[i+2]),stod(measurement[i+3]),stod(measurement[i+4]));
                q.normalize();
                p << stod(measurement[i+5]),stod(measurement[i+6]),stod(measurement[i+7]);
                H.block<3,3>(0,0) = q.toRotationMatrix();
                H.block<3,1>(0,3) = p;
                for (int j=0; j<6; ++j) {
                    for (int k=0; k<6; ++k) {
                        Cov(j,k) = stod(measurement[i+8 + j*6+k]);
                    }
                }
                measured_kinematics.push_back(tuple<int,Eigen::Matrix4d,Eigen::Matrix<double,6,6>> (id, H, Cov));
            }
            // Correct state using kinematic measurements
            filter.CorrectKinematics(measured_kinematics);
        }

        // Store previous timestamp
        t_prev = t;
        imu_measurement_prev = imu_measurement;
    }

    // Print final state
    cout << filter.getState() << endl;
    return 0;
}
