/**
 * \file gait_generator.hpp
 *
 * \date Created on sep 1, 2021
 * \author Pablo Romero Sorozábal <pablo.romero@csic.es>
 */

#ifndef GAIT_GENERATOR__HPP_
#define GAIT_GENERATOR__HPP_

#include <cstdio>
#include <memory>
#include <string>

#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>  // for std::min_element, std::max_element
#include <cassert>    // for assert
#include <cmath>

#include <chrono>
#include <thread>

#include "rclcpp/rclcpp.hpp"
#include "rcutils/cmdline_parser.h"

#include <ament_index_cpp/get_package_share_directory.hpp>

#include <std_msgs/msg/string.hpp>
#include <std_msgs/msg/float32_multi_array.hpp>
#include <std_msgs/msg/float64_multi_array.hpp>
#include <std_msgs/msg/int32_multi_array.hpp>

#include <geometry_msgs/msg/pose.hpp>
#include <geometry_msgs/msg/pose_array.hpp>
#include <geometry_msgs/msg/twist.hpp>
#include <geometry_msgs/msg/transform_stamped.hpp>

#include <tf2/LinearMath/Quaternion.h>
#include <tf2_ros/static_transform_broadcaster.h>
#include <tf2_ros/transform_broadcaster.h>

#include <visualization_msgs/msg/marker_array.hpp>

#include <d2w_interfaces/srv/gait_generator_config.hpp>
#include <d2w_interfaces/msg/pelvis_forces.hpp>
#include <d2w_interfaces/msg/pelvis_motors_states.hpp>
#include <d2w_interfaces/msg/pelvis_pose.hpp>
#include <d2w_interfaces/msg/gait_setpoints.hpp>
#include <d2w_interfaces/msg/treadmill_state.hpp>


#include <Eigen/Dense>
#include <vector>

// Show comments
#define DEBUG 1

#define FORCE_SCALE 0.01
#define MOMENTUM_SCALE 0.02

#define TORADIANS 3.1415/180

using namespace std::chrono_literals;

// Lightweight linear interpolation class to replace libInterpolate
class LightweightInterpolator {
public:
  LightweightInterpolator() = default;
  
  void setData(const std::vector<float>& x, const std::vector<float>& y) {
    if (x.size() != y.size() || x.empty()) {
      return;
    }
    x_data_ = x;
    y_data_ = y;
  }
  
  double operator()(double x) const {
    return interpolate(x);
  }
  
  double interpolate(double x) const {
    // Handle empty data case
    if (x_data_.empty()) {
      return 0.0;
    }
    if (x_data_.size() == 1) {
      return y_data_.front();
    }

    const double dx = x_data_[1] - x_data_[0];
    const double period = x_data_.back() - x_data_.front() + dx;
    if (period > 0.0) {
      x = std::fmod(x - x_data_.front(), period);
      if (x < 0.0) {
        x += period;
      }
      x += x_data_.front();
    }
    
    if (x <= x_data_.front()) {
      return y_data_.front();
    }
    if (x > x_data_.back()) {
      const double alpha = (x - x_data_.back()) / dx;
      return y_data_.back() + alpha * (y_data_.front() - y_data_.back());
    }
    
    // Find the interval where x falls
    size_t i = 0;
    while (i < x_data_.size() - 1 && x > x_data_[i + 1]) {
      ++i;
    }
    
    // Linear interpolation formula: y = y1 + (x - x1) * (y2 - y1) / (x2 - x1)
    double x1 = x_data_[i];
    double x2 = x_data_[i + 1];
    double y1 = y_data_[i];
    double y2 = y_data_[i + 1];
    
    return y1 + (x - x1) * (y2 - y1) / (x2 - x1);
  }
  
private:
  std::vector<float> x_data_;
  std::vector<float> y_data_;
};

class GaitGenerator : public rclcpp::Node {
  public:
    GaitGenerator();
    virtual ~GaitGenerator();

    // Callback current treadmill vel
    void treadmillSystemVelCallback(const d2w_interfaces::msg::TreadmillState::SharedPtr msg);

    // Timer
    void timer_callback();

    // Publis pelvis pose
    void gaitSetpointsPub();

    // System service response
    void gaitGeneratorConfigSRV(const std::shared_ptr<d2w_interfaces::srv::GaitGeneratorConfig::Request> request, std::shared_ptr<d2w_interfaces::srv::GaitGeneratorConfig::Response> response);

    void updateGaitParameters(float patient_height, float gait_velocity);

    void updateGaitPercentage();

    void updateBodySegments(float patient_height);

    void updateGaitGainsXYZ(float step_x_gain, float step_y_gain, float step_z_gain);

    void updateDirectGaitKinematics();

    Eigen::Matrix<float, 3, 3> getRotationMatrix(float roll,float pitch,float yaw);

    void updateSpatioTemporalGaitData();

    void updateAnklePelvisPathsXYZ();

    void updatePaths();

    void smoothGaussian(std::vector<Eigen::Matrix<float, 3, 1>>& data, const int window_size);

    void loadAnglesFromFiles();

    void parametrizeAndRemap(
      const std::vector<float>& y,
      const std::vector<float>& t,
      float new_min,
      float new_max,
      std::vector<float>& y_new,
      float& t0_min,
      float& t0_max,
      float& y0_min,
      float& y0_max
    );

    void computeBounds(
      const std::vector<float>& data,
      float gain,
      float& lo,
      float& hi
    );

    void getReferenceAnglesData(const std::string& file_name, std::vector<float>& angle_x, std::vector<float>& angle_y, std::vector<float>& angle_z);

    // call this after loading your file‐based arrays
    void applySagittalScaling();

    // Scale X & Y components of that trajectory by parametrization
    void applyCartesianScaling();

    // Vector module
    float vectorModule(Eigen::Matrix<float, 3, 1> A);

  private:
    rclcpp::Time time_;
    rclcpp::TimerBase::SharedPtr timer_;

    float gait_percentage_;
    std::vector<float> gait_percentage_vector_;
    float gaint_increment_;
    float l_gaint_index_;
    float r_gaint_index_;
    float delta_t_;
    int count_;

    // Replace libInterpolate with our lightweight interpolator
    LightweightInterpolator interp_;
    LightweightInterpolator interp_pelvis_x_;
    LightweightInterpolator interp_pelvis_y_;
    LightweightInterpolator interp_pelvis_z_;
    LightweightInterpolator interp_pelvis_jointRot_roll_;
    LightweightInterpolator interp_pelvis_jointRot_pitch_;
    LightweightInterpolator interp_pelvis_jointRot_yaw_;
    LightweightInterpolator interp_hip_x_;
    LightweightInterpolator interp_hip_y_;
    LightweightInterpolator interp_hip_z_;
    LightweightInterpolator interp_hip_jointRot_roll_;
    LightweightInterpolator interp_hip_jointRot_pitch_;
    LightweightInterpolator interp_hip_jointRot_yaw_;
    LightweightInterpolator interp_hip_dot_x_;
    LightweightInterpolator interp_hip_dot_y_;
    LightweightInterpolator interp_hip_dot_z_;
    LightweightInterpolator interp_knee_x_;
    LightweightInterpolator interp_knee_y_;
    LightweightInterpolator interp_knee_z_;
    LightweightInterpolator interp_knee_jointRot_roll_;
    LightweightInterpolator interp_knee_jointRot_pitch_;
    LightweightInterpolator interp_knee_jointRot_yaw_;
    LightweightInterpolator interp_knee_dot_x_;
    LightweightInterpolator interp_knee_dot_y_;
    LightweightInterpolator interp_knee_dot_z_;
    LightweightInterpolator interp_ankle_x_;
    LightweightInterpolator interp_ankle_y_;
    LightweightInterpolator interp_ankle_z_;
    LightweightInterpolator interp_ankle_jointRot_roll_;
    LightweightInterpolator interp_ankle_jointRot_pitch_;
    LightweightInterpolator interp_ankle_jointRot_yaw_;
    LightweightInterpolator interp_ankle_dot_x_;
    LightweightInterpolator interp_ankle_dot_y_;
    LightweightInterpolator interp_ankle_dot_z_;
    LightweightInterpolator interp_ankle_ddot_x_;
    LightweightInterpolator interp_ankle_ddot_y_;
    LightweightInterpolator interp_ankle_ddot_z_;

    int toe_off_index_;

    float patient_height_;

    Eigen::Matrix<float, 3, 1> leg_dimentions_;
    Eigen::Matrix<float, 3, 1> pelvis_dimentions_;
    Eigen::Matrix<float, 3, 1> shank_dimentions_;
    Eigen::Matrix<float, 3, 1> thigh_dimentions_;
    Eigen::Matrix<float, 3, 1> foot_dimentions_;
    Eigen::Matrix<float, 3, 1> heel_position_;
    Eigen::Matrix<float, 3, 1> leg_proportion_;
    Eigen::Matrix<float, 3, 1> pelvis_proportion_;
    Eigen::Matrix<float, 3, 1> shank_proportion_;
    Eigen::Matrix<float, 3, 1> thigh_proportion_;
    Eigen::Matrix<float, 3, 1> foot_proportion_;
    Eigen::Matrix<float, 3, 1> heel_proportion_;


    float gait_velocity_;
    float step_length_;
    float step_period_;
    float t_increment_;
    int gait_n_; // Gait Samples

    int control_amplitudes_;
    int control_mode_;
    float hip_amplitude_;
    float knee_amplitude_;
    float ankle_amplitude_;

    float gait_x_gain_;
    float gait_y_gain_;
    float gait_z_gain_;


    // Gains you set elsewhere:
    float ankleTraj_gain_x_{1.0f}, ankleTraj_gain_y_{1.0f};

    float pelvis_jointRot_roll_gain_;
    float pelvis_jointRot_pitch_gain_;
    float pelvis_jointRot_yaw_gain_;
    Eigen::Matrix<float, 3, 3> R_jointRot_pelvis_;
    Eigen::Matrix<float, 3, 3> R_pelvis_absolute_;

    float hip_jointRot_roll_gain_;
    float hip_jointRot_pitch_gain_;
    float hip_jointRot_yaw_gain_;
    Eigen::Matrix<float, 3, 3> R_jointRot_hip_;
    Eigen::Matrix<float, 3, 3> R_hip_absolute_;

    float knee_jointRot_roll_gain_;
    float knee_jointRot_pitch_gain_;
    float knee_jointRot_yaw_gain_;
    Eigen::Matrix<float, 3, 3> R_jointRot_knee_;
    Eigen::Matrix<float, 3, 3> R_knee_absolute_;

    float ankle_jointRot_roll_gain_;
    float ankle_jointRot_pitch_gain_;
    float ankle_jointRot_yaw_gain_;
    Eigen::Matrix<float, 3, 3> R_jointRot_ankle_;
    Eigen::Matrix<float, 3, 3> R_ankle_absolute_;

    // storage for scaled (remapped) sagittal angles
    std::vector<double> hip_pitch_scaled_,
                        knee_pitch_scaled_,
                        ankle_pitch_scaled_,
                        pelvis_yaw_scaled_;

    // extrema and phase locations (optional diagnostics)
    double hip_phaseMin_, hip_phaseMax_, hip_valMin_, hip_valMax_;
    double knee_phaseMin_, knee_phaseMax_, knee_valMin_, knee_valMax_;
    double ankle_phaseMin_, ankle_phaseMax_, ankle_valMin_, ankle_valMax_;
    double psi_phaseMin_, psi_phaseMax_, psi_valMin_, psi_valMax_;


    std::vector<float> pelvis_jointRot_roll_;
    std::vector<float> pelvis_jointRot_pitch_;
    std::vector<float> pelvis_jointRot_yaw_;
    std::vector<float> pelvis_jointRot_roll_file_;
    std::vector<float> pelvis_jointRot_pitch_file_;
    std::vector<float> pelvis_jointRot_yaw_file_;

    std::vector<float> hip_jointRot_roll_;
    std::vector<float> hip_jointRot_pitch_;
    std::vector<float> hip_jointRot_yaw_;
    std::vector<float> hip_jointRot_roll_file_;
    std::vector<float> hip_jointRot_pitch_file_;
    std::vector<float> hip_jointRot_yaw_file_;

    std::vector<float> knee_jointRot_roll_;
    std::vector<float> knee_jointRot_pitch_;
    std::vector<float> knee_jointRot_yaw_;
    std::vector<float> knee_jointRot_roll_file_;
    std::vector<float> knee_jointRot_pitch_file_;
    std::vector<float> knee_jointRot_yaw_file_;

    std::vector<float> ankle_jointRot_roll_;
    std::vector<float> ankle_jointRot_pitch_;
    std::vector<float> ankle_jointRot_yaw_;
    std::vector<float> ankle_jointRot_roll_file_;
    std::vector<float> ankle_jointRot_pitch_file_;
    std::vector<float> ankle_jointRot_yaw_file_;

    std::vector<float> pelvis_jointRot_roll_absolute_;
    std::vector<float> pelvis_jointRot_pitch_absolute_;
    std::vector<float> pelvis_jointRot_yaw_absolute_;

    std::vector<float> hip_jointRot_roll_absolute_;
    std::vector<float> hip_jointRot_pitch_absolute_;
    std::vector<float> hip_jointRot_yaw_absolute_;

    std::vector<float> knee_jointRot_roll_absolute_;
    std::vector<float> knee_jointRot_pitch_absolute_;
    std::vector<float> knee_jointRot_yaw_absolute_;

    std::vector<float> ankle_jointRot_roll_absolute_;
    std::vector<float> ankle_jointRot_pitch_absolute_;
    std::vector<float> ankle_jointRot_yaw_absolute_;

    std::vector<Eigen::Matrix<float, 3, 1>> K_p_absolute_xyz_;

    std::vector<Eigen::Matrix<float, 3, 1>> hip_relative_xyz_;
    std::vector<Eigen::Matrix<float, 3, 1>> knee_relative_xyz_;
    std::vector<Eigen::Matrix<float, 3, 1>> ankle_relative_xyz_;
    std::vector<Eigen::Matrix<float, 3, 1>> front_foot_relative_xyz_;
    std::vector<Eigen::Matrix<float, 3, 1>> heel_foot_relative_xyz_;

    std::vector<Eigen::Matrix<float, 3, 1>> hip_absolute_xyz_;
    std::vector<Eigen::Matrix<float, 3, 1>> knee_absolute_xyz_;
    std::vector<Eigen::Matrix<float, 3, 1>> ankle_absolute_xyz_;
    std::vector<Eigen::Matrix<float, 3, 1>> front_foot_absolute_xyz_;
    std::vector<Eigen::Matrix<float, 3, 1>> heel_foot_absolute_xyz_;

    std::vector<std::vector<float>> pelvis_path_xyz_;
    std::vector<std::vector<float>> hip_path_xyz_;
    std::vector<std::vector<float>> knee_path_xyz_;
    std::vector<std::vector<float>> ankle_path_xyz_;

    std::vector<std::vector<float>> hip_path_xyz_dot_;
    std::vector<std::vector<float>> knee_path_xyz_dot_;
    std::vector<std::vector<float>> ankle_path_xyz_dot_;
    std::vector<std::vector<float>> ankle_path_xyz_ddot_;

    std::vector<Eigen::Matrix<float, 3, 1>> r_hip_relative_xyz_;
    std::vector<Eigen::Matrix<float, 3, 1>> r_knee_relative_xyz_;
    std::vector<Eigen::Matrix<float, 3, 1>> r_ankle_relative_xyz_;
    std::vector<Eigen::Matrix<float, 3, 1>> r_front_foot_relative_xyz_;
    std::vector<Eigen::Matrix<float, 3, 1>> r_heel_foot_relative_xyz_;

    std::vector<Eigen::Matrix<float, 3, 1>> l_hip_relative_xyz_;
    std::vector<Eigen::Matrix<float, 3, 1>> l_knee_relative_xyz_;
    std::vector<Eigen::Matrix<float, 3, 1>> l_ankle_relative_xyz_;
    std::vector<Eigen::Matrix<float, 3, 1>> l_front_foot_relative_xyz_;
    std::vector<Eigen::Matrix<float, 3, 1>> l_heel_foot_relative_xyz_;

    Eigen::Matrix<float, 3, 1> pelvis_position_ref_;
    Eigen::Matrix<float, 3, 1> pelvis_orientation_ref_;

    Eigen::Matrix<float, 3, 1> r_hip_position_ref_;
    Eigen::Matrix<float, 3, 1> r_hip_orientation_ref_;
    Eigen::Matrix<float, 3, 1> l_hip_position_ref_;
    Eigen::Matrix<float, 3, 1> l_hip_orientation_ref_;
    Eigen::Matrix<float, 3, 1> r_dot_hip_position_ref_;
    Eigen::Matrix<float, 3, 1> l_dot_hip_position_ref_;

    Eigen::Matrix<float, 3, 1> r_knee_position_ref_;
    Eigen::Matrix<float, 3, 1> r_knee_orientation_ref_;
    Eigen::Matrix<float, 3, 1> l_knee_position_ref_;
    Eigen::Matrix<float, 3, 1> l_knee_orientation_ref_;
    Eigen::Matrix<float, 3, 1> r_dot_knee_position_ref_;
    Eigen::Matrix<float, 3, 1> l_dot_knee_position_ref_;

    Eigen::Matrix<float, 3, 1> r_ankle_position_ref_;
    Eigen::Matrix<float, 3, 1> r_ankle_orientation_ref_;
    Eigen::Matrix<float, 3, 1> l_ankle_position_ref_;
    Eigen::Matrix<float, 3, 1> l_ankle_orientation_ref_;
    Eigen::Matrix<float, 3, 1> r_dot_ankle_position_ref_;
    Eigen::Matrix<float, 3, 1> l_dot_ankle_position_ref_;
    Eigen::Matrix<float, 3, 1> r_ddot_ankle_position_ref_;
    Eigen::Matrix<float, 3, 1> l_ddot_ankle_position_ref_;

    rclcpp::Publisher<d2w_interfaces::msg::GaitSetpoints>::SharedPtr gait_setpoints_pub_; // Cables Tension publisher

    rclcpp::Subscription<d2w_interfaces::msg::TreadmillState>::SharedPtr treasmill_system_state_;    // Pelvis system positions subscriber
    rclcpp::Subscription<d2w_interfaces::msg::PelvisForces>::SharedPtr pelvis_system_forces_;  // Pelvis cables system forces subscriber
    // Services
    rclcpp::Service<d2w_interfaces::srv::GaitGeneratorConfig>::SharedPtr gait_generator_config_srv_;

    Eigen::Matrix<float, 3, 1> K_O_;   // World reference
    Eigen::Matrix<float, 3, 1> K_p_;   // Pelvis reference (position respect world)
    Eigen::Matrix<float, 3, 1> r_;     // World and pelvis reference relation

    Eigen::Matrix<float, 3, 3> rotationMatrix_ref_;
    Eigen::Matrix<float, 3, 3> rotationMatrix_real_;

    std::vector<Eigen::Matrix<float, 3, 1>>  frame_position_;                // Pelvis frame position
    std::vector<Eigen::Matrix<float, 3, 1>>  pelvis_vertex_position_;        // Vertex pelvis frame position
    geometry_msgs::msg::Pose pelvis_position_real_;

    std::vector<Eigen::Matrix<float, 3, 1>> cables_vector_ref_;
    std::vector<Eigen::Matrix<float, 3, 1>> cables_vector_real_;

    float cables_length_ref_[4];
    float cables_length_real_[4];

    Eigen::Matrix<float, 4, 1> F_r_ref_;     // Desired Force
    Eigen::Matrix<float, 4, 1> F_r_real_; // Current Force
    Eigen::Matrix<float, 3, 1> F_n_ref_;     // Desired Net Force
    Eigen::Matrix<float, 3, 1> F_n_real_; // Current Force
    Eigen::Matrix<float, 3, 1> M_n_ref_;     // Desired Net Momentum
    Eigen::Matrix<float, 3, 1> M_n_real_; // Current Force


    std::vector<std::vector<double>> joints_b_matrix_hip_x = {
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {1, 100, 8.8646, 25.2285, 63.3460, 66.6465, 0.0720, 0.0710, 0.0606, 0.0676, 0.0156, 0.0002},
        {0, 0, 5.3787, -8.6110, -21.8429, -3.8830, 0.0334, 0.0348, 0.0310, 0.0106, 0.0154, 0.0437},
        {0, 0, -1.9211, 2.9896, 7.3256, -0.2061, -0.0057, -0.0061, -0.0067, -0.0003, -0.0179, -0.0226},
        {0, 0, -6.0606, -5.0543, 7.8859, 15.3166, -0.0913, -0.0919, -0.0733, -0.0684, -0.0439, -0.0246}
    };

    std::vector<std::vector<double>> joints_b_matrix_hip_y = {
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {1, 100, 5.9468, 31.4805, 45.6573, 65.1239, 83.1363, 95.3359, -0.1928, -0.1933, -0.1888, -0.1840, -0.1823, -0.2012, -0.2095, -0.2104},
        {0, 0, 11.2694, -5.5060, -12.8955, -4.9300, 1.7166, -7.2327, 0.0059, 0.0064, -0.0055, -0.0073, -0.0111, -0.0135, 0.0094, 0.0056},
        {0, 0, -3.0070, 3.4434, 6.5973, 2.9452, -0.3467, 3.5352, -0.0018, -0.0020, 0.0003, 0.0012, 0.0028, 0.0056, -0.0013, -0.0010},
        {0, 0, -3.9045, -5.8303, -3.3970, -12.3372, -19.0547, -5.3547, 0, 0.0001, 0.0029, 0.0037, 0.0086, 0.0272, 0.0227, 0.0162}
    };

    std::vector<std::vector<double>> joints_b_matrix_hip_z = {
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {1, 11.7, 37.0154, 45.8875, 54.7139, 78.3177, 86.5573, 96.3023, 100, -0.0643, -0.0433, -0.0545, -0.0402, -0.0351, -0.0582, -0.0553, -0.0522, -0.0636},
        {0, -3.4047, -1.6027, -8.0697, 0.7883, -17.0165, -6.1448, 1.8475, 0, -0.0073, 0.0069, 0.0073, 0.0006, 0.0142, -0.0158, -0.0115, -0.0039, -0.0083},
        {0, 1.0898, -0.2568, 3.3916, -1.1818, 7.2416, 2.9204, -0.1356, 0, 0.0029, 0.0002, -0.0021, 0.0034, -0.0018, 0.0022, 0.0011, -0.0004, 0.0033},
        {0, 2.8548, -2.1529, -1.2473, -2.9856, -0.9579, -5.7899, -11.4876, 0, -0.1040, -0.1334, -0.1266, -0.1386, -0.1472, -0.1142, -0.1141, -0.1188, -0.1045}
    };

    std::vector<std::vector<double>> joints_b_matrix_knee_x = {
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {1, 100, 3, 105.7534, 62.3725, 31.5462, 79.9849, 0.0208, 0.0247, 0.0082, 0.0480, -0.0646, -0.0592, 0.0248},
        {0, 0, 0, -12.8941, -8.7701, -2.9197, -6.0040, 0.1046, 0.0995, 0.1158, 0.0918, -0.0300, 0.0379, 0.0416},
        {0, 0, 0, 8.2491, 2.5472, 0.8830, 1.2487, -0.0199, -0.0182, -0.0229, -0.0205, -0.0154, -0.0215, -0.0334},
        {0, 0, 0, -14.4223, -3.1504, -0.5634, -9.6609, 0.0555, 0.0549, 0.0609, 0.0537, 0.0057, 0.0646, -0.0285}
    };

    std::vector<std::vector<double>> joints_b_matrix_knee_y = {
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {1, 100, 35.3133, 64.2554, 40.9716, -0.1323, -0.1337, -0.1018, -0.1560, -0.1848},
        {0, 0, 1.5127, -5.9552, 14.6156, 0.0234, 0.0242, -0.0262, 0.0087, 0.0286},
        {0, 0, -1.5240, -0.6017, -7.6236, -0.0070, -0.0071, 0.0069, -0.0018, -0.0058},
        {0, 0, 0.7373, 2.9130, 38.5794, -0.0361, -0.0359, -0.0056, -0.0026, -0.0087}
    };

    std::vector<std::vector<double>> joints_b_matrix_knee_z = {
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {1, 100, 3.6658, 7.8639, 17.7545, 58.1538, 89.7454, 25.5931, 74.4786, 68.4865, -0.0610, -0.0592, -0.0691, -0.0601, -0.0486, 0.0016, -0.0506, -0.0487, -0.0245, -0.0151},
        {0, 0, 3.2288, 6.0583, 1.1279, -2.3664, -2.5824, -4.2946, -5.9368, -2.7077, 0.0056, 0.0039, 0.0158, 0.0200, 0.0132, 0.0127, 0.0155, -0.0004, -0.0167, -0.0095},
        {0, 0, -0.8010, -1.4178, -0.3918, 0.1448, 2.1660, 0.8356, 2.2208, 0.5881, 0.0042, 0.0048, 0.0004, 0.0020, -0.0018, 0.0027, -0.0041, -0.0020, 0.0033, 0.0050},
        {0, 0, -2.9177, -6.8727, 0.9488, -2.3826, -3.8450, 10.9049, -5.8527, -6.9611, -0.5294, -0.5301, -0.5268, -0.5405, -0.5596, -0.6109, -0.5315, -0.5565, -0.5715, -0.5841}
    };

    std::vector<std::vector<double>> joints_b_matrix_ankle_x = {
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {1, 100, 94.1282, 67.8849, 34.0867, 80.8377, 0.0688, 0.0696, 0.0876, 0.0219, 0.0318, -0.0486},
        {0, 0, 2.0552, -11.1419, -5.7542, 8.7105, 0.1880, 0.1888, 0.1865, -0.1133, 0.0379, 0.2926},
        {0, 0, -0.5700, 2.8409, 1.4942, -2.9021, -0.0399, -0.0411, -0.0412, -0.0031, -0.0279, -0.0881},
        {0, 0, 1.1641, 0.0859, 0.2834, -0.7210, 0.0372, 0.0416, 0.0466, -0.2813, -0.1129, -0.0702}
    };

    std::vector<std::vector<double>> joints_b_matrix_ankle_y = {
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {1, 100, 25.6709, 71.4560, 97.0288, 12.7549, 87.2338, -0.1490, -0.1518, -0.0547, -0.1796, -0.1669, -0.0805, -0.1512},
        {0, 0, 2.1959, -7.0469, 4.2048, 1.3383, -2.1021, 0.0456, 0.0468, -0.0342, 0.0425, 0.0502, -0.0143, 0.0563},
        {0, 0, -0.0796, 0.2736, -1.1993, -0.1675, -0.2167, -0.0136, -0.0137, 0.0062, -0.0127, -0.0136, 0.0009, -0.0125},
        {0, 0, 10.5549, 2.2769, -6.9066, 5.5219, -4.9807, -0.0106, -0.0096, 0.0012, 0.0214, -0.0118, -0.0011, -0.0215}
    };

    std::vector<std::vector<double>> joints_b_matrix_ankle_z = {
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {1, 100, 16.4272, 31.8741, 66.4170, 72.5524, 74.4771, 91.4997, -0.0188, -0.0177, -0.0291, -0.0406, -0.0703, -0.0769, -0.0322, -0.0220},
        {0, 0, -0.8210, -1.2131, -4.7483, -4.5847, 4.3339, -1.7623, 0.0316, 0.0325, 0.0191, 0.0007, 0.1171, 0.1107, 0.0341, 0.0121},
        {0, 0, 0.2948, 0.4183, -1.2272, -0.0272, -2.7957, 0.0867, -0.0045, -0.0050, 0.0013, -0.0012, -0.0355, -0.0228, -0.0036, -0.0035},
        {0, 0, -1.6838, -3.0703, 2.7522, 1.9151, 4.8419, 0.0645, -1.0498, -1.0504, -1.0484, -1.0263, -0.9421, -0.9035, -0.9719, -1.0280}
    };

    // New subscriber for adaptive oscillator's gait cycle percentage
    rclcpp::Subscription<std_msgs::msg::Float64MultiArray>::SharedPtr adaptive_gait_cycle_sub_;
  
    // Variables to store adaptive gait cycle percentages
    double adaptive_gait_cycle_left_;
    double adaptive_gait_cycle_right_;
    double adaptive_estimated_freq_left_;
    double adaptive_estimated_freq_right_;

    // New callback function for adaptive gait cycle percentage
    void adaptiveGaitCycleCallback(const std_msgs::msg::Float64MultiArray::SharedPtr msg);

    // Add these new member variables
    double ankle_x_offset_;      // Cached ankle x offset
    float last_gait_percentage_; // Last processed gait percentage
    float last_r_gaint_index_;   // Last processed right gait index
    float last_l_gaint_index_;   // Last processed left gait index;
};

#endif  // GAIT_GENERATOR_HPP_
