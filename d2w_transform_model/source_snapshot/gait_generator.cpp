/**
 * \file gait_generator.cpp
 *
 * \date Created on sep 1, 2021
 * \author Pablo Romero Sorozábal <pablo.romero@csic.es>
 */

#include <cstdio>
#include <cmath>
#include <memory>
#include <string>
#include <iostream>

#include <math.h>
#include <Eigen/Geometry>

#include "rclcpp/rclcpp.hpp"
#include "rcutils/cmdline_parser.h"
#include <tf2/LinearMath/Quaternion.h>
#include <tf2/LinearMath/Matrix3x3.h>
#include <tf2_geometry_msgs/tf2_geometry_msgs.hpp>


#include "gait_generator/gait_generator.hpp"

using Eigen::MatrixXd;


GaitGenerator::GaitGenerator() : Node("gait_generator") {
  RCLCPP_INFO(this->get_logger(), "[GAIT GENERATOR] Run pelvis control node");

  // QoS ROS2 settings:
  int8_t qos_depth = 1;
  this->declare_parameter("qos_depth", 10);
  this->get_parameter("qos_depth", qos_depth);
  const auto QOS_RKL10V = rclcpp::QoS(rclcpp::KeepLast(qos_depth)).best_effort().durability_volatile();


  // Initilize xyz positions
  Eigen::Matrix<float, 3, 1> zero_xyz;
  zero_xyz << 0,0,0;

  control_mode_ = 0;
  control_amplitudes_ = 0;
  hip_amplitude_ = 0.30F;
  knee_amplitude_ = 0.01F;
  ankle_amplitude_ = 0.20F;

  std::vector<float> zero_vector_float;
  for (int i = 0; i < 100; i++)
    zero_vector_float.push_back(0);

  for (int i = 0; i < 100; i++) {
    // Get direct kinematics of each joint
    hip_relative_xyz_.push_back(zero_xyz);
    knee_relative_xyz_.push_back(zero_xyz);
    ankle_relative_xyz_.push_back(zero_xyz);
    front_foot_relative_xyz_.push_back(zero_xyz);
    heel_foot_relative_xyz_.push_back(zero_xyz);
    // Get absolut position
    K_p_absolute_xyz_.push_back(zero_xyz);
    hip_absolute_xyz_.push_back(zero_xyz);
    knee_absolute_xyz_.push_back(zero_xyz);
    ankle_absolute_xyz_.push_back(zero_xyz);
    front_foot_absolute_xyz_.push_back(zero_xyz);
    heel_foot_absolute_xyz_.push_back(zero_xyz);
    //
    pelvis_jointRot_roll_.push_back(0);
    pelvis_jointRot_pitch_.push_back(0);
    pelvis_jointRot_yaw_.push_back(0);
    pelvis_jointRot_roll_file_.push_back(0);
    pelvis_jointRot_pitch_file_.push_back(0);
    pelvis_jointRot_yaw_file_.push_back(0);
    pelvis_jointRot_roll_absolute_.push_back(0);
    pelvis_jointRot_pitch_absolute_.push_back(0);
    pelvis_jointRot_yaw_absolute_.push_back(0);

    hip_jointRot_roll_.push_back(0);
    hip_jointRot_pitch_.push_back(0);
    hip_jointRot_yaw_.push_back(0);
    hip_jointRot_roll_file_.push_back(0);
    hip_jointRot_pitch_file_.push_back(0);
    hip_jointRot_yaw_file_.push_back(0);
    hip_jointRot_roll_absolute_.push_back(0);
    hip_jointRot_pitch_absolute_.push_back(0);
    hip_jointRot_yaw_absolute_.push_back(0);

    knee_jointRot_roll_.push_back(0);
    knee_jointRot_pitch_.push_back(0);
    knee_jointRot_yaw_.push_back(0);
    knee_jointRot_roll_file_.push_back(0);
    knee_jointRot_pitch_file_.push_back(0);
    knee_jointRot_yaw_file_.push_back(0);
    knee_jointRot_roll_absolute_.push_back(0);
    knee_jointRot_pitch_absolute_.push_back(0);
    knee_jointRot_yaw_absolute_.push_back(0);

    ankle_jointRot_roll_.push_back(0);
    ankle_jointRot_pitch_.push_back(0);
    ankle_jointRot_yaw_.push_back(0);
    ankle_jointRot_roll_file_.push_back(0);
    ankle_jointRot_pitch_file_.push_back(0);
    ankle_jointRot_yaw_file_.push_back(0);
    ankle_jointRot_roll_absolute_.push_back(0);
    ankle_jointRot_pitch_absolute_.push_back(0);
    ankle_jointRot_yaw_absolute_.push_back(0);

    pelvis_path_xyz_.push_back(zero_vector_float);
    hip_path_xyz_.push_back(zero_vector_float);
    knee_path_xyz_.push_back(zero_vector_float);
    ankle_path_xyz_.push_back(zero_vector_float);

    hip_path_xyz_dot_.push_back(zero_vector_float);
    knee_path_xyz_dot_.push_back(zero_vector_float);
    ankle_path_xyz_dot_.push_back(zero_vector_float);
    ankle_path_xyz_ddot_.push_back(zero_vector_float);  // Initialize ankle acceleration vector
  }

  // Start in gait percentage 0%
  gait_percentage_ = 0; // 0-100%
  gaint_increment_ = 0;
  step_length_ = 0;
  count_ = 0;
  gait_velocity_ = 0;
  gait_n_ = 100;
  step_period_ = 0.0F;
  t_increment_ = 0.0F;
  for (int i = 0; i < 100; i++)
    gait_percentage_vector_.push_back(i);


  // Load reference angles
  RCLCPP_INFO(this->get_logger(), "[GAIT GENERATOR] Get angles from files");
  loadAnglesFromFiles();

  // Set gait gains
  RCLCPP_INFO(this->get_logger(), "[GAIT GENERATOR] Update Gait Gains XYZ");
  updateGaitGainsXYZ(0.001,0.001,0.001);

  // Kinematics initialization
  // Reference systems
  K_O_ << 0,0,0.15; // World reference origin
  K_p_ << 0,0,0.44; // Platform reference origin (start in the same position)

  // Inital patient gait:
  patient_height_ = 1.2;

  // Segments relative dimentions to height:
  leg_proportion_     <<  0         ,0        ,-0.530;
  pelvis_proportion_  <<  0         ,0.191/2  ,-0.01 ;
  shank_proportion_   <<  0         ,0        ,-0.246;
  thigh_proportion_   <<  0         ,0        ,-0.246;
  foot_proportion_    << -0.122/1.5 ,0        ,-0.05;
  heel_proportion_      <<  0.03/1.5  ,0        ,-0.05;
 
  // Update body segments
  RCLCPP_INFO(this->get_logger(), "[GAIT GENERATOR] Update body segments");
  updateBodySegments(patient_height_);

  // Update initial gait kinematics
  RCLCPP_INFO(this->get_logger(), "[GAIT GENERATOR] Update direct gait kinematics");
  updateDirectGaitKinematics();
 
  // Update initial temporal gait data
  RCLCPP_INFO(this->get_logger(), "[GAIT GENERATOR] Get spatio temporal gait data");
  updateSpatioTemporalGaitData();

  // Update initial paths
  RCLCPP_INFO(this->get_logger(), "[GAIT GENERATOR] Update ankle pelvis paths XYZ");
  updateAnklePelvisPathsXYZ();


  RCLCPP_INFO(this->get_logger(), "[GAIT GENERATOR] Initilize Gait Generator");
  // Initilize subscriptors:
  treasmill_system_state_ = this->create_subscription<d2w_interfaces::msg::TreadmillState>("/d2w/treadmill/state", QOS_RKL10V, std::bind(&GaitGenerator::treadmillSystemVelCallback, this, std::placeholders::_1));

  // New subscriber for adaptive oscillator's gait cycle percentage
  adaptive_gait_cycle_sub_ = this->create_subscription<std_msgs::msg::Float64MultiArray>(
      "/d2w/gait_generator/adapted_gait_cycle", 10, std::bind(&GaitGenerator::adaptiveGaitCycleCallback, this, std::placeholders::_1));

  // Initilize publishers:
  gait_setpoints_pub_ = this->create_publisher<d2w_interfaces::msg::GaitSetpoints>("/d2w/gait_generator/setpoints", 10);
  
  // Services
  gait_generator_config_srv_ = this->create_service<d2w_interfaces::srv::GaitGeneratorConfig>("/d2w/gait_generator/config", std::bind(&GaitGenerator::gaitGeneratorConfigSRV, this, std::placeholders::_1, std::placeholders::_2));
  
  // Initilize timers:
  delta_t_ = 1/ 100.0; // [s]
  auto timer_period = std::chrono::duration<float, std::milli>(delta_t_ * 1000);
  timer_ = this->create_wall_timer(timer_period, std::bind(&GaitGenerator::timer_callback, this));
  
  // Initialize adaptive gait cycle percentages
  adaptive_gait_cycle_left_ = 0.0;
  adaptive_gait_cycle_right_ = 0.0;
  adaptive_estimated_freq_left_ = 0.5;
  adaptive_estimated_freq_right_ = 0.5;
}

GaitGenerator::~GaitGenerator() {
}


void GaitGenerator::treadmillSystemVelCallback(const d2w_interfaces::msg::TreadmillState::SharedPtr msg) {
  const float update_threshold = 0.01F; // [m/s]
  const float measured_speed = std::abs(msg->real_speed) < update_threshold
      ? 0.0F
      : msg->real_speed;

  if (std::abs(gait_velocity_ - measured_speed) <= update_threshold) {
    return;
  }

  gait_velocity_ = measured_speed;
  RCLCPP_INFO(this->get_logger(), "[GAIT GENERATOR] Treadmill speed updated: %.3f m/s", gait_velocity_);
  if (control_amplitudes_ == 1) {
    updateSpatioTemporalGaitData();
  } else {
    updateGaitParameters(patient_height_, gait_velocity_);
  }
}

void GaitGenerator::adaptiveGaitCycleCallback(const std_msgs::msg::Float64MultiArray::SharedPtr msg) {
  if (msg->data.size() < 4) {
    RCLCPP_WARN_THROTTLE(this->get_logger(), *this->get_clock(), 2000,
      "Ignoring adapted gait cycle message with %zu values; expected at least 4.",
      msg->data.size());
    return;
  }
  adaptive_gait_cycle_left_ = msg->data[0];
  adaptive_gait_cycle_right_ = msg->data[1];
  adaptive_estimated_freq_left_ = msg->data[2];
  adaptive_estimated_freq_right_ = msg->data[3];
}

void GaitGenerator::gaitGeneratorConfigSRV(const std::shared_ptr<d2w_interfaces::srv::GaitGeneratorConfig::Request> request, std::shared_ptr<d2w_interfaces::srv::GaitGeneratorConfig::Response> response) {
  auto finite_or = [this](double value, float fallback, const char* name) -> float {
    if (!std::isfinite(value)) {
      RCLCPP_WARN(
        this->get_logger(),
        "[GAIT GENERATOR] Ignoring non-finite %s from config request; keeping %.3f.",
        name, fallback);
      return fallback;
    }
    return static_cast<float>(value);
  };

  control_mode_ = static_cast<int>(request->control_generator);
  if (control_mode_ < 0 || control_mode_ > 2) {
    RCLCPP_WARN(
      this->get_logger(),
      "[GAIT GENERATOR] Unsupported control_generator=%ld; using standard mode 0.",
      request->control_generator);
    control_mode_ = 0;
  }

  control_amplitudes_ = static_cast<int>(request->control_amplitudes);
  if (control_amplitudes_ < 0 || control_amplitudes_ > 2) {
    RCLCPP_WARN(
      this->get_logger(),
      "[GAIT GENERATOR] Unsupported control_amplitudes=%ld; using mode 0.",
      request->control_amplitudes);
    control_amplitudes_ = 0;
  }

  hip_amplitude_ = finite_or(request->hip_amplitude, hip_amplitude_, "hip_amplitude");
  knee_amplitude_ = finite_or(request->knee_amplitude, knee_amplitude_, "knee_amplitude");
  ankle_amplitude_ = finite_or(request->ankle_amplitude, ankle_amplitude_, "ankle_amplitude");

  if (hip_amplitude_ <= 1.0e-4F) {
    RCLCPP_WARN(this->get_logger(), "[GAIT GENERATOR] hip_amplitude %.4f is invalid; using 0.30 m.", hip_amplitude_);
    hip_amplitude_ = 0.30F;
  }
  if (knee_amplitude_ < 0.0F) {
    RCLCPP_WARN(this->get_logger(), "[GAIT GENERATOR] knee_amplitude %.4f is below zero; using 0.01 m.", knee_amplitude_);
    knee_amplitude_ = 0.01F;
  }
  if (ankle_amplitude_ <= 0.0F) {
    RCLCPP_WARN(this->get_logger(), "[GAIT GENERATOR] ankle_amplitude %.4f is invalid; using 0.20 m.", ankle_amplitude_);
    ankle_amplitude_ = 0.20F;
  }
  if (knee_amplitude_ > ankle_amplitude_) {
    RCLCPP_WARN(
      this->get_logger(),
      "[GAIT GENERATOR] knee_amplitude %.4f is above ankle_amplitude %.4f; swapping vertical bounds.",
      knee_amplitude_, ankle_amplitude_);
    std::swap(knee_amplitude_, ankle_amplitude_);
  }

  const float requested_height = finite_or(request->patient_height, patient_height_, "patient_height");
  const float safe_height = requested_height > 0.3F && requested_height < 2.5F
      ? requested_height
      : patient_height_;
  if (safe_height != requested_height) {
    RCLCPP_WARN(
      this->get_logger(),
      "[GAIT GENERATOR] patient_height %.3f m is outside the valid range; keeping %.3f m.",
      requested_height, safe_height);
  }

  float safe_velocity = finite_or(request->gait_velocity, gait_velocity_, "gait_velocity");
  if (safe_velocity < 0.0F) {
    RCLCPP_WARN(this->get_logger(), "[GAIT GENERATOR] Negative gait_velocity %.3f m/s; using 0.", safe_velocity);
    safe_velocity = 0.0F;
  }

  if (control_mode_ >= 1)
    updateGaitParameters(safe_height, safe_velocity);
  response->subsystem_ok = 1;
}


void GaitGenerator::gaitSetpointsPub() {

  auto gait_setpoints_msg = d2w_interfaces::msg::GaitSetpoints();

  tf2::Quaternion q;
  q.setRPY(pelvis_orientation_ref_(0), pelvis_orientation_ref_(1), pelvis_orientation_ref_(2));

  gait_setpoints_msg.header.frame_id = "Gait Generator Setpoints";
  gait_setpoints_msg.header.stamp = rclcpp::Clock().now();

  gait_setpoints_msg.r_gait_percentage.data = r_gaint_index_;
  gait_setpoints_msg.l_gait_percentage.data = l_gaint_index_;

  // Publish position and orientation setpoints
  // Pelvis position and orientation
  gait_setpoints_msg.pelvis_ref.position.x = pelvis_position_ref_(0);
  gait_setpoints_msg.pelvis_ref.position.y = pelvis_position_ref_(1);
  gait_setpoints_msg.pelvis_ref.position.z = pelvis_position_ref_(2);
  gait_setpoints_msg.pelvis_ref.orientation = tf2::toMsg(q);

  // Left hip position and orientation
  gait_setpoints_msg.l_hip_ref.position.x = l_hip_position_ref_(0);
  gait_setpoints_msg.l_hip_ref.position.y = l_hip_position_ref_(1);
  gait_setpoints_msg.l_hip_ref.position.z = l_hip_position_ref_(2);
  // Transform the absolute matrix rotation into a quaternion
  q.setRPY(l_hip_orientation_ref_(0), l_hip_orientation_ref_(1), l_hip_orientation_ref_(2));
  gait_setpoints_msg.l_hip_ref.orientation = tf2::toMsg(q);
  // Left hip velocity
  gait_setpoints_msg.l_hip_dot_ref.linear.x = l_dot_hip_position_ref_(0);
  gait_setpoints_msg.l_hip_dot_ref.linear.y = l_dot_hip_position_ref_(1);
  gait_setpoints_msg.l_hip_dot_ref.linear.z = l_dot_hip_position_ref_(2);

  // Right hip position and orientation
  gait_setpoints_msg.r_hip_ref.position.x = r_hip_position_ref_(0);
  gait_setpoints_msg.r_hip_ref.position.y = r_hip_position_ref_(1);
  gait_setpoints_msg.r_hip_ref.position.z = r_hip_position_ref_(2);
  q.setRPY(r_hip_orientation_ref_(0), r_hip_orientation_ref_(1), r_hip_orientation_ref_(2));
  gait_setpoints_msg.r_hip_ref.orientation = tf2::toMsg(q);
  // Right hip velocity
  gait_setpoints_msg.r_hip_dot_ref.linear.x = r_dot_hip_position_ref_(0);
  gait_setpoints_msg.r_hip_dot_ref.linear.y = r_dot_hip_position_ref_(1);
  gait_setpoints_msg.r_hip_dot_ref.linear.z = r_dot_hip_position_ref_(2);

  // Left knee position and orientation
  gait_setpoints_msg.l_knee_ref.position.x = l_knee_position_ref_(0);
  gait_setpoints_msg.l_knee_ref.position.y = l_knee_position_ref_(1);
  gait_setpoints_msg.l_knee_ref.position.z = l_knee_position_ref_(2);
  q.setRPY(l_knee_orientation_ref_(0), l_knee_orientation_ref_(1), l_knee_orientation_ref_(2));
  gait_setpoints_msg.l_knee_ref.orientation = tf2::toMsg(q);
  // Left knee velocity
  gait_setpoints_msg.l_knee_dot_ref.linear.x = l_dot_knee_position_ref_(0);
  gait_setpoints_msg.l_knee_dot_ref.linear.y = l_dot_knee_position_ref_(1);
  gait_setpoints_msg.l_knee_dot_ref.linear.z = l_dot_knee_position_ref_(2);

  // Right knee position and orientation
  gait_setpoints_msg.r_knee_ref.position.x = r_knee_position_ref_(0);
  gait_setpoints_msg.r_knee_ref.position.y = r_knee_position_ref_(1);
  gait_setpoints_msg.r_knee_ref.position.z = r_knee_position_ref_(2);
  q.setRPY(r_knee_orientation_ref_(0), r_knee_orientation_ref_(1), r_knee_orientation_ref_(2));
  gait_setpoints_msg.r_knee_ref.orientation = tf2::toMsg(q);
  // Right knee velocity
  gait_setpoints_msg.r_knee_dot_ref.linear.x = r_dot_knee_position_ref_(0);
  gait_setpoints_msg.r_knee_dot_ref.linear.y = r_dot_knee_position_ref_(1);
  gait_setpoints_msg.r_knee_dot_ref.linear.z = r_dot_knee_position_ref_(2);

  // Left ankle position and orientation
  gait_setpoints_msg.l_ankle_ref.position.x = l_ankle_position_ref_(0);
  gait_setpoints_msg.l_ankle_ref.position.y = l_ankle_position_ref_(1);
  gait_setpoints_msg.l_ankle_ref.position.z = l_ankle_position_ref_(2);
  q.setRPY(l_ankle_orientation_ref_(0), l_ankle_orientation_ref_(1), l_ankle_orientation_ref_(2));
  gait_setpoints_msg.l_ankle_ref.orientation = tf2::toMsg(q);
  // Left ankle velocity
  gait_setpoints_msg.l_ankle_dot_ref.linear.x = l_dot_ankle_position_ref_(0);
  gait_setpoints_msg.l_ankle_dot_ref.linear.y = l_dot_ankle_position_ref_(1);
  gait_setpoints_msg.l_ankle_dot_ref.linear.z = l_dot_ankle_position_ref_(2);

  // Right ankle position and orientation
  gait_setpoints_msg.r_ankle_ref.position.x = r_ankle_position_ref_(0);
  gait_setpoints_msg.r_ankle_ref.position.y = r_ankle_position_ref_(1);
  gait_setpoints_msg.r_ankle_ref.position.z = r_ankle_position_ref_(2);
  q.setRPY(r_ankle_orientation_ref_(0), r_ankle_orientation_ref_(1), r_ankle_orientation_ref_(2));
  gait_setpoints_msg.r_ankle_ref.orientation = tf2::toMsg(q);

  // Right ankle velocity
  gait_setpoints_msg.r_ankle_dot_ref.linear.x = r_dot_ankle_position_ref_(0);
  gait_setpoints_msg.r_ankle_dot_ref.linear.y = r_dot_ankle_position_ref_(1);
  gait_setpoints_msg.r_ankle_dot_ref.linear.z = r_dot_ankle_position_ref_(2);

  // Add ankle accelerations to the message
  // Right ankle acceleration
  gait_setpoints_msg.r_ankle_ddot_ref.linear.x = r_ddot_ankle_position_ref_(0);
  gait_setpoints_msg.r_ankle_ddot_ref.linear.y = r_ddot_ankle_position_ref_(1);
  gait_setpoints_msg.r_ankle_ddot_ref.linear.z = r_ddot_ankle_position_ref_(2);
  
  // Left ankle acceleration
  gait_setpoints_msg.l_ankle_ddot_ref.linear.x = l_ddot_ankle_position_ref_(0);
  gait_setpoints_msg.l_ankle_ddot_ref.linear.y = l_ddot_ankle_position_ref_(1);
  gait_setpoints_msg.l_ankle_ddot_ref.linear.z = l_ddot_ankle_position_ref_(2);

  // Publish linear velocity
  gait_setpoints_msg.vel_ref.linear.x = gait_velocity_;

  gait_setpoints_pub_->publish(gait_setpoints_msg);
}

void GaitGenerator::timer_callback() {
  updateGaitPercentage();
  gaitSetpointsPub();
}


float GaitGenerator::vectorModule(Eigen::Matrix<float, 3, 1> A) {
  return sqrt(A(0) * A(0) + A(1)*A(1) + A(2)*A(2));
}

void GaitGenerator::updateDirectGaitKinematics() {
  if (DEBUG)
    std::cout << "[GAIT GENERATOR] Get Direct Gait Kinematics" << std::endl;

  // Get XYZ joints positions relatives to each joint
  Eigen::Matrix<float, 3, 1> K_p;
  K_p << -leg_dimentions_;

  for (int i = 0; i < 100; i++) {
    // Get rotation matrixes for each joint
    R_jointRot_pelvis_ = getRotationMatrix(pelvis_jointRot_roll_[i], pelvis_jointRot_pitch_[i], pelvis_jointRot_yaw_[i]);
    R_jointRot_hip_ = getRotationMatrix(hip_jointRot_roll_[i], hip_jointRot_pitch_[i], hip_jointRot_yaw_[i]);
    R_jointRot_knee_ = getRotationMatrix(knee_jointRot_roll_[i], knee_jointRot_pitch_[i], knee_jointRot_yaw_[i]);
    R_jointRot_ankle_ = getRotationMatrix(ankle_jointRot_roll_[i], ankle_jointRot_pitch_[i], ankle_jointRot_yaw_[i]);

    R_pelvis_absolute_ = R_jointRot_pelvis_;
    R_hip_absolute_ = R_pelvis_absolute_ * R_jointRot_hip_;
    R_knee_absolute_ = R_hip_absolute_ * R_jointRot_knee_;
    R_ankle_absolute_ = R_knee_absolute_ * R_jointRot_ankle_;

    // Get direct kinematics of each joint
    hip_relative_xyz_[i]        = R_jointRot_pelvis_ * pelvis_dimentions_;
    knee_relative_xyz_[i]       = R_jointRot_pelvis_ * R_jointRot_hip_ * thigh_dimentions_;
    ankle_relative_xyz_[i]      = R_jointRot_pelvis_ * R_jointRot_hip_ * R_jointRot_knee_ * shank_dimentions_;
    front_foot_relative_xyz_[i] = R_jointRot_pelvis_ * R_jointRot_hip_ * R_jointRot_knee_ * R_jointRot_ankle_ * foot_dimentions_;
    heel_foot_relative_xyz_[i]  = R_jointRot_pelvis_ * R_jointRot_hip_ * R_jointRot_knee_ * R_jointRot_ankle_ * heel_position_;

    // Get absolut position
    hip_absolute_xyz_[i]        = K_p + hip_relative_xyz_[i];
    knee_absolute_xyz_[i]       = hip_absolute_xyz_[i] + knee_relative_xyz_[i];
    ankle_absolute_xyz_[i]      = knee_absolute_xyz_[i] + ankle_relative_xyz_[i];
    front_foot_absolute_xyz_[i] = ankle_absolute_xyz_[i] + front_foot_relative_xyz_[i];
    heel_foot_absolute_xyz_[i]  = ankle_absolute_xyz_[i] + heel_foot_relative_xyz_[i];
  }
}

void GaitGenerator::updateSpatioTemporalGaitData() {
  if (DEBUG)
    std::cout << "[GAIT GENERATOR] Calculate Spatio Temporal Data" << std::endl;

  // Get step length
  // Add al x increments
  float computed_step_length = 0.0F;
  bool valid_ankle_path = true;
  for (int i = 1; i < static_cast<int>(ankle_absolute_xyz_.size()); i++) {
      const float x_now = ankle_absolute_xyz_[i](0);
      const float x_prev = ankle_absolute_xyz_[i-1](0);
      if (!std::isfinite(x_now) || !std::isfinite(x_prev)) {
        valid_ankle_path = false;
        break;
      }
      const float x_increment = std::abs(x_now - x_prev);
      computed_step_length += x_increment;
  }

  if (!valid_ankle_path || !std::isfinite(computed_step_length) || computed_step_length <= 1.0e-6F) {
    RCLCPP_WARN_THROTTLE(
      this->get_logger(), *this->get_clock(), 2000,
      "[GAIT GENERATOR] Invalid ankle path for spatio-temporal data "
      "(valid=%d, computed_step_length=%.6f). Freezing gait timing.",
      valid_ankle_path ? 1 : 0, computed_step_length);
    step_length_ = 0.0F;
    step_period_ = 0.0F;
    t_increment_ = 0.0F;
    gaint_increment_ = 0.0F;
  } else {
    step_length_ = computed_step_length;
  }
  
  // Get gait period. During startup the treadmill and gait generator can be
  // configured a few seconds apart, so gait_velocity_ may legitimately be zero.
  // Keep the gait frozen instead of producing inf/NaN timing values.
  if (step_length_ <= 1.0e-6F || std::abs(gait_velocity_) < 1.0e-6F || gait_n_ <= 0) {
    step_period_ = 0.0F;
    t_increment_ = 0.0F;
    gaint_increment_ = 0.0F;
  } else {
    step_period_ = step_length_ / gait_velocity_;
    t_increment_ = step_period_ / gait_n_;
    const float loop_rate = delta_t_; // In seconds
    gaint_increment_ = 100.00F / (step_period_ / loop_rate);
    if (!std::isfinite(step_period_) || !std::isfinite(t_increment_) || !std::isfinite(gaint_increment_)) {
      RCLCPP_WARN_THROTTLE(
        this->get_logger(), *this->get_clock(), 2000,
        "[GAIT GENERATOR] Non-finite gait timing computed "
        "(length=%.6f, velocity=%.6f, n=%d). Freezing gait timing.",
        step_length_, gait_velocity_, gait_n_);
      step_period_ = 0.0F;
      t_increment_ = 0.0F;
      gaint_increment_ = 0.0F;
    }
  }
  std::cout << "Gait increment: " << gaint_increment_ << std::endl;
  std::cout << "Gait period: " << step_period_ << std::endl;
  std::cout << "Gait length: " << step_length_ << std::endl;
  std::cout << "Gait velocity: " << gait_velocity_ << std::endl;
  std::cout << "Gait n: " << gait_n_ << std::endl;
}

Eigen::Matrix<float, 3, 3> GaitGenerator::getRotationMatrix(float roll,float pitch,float yaw) {
  Eigen::Matrix<float, 3, 3> R_x;
  Eigen::Matrix<float, 3, 3> R_y;
  Eigen::Matrix<float, 3, 3> R_z;
  Eigen::Matrix<float, 3, 3> R_general;

  R_x << 1          ,     0     ,    0     ,
         0          , cos(roll) ,-sin(roll),
         0          , sin(roll) , cos(roll);

  R_y << cos(pitch) ,     0     , sin(pitch),
             0      ,     1     ,     0     ,
        -sin(pitch) ,     0     , cos(pitch);

  R_z << cos(yaw)   , -sin(yaw) ,     0     ,
         sin(yaw)   ,  cos(yaw) ,     0     ,
            0       ,     0     ,     1     ;

  R_general = R_z * R_y * R_x;

  // Generate rotation matrix from quaternion:
   return R_general;
}

void GaitGenerator::getReferenceAnglesData(const std::string& file_name, std::vector<float>& angle_x, std::vector<float>& angle_y, std::vector<float>& angle_z) {
  std::ifstream anglesFile(file_name);
  std::cout << "Read: " << file_name << std::endl;
  if (anglesFile.is_open()) {
      std::string line;
      std::string token;
      auto read_axis = [this, &file_name](std::ifstream& file, std::vector<float>& target, const char* axis_name) {
        std::string axis_line;
        std::string axis_token;
        if (!std::getline(file, axis_line)) {
          RCLCPP_WARN(
            this->get_logger(),
            "[GAIT GENERATOR] Reference angle file %s is missing axis %s.",
            file_name.c_str(), axis_name);
          return;
        }
        std::istringstream stream(axis_line);
        size_t count = 0;
        while (std::getline(stream, axis_token, ',') && count < target.size()) {
          try {
            target[count] = static_cast<float>(std::stof(axis_token) * TORADIANS);
            count++;
          } catch (const std::exception& e) {
            RCLCPP_WARN(
              this->get_logger(),
              "[GAIT GENERATOR] Invalid token in %s axis %s at sample %zu: %s",
              file_name.c_str(), axis_name, count, axis_token.c_str());
          }
        }
        if (count != target.size()) {
          RCLCPP_WARN(
            this->get_logger(),
            "[GAIT GENERATOR] Reference angle file %s axis %s loaded %zu samples; expected %zu.",
            file_name.c_str(), axis_name, count, target.size());
        }
      };
      // Get X (angle_x)
      read_axis(anglesFile, angle_x, "X");
      // Get Y (angle_y)
      read_axis(anglesFile, angle_y, "Y");
      // Get Z (angle_z)
      read_axis(anglesFile, angle_z, "Z");

  } else {
    RCLCPP_ERROR(this->get_logger(), "[GAIT GENERATOR] Could NOT read reference angle file: %s", file_name.c_str());
  }
}

void GaitGenerator::updateGaitGainsXYZ(float step_x_gain, float step_y_gain, float step_z_gain) {
  if (DEBUG)
    std::cout << "[GAIT GENERATOR] Update Gains" << std::endl;
  // Set XYZ gains
  gait_x_gain_ = step_x_gain;
  gait_y_gain_ = step_y_gain;
  gait_z_gain_ = step_z_gain;

  // Corresponding roll pitch yaw:
  pelvis_jointRot_roll_gain_ = gait_y_gain_;
  pelvis_jointRot_pitch_gain_ = gait_x_gain_;
  pelvis_jointRot_yaw_gain_ = gait_y_gain_;


  if (control_amplitudes_ == 0) {
    float thresh_hold = 0.3;
    ankle_jointRot_roll_gain_ = gait_y_gain_;
    ankle_jointRot_pitch_gain_ = gait_x_gain_;
    ankle_jointRot_yaw_gain_ = gait_z_gain_;
    if (gait_x_gain_ > 0.14 && gait_x_gain_ < 0.5) {
      hip_jointRot_roll_gain_ = gait_y_gain_ + thresh_hold;
      hip_jointRot_pitch_gain_ = gait_x_gain_ + thresh_hold;
      hip_jointRot_yaw_gain_ = gait_z_gain_ + thresh_hold;
    } else {
      if(gait_x_gain_ < 0.14) {
        hip_jointRot_roll_gain_ = gait_y_gain_ + thresh_hold;
        hip_jointRot_pitch_gain_ = gait_x_gain_ + thresh_hold;
        hip_jointRot_yaw_gain_ = gait_z_gain_ + thresh_hold;
      } else {
        hip_jointRot_roll_gain_ = 0.5+ thresh_hold;
        hip_jointRot_pitch_gain_ = 0.5+ thresh_hold;
        hip_jointRot_yaw_gain_ =  0.5+ thresh_hold;
      }
    }

    if (gait_x_gain_ > 0.3 && gait_x_gain_ < 0.6) {
      knee_jointRot_roll_gain_ = gait_y_gain_ + thresh_hold;
      knee_jointRot_pitch_gain_ = gait_x_gain_ + thresh_hold;
      knee_jointRot_yaw_gain_ = gait_z_gain_ + thresh_hold;
    } else {
      if (gait_x_gain_ < 0.3) {
        knee_jointRot_roll_gain_ = gait_y_gain_ + thresh_hold;
        knee_jointRot_pitch_gain_ = gait_x_gain_ + thresh_hold;
        knee_jointRot_yaw_gain_ = gait_z_gain_ + thresh_hold;
      } else {
        knee_jointRot_roll_gain_ = 0.6;
        knee_jointRot_pitch_gain_ = 0.6;
        knee_jointRot_yaw_gain_ = 0.6;
      }
    }
  } else {
    hip_jointRot_roll_gain_ = 1;
    hip_jointRot_pitch_gain_ = 1;
    hip_jointRot_yaw_gain_ = 1;

    knee_jointRot_roll_gain_ = 1;
    knee_jointRot_pitch_gain_ = 1;
    knee_jointRot_yaw_gain_ = 1;

    ankle_jointRot_roll_gain_ = 1;
    ankle_jointRot_pitch_gain_ = 1;
    ankle_jointRot_yaw_gain_ = 1;
  }

  // Update gains for loaded reference angles
  for (int i = 0; i < 100; i++)
  {
    pelvis_jointRot_roll_[i] = pelvis_jointRot_roll_file_[i];
    pelvis_jointRot_pitch_[i] = pelvis_jointRot_pitch_file_[i] ;
    pelvis_jointRot_yaw_[i] = pelvis_jointRot_yaw_file_[i];

    hip_jointRot_roll_[i] = hip_jointRot_roll_file_[i];
    hip_jointRot_pitch_[i] = (hip_jointRot_pitch_file_[i] + 4) ;
    hip_jointRot_yaw_[i] = hip_jointRot_yaw_file_[i];

    knee_jointRot_roll_[i] = knee_jointRot_roll_file_[i];
    knee_jointRot_pitch_[i] = (knee_jointRot_pitch_file_[i]);
    knee_jointRot_yaw_[i] = knee_jointRot_yaw_file_[i];

    ankle_jointRot_roll_[i] = ankle_jointRot_roll_file_[i];
    ankle_jointRot_pitch_[i] = ankle_jointRot_pitch_file_[i];
    ankle_jointRot_yaw_[i] = ankle_jointRot_yaw_file_[i];

  }
  applySagittalScaling();


  // Update gains for loaded reference angles
  for (int i = 0; i < 100; i++)
  {    
    pelvis_jointRot_roll_absolute_[i] = pelvis_jointRot_roll_[i];
    pelvis_jointRot_pitch_absolute_[i] = pelvis_jointRot_pitch_[i];
    pelvis_jointRot_yaw_absolute_[i] = pelvis_jointRot_yaw_[i];

    hip_jointRot_roll_absolute_[i] = pelvis_jointRot_roll_absolute_[i] + hip_jointRot_roll_[i];
    hip_jointRot_pitch_absolute_[i] = pelvis_jointRot_pitch_absolute_[i] + hip_jointRot_pitch_[i];
    hip_jointRot_yaw_absolute_[i] = pelvis_jointRot_yaw_absolute_[i] + hip_jointRot_yaw_[i];

    knee_jointRot_roll_absolute_[i] = hip_jointRot_roll_absolute_[i] + knee_jointRot_roll_[i];
    knee_jointRot_pitch_absolute_[i] = hip_jointRot_pitch_absolute_[i] + knee_jointRot_pitch_[i];
    knee_jointRot_yaw_absolute_[i] = hip_jointRot_yaw_absolute_[i] + knee_jointRot_yaw_[i];

    ankle_jointRot_roll_absolute_[i] = knee_jointRot_roll_absolute_[i] + ankle_jointRot_roll_[i];
    ankle_jointRot_pitch_absolute_[i] = knee_jointRot_pitch_absolute_[i] + ankle_jointRot_pitch_[i];
    ankle_jointRot_yaw_absolute_[i] = knee_jointRot_yaw_absolute_[i] + ankle_jointRot_yaw_[i];
  }
}

void GaitGenerator::computeBounds(
  const std::vector<float>& data,
  float gain,
  float& lo,
  float& hi
) {
  
  if (data.empty()) {
    lo = hi = 0.0f;
    return;
  }
  
  // find original extrema
  auto mm = std::minmax_element(data.begin(), data.end());
  float origMin = *mm.first;
  float origMax = *mm.second;
  
  
  // compute midpoint and scaled range
  float mid = 0.5f * (origMin + origMax);
  float range = origMax - origMin;
  float newRange = range * gain;
  

  
  // set new bounds symmetrically about mid
  lo = mid - 0.5f * newRange;
  hi = mid + 0.5f * newRange;
  
  
  // Check for potential issues
  if (!std::isfinite(lo) || !std::isfinite(hi)) {
  }
  
  if (lo >= hi) {
  }
}

void GaitGenerator::applySagittalScaling() {
  std::cout << "[GAIT GENERATOR] Apply Sagittal Scaling" << std::endl;
  const size_t N = 100;
  std::vector<float> t;
  for (int i = 0; i < 100; i++)
   t.push_back(i);

  
  std::cout << "[GAIT GENERATOR] Apply Sagittal Scaling: t size: " << t.size() << std::endl;

  // --- HIP pitch remap ---
  {
    float lo, hi;
    computeBounds(hip_jointRot_pitch_file_, hip_jointRot_pitch_gain_, lo, hi);
    std::cout << "[GAIT GENERATOR] Apply Sagittal Scaling: lo: " << lo << " hi: " << hi << std::endl;
    // these four will be populated by parametrizeAndRemap:
    float t0_min, t0_max, y0_min, y0_max;

    parametrizeAndRemap(
      hip_jointRot_pitch_file_,
      t,
      lo, hi,
      hip_jointRot_pitch_,
      t0_min, t0_max,
      y0_min, y0_max
    );
    std::cout << "[GAIT GENERATOR] Apply Sagittal Scaling: t0_min: " << t0_min << " t0_max: " << t0_max << std::endl;
  }

  // --- KNEE pitch remap ---
  {
    float lo, hi;
    computeBounds(knee_jointRot_pitch_file_, knee_jointRot_pitch_gain_, lo, hi);
    float t0_min, t0_max, y0_min, y0_max;
    parametrizeAndRemap(
      knee_jointRot_pitch_file_,
      t,
      lo, hi,
      knee_jointRot_pitch_,
      t0_min, t0_max,
      y0_min, y0_max
    );
  }

  // --- ANKLE pitch remap ---
  {
    float lo, hi;
    computeBounds(ankle_jointRot_pitch_file_, ankle_jointRot_pitch_gain_, lo, hi);
    float t0_min, t0_max, y0_min, y0_max;
    parametrizeAndRemap(
      ankle_jointRot_pitch_file_,
      t,
      lo, hi,
      ankle_jointRot_pitch_,
      t0_min, t0_max,
      y0_min, y0_max
    );
  }

  // --- PELVIS yaw remap ---
  {
    float lo, hi;
    computeBounds(pelvis_jointRot_yaw_file_, pelvis_jointRot_yaw_gain_, lo, hi);
    float t0_min, t0_max, y0_min, y0_max;
    parametrizeAndRemap(
      pelvis_jointRot_yaw_file_,
      t,
      lo, hi,
      pelvis_jointRot_yaw_,
      t0_min, t0_max,
      y0_min, y0_max
    );
  }

}


// in gait_generator.cpp
void GaitGenerator::applyCartesianScaling() {
  std::cout << "[GAIT GENERATOR] Apply Cartesian Remap (X–Z)" << std::endl;

  const size_t N = ankle_absolute_xyz_.size();
  if (N == 0) return;

  // 1) Build a “phase” vector t = [0,1,2,…,N-1]
  std::vector<float> t(N);
  for (size_t i = 0; i < N; ++i) t[i] = float(i);

  // 2) Extract original X and Z trajectories
  std::vector<float> x_orig(N), z_orig(N);
  // extract X/Y into simple arrays
  for (size_t i=0; i<N; ++i) {
    x_orig[i] = ankle_absolute_xyz_[i](0);  // X coordinate
    z_orig[i] = ankle_absolute_xyz_[i](2);  // Z coordinate (height, using index 2)
  }

  if (!std::isfinite(hip_amplitude_) || !std::isfinite(knee_amplitude_) || !std::isfinite(ankle_amplitude_)) {
    RCLCPP_WARN(
      this->get_logger(),
      "[GAIT GENERATOR] Skipping Cartesian remap because amplitude bounds are non-finite "
      "(hip=%.6f, knee=%.6f, ankle=%.6f).",
      hip_amplitude_, knee_amplitude_, ankle_amplitude_);
    return;
  }

  ankleTraj_gain_x_ = hip_amplitude_;
  ankleTraj_gain_y_ = knee_amplitude_;

  float x_new_min = -hip_amplitude_ /2;
  float x_new_max = hip_amplitude_ /2;
  float z_new_min = knee_amplitude_;
  float z_new_max = ankle_amplitude_;

  // now directly remap into your desired [min,max] ranges:
  std::vector<float> x_new, z_new;
  float t0min, t0max, y0min, y0max;

  // 1) X remap using user-provided bounds
  parametrizeAndRemap(
    /*y=*/       x_orig,
    /*t=*/       t,
    /*new_min=*/ x_new_min,
    /*new_max=*/ x_new_max,
    /*y_out=*/   x_new,
    /*t0_min=*/  t0min,
    /*t0_max=*/  t0max,
    /*y0_min=*/  y0min,
    /*y0_max=*/  y0max
  );

  // 2) Z remap using user-provided bounds
  parametrizeAndRemap(
    /*y=*/       z_orig,
    /*t=*/       t,
    /*new_min=*/ z_new_min,
    /*new_max=*/ z_new_max,
    /*y_out=*/   z_new,
    /*t0_min=*/  t0min,
    /*t0_max=*/  t0max,
    /*y0_min=*/  y0min,
    /*y0_max=*/  y0max
  );

  if (x_new.size() != N || z_new.size() != N) {
    RCLCPP_WARN(
      this->get_logger(),
      "[GAIT GENERATOR] Cartesian remap returned invalid sizes x=%zu z=%zu expected=%zu; keeping previous trajectory.",
      x_new.size(), z_new.size(), N);
    return;
  }

  // write back into your ankle_absolute_xyz_
  for (size_t i=0; i<N; ++i) {
    ankle_absolute_xyz_[i](0) = x_new[i];
    ankle_absolute_xyz_[i](2) = z_new[i] + 1 *foot_dimentions_(2);
  }
}

void GaitGenerator::parametrizeAndRemap(
  const std::vector<float>& y,
  const std::vector<float>& t,
  float new_min,
  float new_max,
  std::vector<float>& y_new,
  float& t0_min,
  float& t0_max,
  float& y0_min,
  float& y0_max
)  {
  size_t N = y.size();
  y_new.resize(N);
  if (N == 0 || t.size() != N) {
    y_new.clear();
    RCLCPP_WARN(
      this->get_logger(),
      "[GAIT GENERATOR] Cannot remap trajectory with invalid sizes y=%zu t=%zu.",
      y.size(), t.size());
    return;
  }

  if (!std::isfinite(new_min) || !std::isfinite(new_max)) {
    RCLCPP_WARN(
      this->get_logger(),
      "[GAIT GENERATOR] Non-finite remap bounds new_min=%.6f new_max=%.6f; copying finite source trajectory.",
      new_min, new_max);
    for (size_t i = 0; i < N; ++i) {
      y_new[i] = std::isfinite(y[i]) ? y[i] : 0.0F;
    }
    t0_min = t0_max = y0_min = y0_max = 0.0F;
    return;
  }

  if (new_min > new_max) {
    std::swap(new_min, new_max);
  }

  // 1) Find original extrema indices, ignoring non-finite source samples.
  bool have_finite_sample = false;
  size_t i_min = 0;
  size_t i_max = 0;
  y0_min = 0.0F;
  y0_max = 0.0F;
  for (size_t i = 0; i < N; ++i) {
    if (!std::isfinite(y[i])) {
      continue;
    }
    if (!have_finite_sample) {
      have_finite_sample = true;
      i_min = i_max = i;
      y0_min = y0_max = y[i];
      continue;
    }
    if (y[i] < y0_min) {
      y0_min = y[i];
      i_min = i;
    }
    if (y[i] > y0_max) {
      y0_max = y[i];
      i_max = i;
    }
  }
  if (!have_finite_sample) {
    const float target_mid = 0.5F * (new_min + new_max);
    RCLCPP_WARN(
      this->get_logger(),
      "[GAIT GENERATOR] Cannot remap trajectory with no finite samples; filling with %.6f.",
      target_mid);
    std::fill(y_new.begin(), y_new.end(), target_mid);
    t0_min = t0_max = y0_min = y0_max = 0.0F;
    return;
  }

  // 2) Extract phase locations
  t0_min = t[i_min];
  t0_max = t[i_max];

  // 3) Build shape function f = (y - y0_min)/(y0_max - y0_min)
  float denom = (y0_max - y0_min);
  if (!std::isfinite(y0_min) || !std::isfinite(y0_max) || !std::isfinite(denom) || std::abs(denom) < 1.0e-6F) {
    const float target_mid = 0.5F * (new_min + new_max);
    RCLCPP_WARN_THROTTLE(
      this->get_logger(), *this->get_clock(), 2000,
      "[GAIT GENERATOR] Degenerate trajectory remap "
      "(y0_min=%.6f, y0_max=%.6f, new_min=%.6f, new_max=%.6f). Filling with %.6f.",
      y0_min, y0_max, new_min, new_max, target_mid);
    std::fill(y_new.begin(), y_new.end(), target_mid);
    return;
  }
  for (size_t i = 0; i < N; ++i) {
    const float source_value = std::isfinite(y[i]) ? y[i] : y0_min;
    float f = (source_value - y0_min) / denom;
    if (!std::isfinite(f)) {
      f = 0.0F;
    }
    // 4) Affine remap
    y_new[i] = new_min + f * (new_max - new_min);
  }
}

void GaitGenerator::loadAnglesFromFiles() {
  if (DEBUG)
    std::cout << "[GAIT GENERATOR] Load Angles From Files" << std::endl;

  std::string package_share_directory = ament_index_cpp::get_package_share_directory("gait_generator");
  const std::string reference_dir = package_share_directory + "/reference_angles/";
  getReferenceAnglesData(reference_dir + "P_XYZ.txt", pelvis_jointRot_roll_file_, pelvis_jointRot_pitch_file_, pelvis_jointRot_yaw_file_);
  getReferenceAnglesData(reference_dir + "H_XYZ.txt", hip_jointRot_roll_file_, hip_jointRot_pitch_file_, hip_jointRot_yaw_file_);
  getReferenceAnglesData(reference_dir + "K_XYZ.txt", knee_jointRot_roll_file_, knee_jointRot_pitch_file_, knee_jointRot_yaw_file_);
  getReferenceAnglesData(reference_dir + "A_XYZ.txt", ankle_jointRot_roll_file_, ankle_jointRot_pitch_file_, ankle_jointRot_yaw_file_);

  // Find toe off as the minimun value in ankle angle trajectory
  toe_off_index_ = 50;

}

void GaitGenerator::updateGaitParameters(float patient_height, float gait_velocity) {
  // Update gait atributes
  updateBodySegments(patient_height);
  updateGaitGainsXYZ(gait_velocity/1, gait_velocity/1, gait_velocity/1);
  gait_velocity_ = gait_velocity;

  // Update gait kinematics
  updateDirectGaitKinematics();

  // Update paths
  updateAnklePelvisPathsXYZ();
}

void GaitGenerator::updateGaitPercentage() {
  double t_scale_right = 1;
  double t_scale_left = 1;
  auto wrap_phase = [](double phase) {
    phase = fmod(phase, 100);
    return phase < 0 ? phase + 100 : phase;
  };
  if (control_mode_ == 2) {
    // Use adaptive oscillator's gait cycle percentage
    gait_percentage_ = wrap_phase(adaptive_gait_cycle_right_);
    r_gaint_index_ = wrap_phase(adaptive_gait_cycle_right_);
    l_gaint_index_ = wrap_phase(adaptive_gait_cycle_left_);

    double t_adapted_increment_left = (1/adaptive_estimated_freq_left_) / gait_n_;
    double t_adapted_increment_right = (1/adaptive_estimated_freq_right_) / gait_n_;
    if (std::isfinite(t_increment_) && std::isfinite(t_adapted_increment_right) && t_adapted_increment_right != 0) {
      t_scale_right = t_increment_ / t_adapted_increment_right;
    } else {
      t_scale_right = 1;
    }

    if (std::isfinite(t_increment_) && std::isfinite(t_adapted_increment_left) && t_adapted_increment_left != 0) {
      t_scale_left = t_increment_ / t_adapted_increment_left;
    } else {
      t_scale_left = 1;
    }
  } else {
    // Default behavior
    gait_percentage_ = wrap_phase(gait_percentage_ + gaint_increment_); // 0-100%
    r_gaint_index_ = gait_percentage_;
    l_gaint_index_ = wrap_phase(gait_percentage_ + 50);
  }

  // Compute x offset
  double x_offset = ( *std::max_element(ankle_path_xyz_[0].begin(), ankle_path_xyz_[0].end()) + 
                      *std::min_element(ankle_path_xyz_[0].begin(), ankle_path_xyz_[0].end()) ) / 2;
  // Update data form splines
  // Pelvis
  pelvis_position_ref_(0) = interp_pelvis_x_(gait_percentage_) - 0.03 * patient_height_;
  pelvis_position_ref_(1) = interp_pelvis_y_(gait_percentage_);
  pelvis_position_ref_(2) = interp_pelvis_z_(gait_percentage_);
  pelvis_orientation_ref_(0) = interp_pelvis_jointRot_roll_(gait_percentage_);
  pelvis_orientation_ref_(1) = interp_pelvis_jointRot_pitch_(gait_percentage_);
  pelvis_orientation_ref_(2) = interp_pelvis_jointRot_yaw_(gait_percentage_);

  // Hip
  // Position
  r_hip_position_ref_(0) = interp_hip_x_(r_gaint_index_);
  r_hip_position_ref_(1) = interp_hip_y_(r_gaint_index_);
  r_hip_position_ref_(2) = interp_hip_z_(r_gaint_index_);
  l_hip_position_ref_(0) = interp_hip_x_(l_gaint_index_);
  l_hip_position_ref_(1) = -interp_hip_y_(l_gaint_index_);
  l_hip_position_ref_(2) = interp_hip_z_(l_gaint_index_);

  // Orientation
  r_hip_orientation_ref_(0) = interp_hip_jointRot_roll_(r_gaint_index_);
  r_hip_orientation_ref_(1) = interp_hip_jointRot_pitch_(r_gaint_index_);
  r_hip_orientation_ref_(2) = interp_hip_jointRot_yaw_(r_gaint_index_);
  l_hip_orientation_ref_(0) = interp_hip_jointRot_roll_(l_gaint_index_);
  l_hip_orientation_ref_(1) = interp_hip_jointRot_pitch_(l_gaint_index_);
  l_hip_orientation_ref_(2) = interp_hip_jointRot_yaw_(l_gaint_index_);

  // Velocity
  r_dot_hip_position_ref_(0) = interp_hip_dot_x_(r_gaint_index_) * t_scale_right;
  r_dot_hip_position_ref_(1) = interp_hip_dot_y_(r_gaint_index_) * t_scale_right;
  r_dot_hip_position_ref_(2) = interp_hip_dot_z_(r_gaint_index_) * t_scale_right;
  l_dot_hip_position_ref_(0) = interp_hip_dot_x_(l_gaint_index_) * t_scale_right;
  l_dot_hip_position_ref_(1) = -interp_hip_dot_y_(l_gaint_index_) * t_scale_right;
  l_dot_hip_position_ref_(2) = interp_hip_dot_z_(l_gaint_index_) * t_scale_right;

  // Knee
  // Position
  r_knee_position_ref_(0) = interp_knee_x_(r_gaint_index_);
  r_knee_position_ref_(1) = interp_knee_y_(r_gaint_index_);
  r_knee_position_ref_(2) = interp_knee_z_(r_gaint_index_);
  l_knee_position_ref_(0) = interp_knee_x_(l_gaint_index_);
  l_knee_position_ref_(1) = -interp_knee_y_(l_gaint_index_);
  l_knee_position_ref_(2) = interp_knee_z_(l_gaint_index_);
  // Orientation
  r_knee_orientation_ref_(0) = interp_knee_jointRot_roll_(r_gaint_index_);
  r_knee_orientation_ref_(1) = interp_knee_jointRot_pitch_(r_gaint_index_);
  r_knee_orientation_ref_(2) = interp_knee_jointRot_yaw_(r_gaint_index_);
  l_knee_orientation_ref_(0) = interp_knee_jointRot_roll_(l_gaint_index_);
  l_knee_orientation_ref_(1) = interp_knee_jointRot_pitch_(l_gaint_index_);
  l_knee_orientation_ref_(2) = interp_knee_jointRot_yaw_(l_gaint_index_);
  // Velocity
  r_dot_knee_position_ref_(0) = interp_knee_dot_x_(r_gaint_index_) * t_scale_right;
  r_dot_knee_position_ref_(1) = interp_knee_dot_y_(r_gaint_index_) * t_scale_right;
  r_dot_knee_position_ref_(2) = interp_knee_dot_z_(r_gaint_index_) * t_scale_right;
  l_dot_knee_position_ref_(0) = interp_knee_dot_x_(l_gaint_index_) * t_scale_right;
  l_dot_knee_position_ref_(1) = interp_knee_dot_y_(l_gaint_index_) * t_scale_right;
  l_dot_knee_position_ref_(2) = interp_knee_dot_z_(l_gaint_index_) * t_scale_right;

  // Ankle
  // Position
  r_ankle_position_ref_(0) = interp_ankle_x_(r_gaint_index_) - x_offset;
  r_ankle_position_ref_(1) = interp_ankle_y_(r_gaint_index_);
  r_ankle_position_ref_(2) = interp_ankle_z_(r_gaint_index_) - 1 *foot_dimentions_(2);
  l_ankle_position_ref_(0) = interp_ankle_x_(l_gaint_index_) - x_offset;
  l_ankle_position_ref_(1) = -interp_ankle_y_(l_gaint_index_);
  l_ankle_position_ref_(2) = interp_ankle_z_(l_gaint_index_) - 1 *foot_dimentions_(2);

  // Orientation
  r_ankle_orientation_ref_(0) = interp_ankle_jointRot_roll_(r_gaint_index_);
  r_ankle_orientation_ref_(1) = interp_ankle_jointRot_pitch_(r_gaint_index_);
  r_ankle_orientation_ref_(2) = interp_ankle_jointRot_yaw_(r_gaint_index_);
  l_ankle_orientation_ref_(0) = interp_ankle_jointRot_roll_(l_gaint_index_);
  l_ankle_orientation_ref_(1) = interp_ankle_jointRot_pitch_(l_gaint_index_);
  l_ankle_orientation_ref_(2) = interp_ankle_jointRot_yaw_(l_gaint_index_);
  // Velocity
  r_dot_ankle_position_ref_(0) = interp_ankle_dot_x_(r_gaint_index_) * t_scale_right;
  r_dot_ankle_position_ref_(1) = interp_ankle_dot_y_(r_gaint_index_) * t_scale_right;
  r_dot_ankle_position_ref_(2) = interp_ankle_dot_z_(r_gaint_index_) * t_scale_right;
  l_dot_ankle_position_ref_(0) = interp_ankle_dot_x_(l_gaint_index_) * t_scale_left;
  l_dot_ankle_position_ref_(1) = interp_ankle_dot_y_(l_gaint_index_) * t_scale_left;
  l_dot_ankle_position_ref_(2) = interp_ankle_dot_z_(l_gaint_index_) * t_scale_left;

  // Update accelerations
  r_ddot_ankle_position_ref_(0) = interp_ankle_ddot_x_(r_gaint_index_) * t_scale_right * t_scale_right;
  r_ddot_ankle_position_ref_(1) = interp_ankle_ddot_y_(r_gaint_index_) * t_scale_right * t_scale_right;
  r_ddot_ankle_position_ref_(2) = interp_ankle_ddot_z_(r_gaint_index_) * t_scale_right * t_scale_right;
  l_ddot_ankle_position_ref_(0) = interp_ankle_ddot_x_(l_gaint_index_) * t_scale_left * t_scale_left;
  l_ddot_ankle_position_ref_(1) = interp_ankle_ddot_y_(l_gaint_index_) * t_scale_left * t_scale_left;
  l_ddot_ankle_position_ref_(2) = interp_ankle_ddot_z_(l_gaint_index_) * t_scale_left * t_scale_left;
}

void GaitGenerator::updateBodySegments(float patient_height) {
  patient_height_ = patient_height;
  // Update body segments dimentions:
  leg_dimentions_ = leg_proportion_* patient_height_;
  pelvis_dimentions_ = pelvis_proportion_ * patient_height_ ;
  shank_dimentions_ = shank_proportion_ * patient_height_;
  thigh_dimentions_ = thigh_proportion_ * patient_height_;
  foot_dimentions_ = foot_proportion_ * patient_height_;
  heel_position_ = heel_proportion_ * patient_height_;
}

void GaitGenerator::updateAnklePelvisPathsXYZ() {
  const int N = 100;
  const int shift = N/2;

  std::vector<Eigen::Vector3f> pelvis_trans(N, Eigen::Vector3f::Zero());
  // copy heel-in-pelvis-frame
  std::vector<Eigen::Vector3f> heel1 = ankle_absolute_xyz_;
  std::vector<Eigen::Vector3f> heel2(N);
  for(int i=0;i<N;++i) heel2[i] = heel1[(i+shift)%N];

  // Phase-1: lock left-heel Y,Z
  for(int i=0;i<N;++i){
    pelvis_trans[i].y() = -heel1[i].y();
    pelvis_trans[i].z() = -heel1[i].z();
  }

  // center Y about zero  
  float y_min=pelvis_trans[0].y(), y_max=pelvis_trans[0].y();
  for(auto &p: pelvis_trans){
    y_min = std::min(y_min, p.y());
    y_max = std::max(y_max, p.y());
  }
  float y_mid = 0.5f*(y_min+y_max);
  for(auto &p: pelvis_trans) p.y() -= y_mid;


  // symmetrize second half
  for(int i=0;i<N/2;++i){
    float y0 = pelvis_trans[i].y();
    float z0 = pelvis_trans[i].z();
    float offset = pelvis_trans[0].y() + pelvis_trans[N/2-1].y();
    pelvis_trans[i+shift].y() = -y0 + offset;
    pelvis_trans[i+shift].z() =  z0;
  }

  // center Y about zero
  y_min=pelvis_trans[0].y(), y_max=pelvis_trans[0].y();
  for(auto &p: pelvis_trans){
    y_min = std::min(y_min, p.y());
    y_max = std::max(y_max, p.y());
  }
  y_mid = 0.5f*(y_min+y_max);
  for(auto &p: pelvis_trans) p.y() -= y_mid;
  
  smoothGaussian(pelvis_trans,        30);

  // 3) Now shift ALL joints into the WORLD by adding pelvis_trans
  for(int i=0;i<N;++i){
    K_p_absolute_xyz_[i]        =  pelvis_trans[i];
    hip_absolute_xyz_[i]        += pelvis_trans[i];
    knee_absolute_xyz_[i]       += pelvis_trans[i];
    ankle_absolute_xyz_[i]      += pelvis_trans[i];
    front_foot_absolute_xyz_[i] += pelvis_trans[i];
    heel_foot_absolute_xyz_[i]  += pelvis_trans[i];
  }

  // Correct y and x:
  for (int i = 0; i < 100; i++) {
    K_p_absolute_xyz_[i](0)   = -K_p_absolute_xyz_[i](0);
    hip_absolute_xyz_[i](0)   = -hip_absolute_xyz_[i](0);
    knee_absolute_xyz_[i](0)  = -knee_absolute_xyz_[i](0);
    ankle_absolute_xyz_[i](0) = -ankle_absolute_xyz_[i](0);

    K_p_absolute_xyz_[i](1)   = -K_p_absolute_xyz_[i](1)*0.6;
    hip_absolute_xyz_[i](1)   = -hip_absolute_xyz_[i](1);
    knee_absolute_xyz_[i](1)  = -knee_absolute_xyz_[i](1);
    ankle_absolute_xyz_[i](1) = -ankle_absolute_xyz_[i](1);
  }


  // 4) Final Gaussian smoothing & publish
  smoothGaussian(K_p_absolute_xyz_,        20);
  smoothGaussian(hip_absolute_xyz_,        20);
  smoothGaussian(knee_absolute_xyz_,       20);
  smoothGaussian(ankle_absolute_xyz_,      20);
  smoothGaussian(front_foot_absolute_xyz_, 20);
  smoothGaussian(heel_foot_absolute_xyz_,  20);
  
  if (control_amplitudes_ > 0) {
    applyCartesianScaling();
  }  

  // Update temporal gait data
  updateSpatioTemporalGaitData();

  updatePaths();
  
}

// void GaitGenerator::updateAnklePelvisPathsXYZ() {
//   if (DEBUG)
//     std::cout << "\n[GAIT GENERATOR] Generate Ankle and Pelvis Paths" << std::endl;

//   int i_r;
//   int i_l;

//   int r_heel_strike = 1;
//   int l_heel_strike = 1;

//   int r_heel_offset = 0;
//   int l_heel_offset = 0;

//   // Initilize swing
//   int r_swing = 0;
//   int l_swing = 1;

//   // Initialize stance
//   int r_stance = 1;
//   int l_stance = 0;

//   Eigen::Matrix<float, 3, 1> K_O;   // World reference
//   Eigen::Matrix<float, 3, 1> K_p;   // Pelvis reference (position respect world)
//   Eigen::Matrix<float, 3, 1> r_hip_pos_xyz;
//   Eigen::Matrix<float, 3, 1> r_knee_pos_xyz;
//   Eigen::Matrix<float, 3, 1> r_ankle_pos_xyz;
//   Eigen::Matrix<float, 3, 1> r_heel_foot_pos_xyz;
//   Eigen::Matrix<float, 3, 1> r_front_foot_pos_xyz;

//   Eigen::Matrix<float, 3, 1> l_hip_pos_xyz;
//   Eigen::Matrix<float, 3, 1> l_knee_pos_xyz;
//   Eigen::Matrix<float, 3, 1> l_ankle_pos_xyz;
//   Eigen::Matrix<float, 3, 1> l_heel_foot_pos_xyz;
//   Eigen::Matrix<float, 3, 1> l_front_foot_pos_xyz;

//   // Pelvis reference frame
//   K_O << 0,0,0;
//   K_p = -leg_dimentions_;

//   for (int i = 0; i < 100; i++) {
//       i_r = i;
//       i_l = (i_r + 50) % 100;

//       // Calculate joints positions
//       r_hip_pos_xyz           = K_p           + hip_relative_xyz_[i_r];
//       r_knee_pos_xyz          = r_hip_pos_xyz   + knee_relative_xyz_[i_r];
//       r_ankle_pos_xyz         = r_knee_pos_xyz  + ankle_relative_xyz_[i_r];
//       r_heel_foot_pos_xyz     = r_ankle_pos_xyz + heel_foot_relative_xyz_[i_r];
//       r_front_foot_pos_xyz    = r_ankle_pos_xyz + front_foot_relative_xyz_[i_r];


//       l_hip_pos_xyz         = K_p             + hip_relative_xyz_[i_l];
//       l_knee_pos_xyz        = l_hip_pos_xyz   + knee_relative_xyz_[i_l];
//       l_ankle_pos_xyz       = l_knee_pos_xyz  + ankle_relative_xyz_[i_l];
//       l_heel_foot_pos_xyz   = l_ankle_pos_xyz + heel_foot_relative_xyz_[i_l];
//       l_front_foot_pos_xyz  = l_ankle_pos_xyz + front_foot_relative_xyz_[i_l];

//       // DETECT GAIT PHASES AND UPDATE K_P
//       // Stance
//       if (i_l > 0 && i_l < toe_off_index_)
//           l_stance = 1;
//       else
//           l_stance = 0;

//       if (i_r > 0 && i_r < toe_off_index_)
//           r_stance = 1;
//       else
//           r_stance = 0;

//       // Heel Strice
//       if (r_stance == 1 && r_swing == 1)
//           r_heel_strike = 1;
//       else
//           r_heel_strike = 0;

//       if (l_stance == 1 && l_swing == 1)
//           l_heel_strike = 1;
//       else
//           l_heel_strike = 0;

//       // Swing Phase
//       if (i_r >= toe_off_index_)
//           r_swing = 1;
//       else
//           r_swing = 0;

//       if (i_l >= toe_off_index_)
//           l_swing = 1;
//       else
//           l_swing = 0;

//       if (i_r > 1 && i_r < 50)
//           K_p(1) = -heel_foot_absolute_xyz_[i_r](1) + heel_foot_absolute_xyz_[6](1);

//       if (i_l >= 1 && i_l < 50)
//           K_p(1) = heel_foot_absolute_xyz_[i_l](1) - heel_foot_absolute_xyz_[6](1);


//       if (l_stance == 1 && r_stance == 1) {

//           if (i_r < i_l) {
//             if (i_r < 15)
//               K_p(2) = K_p(2) - r_heel_foot_pos_xyz(2);
//             else
//               K_p(2) = K_p(2) - r_front_foot_pos_xyz(2);
//           }

//           if (i_l < i_r) {
//             if (i_l < 15)
//               K_p(2) = K_p(2) - l_heel_foot_pos_xyz(2);
//             else
//               K_p(2) = K_p(2) - l_front_foot_pos_xyz(2);
//           }
//       } else {
//           if (l_stance)
//               if (i_l < 15)
//                   K_p(2) = K_p(2) - l_heel_foot_pos_xyz(2);
//               else
//                   K_p(2) = K_p(2) - l_front_foot_pos_xyz(2);

//           if (r_stance)
//               if (i_r < 15)
//                   K_p(2) = K_p(2) - r_heel_foot_pos_xyz(2);
//               else
//                   K_p(2) = K_p(2) - r_front_foot_pos_xyz(2);
//       }


//       // Save positions
//       K_p_absolute_xyz_[i]   = K_p;
//       hip_absolute_xyz_[i]   = r_hip_pos_xyz;
//       knee_absolute_xyz_[i]  = r_knee_pos_xyz;
//       ankle_absolute_xyz_[i] = r_ankle_pos_xyz;
//   }

//   // Correct y and x:
//   for (int i = 0; i < 100; i++) {
//     K_p_absolute_xyz_[i](0)   = -K_p_absolute_xyz_[i](0);
//     hip_absolute_xyz_[i](0)   = -hip_absolute_xyz_[i](0);
//     knee_absolute_xyz_[i](0)  = -knee_absolute_xyz_[i](0);
//     ankle_absolute_xyz_[i](0) = -ankle_absolute_xyz_[i](0);

//     K_p_absolute_xyz_[i](1)   = -K_p_absolute_xyz_[i](1);
//     hip_absolute_xyz_[i](1)   = -hip_absolute_xyz_[i](1);
//     knee_absolute_xyz_[i](1)  = -knee_absolute_xyz_[i](1);
//     ankle_absolute_xyz_[i](1) = -ankle_absolute_xyz_[i](1);
//   }

//   smoothGaussian(K_p_absolute_xyz_, 5);
//   smoothGaussian(hip_absolute_xyz_, 5);
//   smoothGaussian(knee_absolute_xyz_, 5);
//   smoothGaussian(ankle_absolute_xyz_, 5);

//   // Generate standar c arrays with Paths
//   updatePaths();
// }

// Smooths the vector "data" using a Gaussian filter of size window_size (odd).
// Wrap-around (circular) boundary conditions ensure the start/end connect smoothly.
void GaitGenerator::smoothGaussian(
  std::vector<Eigen::Matrix<float,3,1>>& data,
  const int window_size
) {
  const int N = data.size();
  if (N == 0) return;

  // Enforce an odd window size
  int ws = window_size % 2 == 1 ? window_size : window_size + 1;
  int half = ws/2;

  // Build Gaussian kernel
  std::vector<float> kernel(ws);
  float sigma = ws/3.0f;
  float sum = 0;
  for(int i = 0; i < ws; ++i) {
    int x = i - half;
    kernel[i] = std::exp(-0.5f * x*x / (sigma*sigma));
    sum += kernel[i];
  }
  for(auto &k : kernel) k /= sum;

  // Convolve with wrap-around indexing
  std::vector<Eigen::Matrix<float,3,1>> out(N, Eigen::Matrix<float,3,1>::Zero());
  for(int i = 0; i < N; ++i) {
    Eigen::Matrix<float,3,1> acc = Eigen::Matrix<float,3,1>::Zero();
    for(int j = -half; j <= half; ++j) {
      int idx = i + j;
      // wrap
      idx = (idx % N + N) % N;
      acc += kernel[j + half] * data[idx];
    }
    out[i] = acc;
  }

  data.swap(out);
}



void GaitGenerator::updatePaths() {
  int l_sign = 1;
  // Generate vector paths
  for (int i = 0; i < 100; i++){
    for (int j = 0; j < 3; j++) {
      if (i == 1)
        l_sign = -1;
      else
        l_sign = 1;
      pelvis_path_xyz_[j][i]  = K_p_absolute_xyz_[i](j) - leg_dimentions_(j)*1.1;
      hip_path_xyz_[j][i]   = hip_absolute_xyz_[i](j);
      knee_path_xyz_[j][i]  = knee_absolute_xyz_[i](j);

      // Apply foot dimension correction only for z-axis
      if (j == 2) { // j=2 is the z-axis
        hip_path_xyz_[j][i]  -= 0.5 * foot_dimentions_(j);
        knee_path_xyz_[j][i] -= 0.5 * foot_dimentions_(j);
      }
      
      ankle_path_xyz_[j][i] = ankle_absolute_xyz_[i](j);
    }
  }

  const bool timing_valid = std::isfinite(t_increment_) && t_increment_ > 1.0e-6F;
  if (!timing_valid) {
    RCLCPP_WARN_THROTTLE(
      this->get_logger(), *this->get_clock(), 2000,
      "[GAIT GENERATOR] Invalid gait timing for derivatives "
      "(t_increment=%.6f, gait_velocity=%.3f). Publishing zero ankle velocities/accelerations.",
      t_increment_, gait_velocity_);
    for (int i = 0; i < 100; i++) {
      for (int j = 0; j < 3; j++) {
        hip_path_xyz_dot_[j][i] = 0.0F;
        knee_path_xyz_dot_[j][i] = 0.0F;
        ankle_path_xyz_dot_[j][i] = 0.0F;
        ankle_path_xyz_ddot_[j][i] = 0.0F;
      }
    }
  } else {
    // Get velocities from the paths. Use circular central differences so
    // the 99->0 phase wrap does not inject a spike at the start of swing.
    for (int i = 0; i < 100; i++) {
      const int im = (i + 99) % 100;
      const int ip = (i + 1) % 100;
      for (int j = 0; j < 3; j++) {
        hip_path_xyz_dot_[j][i] =
          (hip_path_xyz_[j][ip] - hip_path_xyz_[j][im]) / (2 * t_increment_);
        knee_path_xyz_dot_[j][i] =
          (knee_path_xyz_[j][ip] - knee_path_xyz_[j][im]) / (2 * t_increment_);
        ankle_path_xyz_dot_[j][i] =
          (ankle_path_xyz_[j][ip] - ankle_path_xyz_[j][im]) / (2 * t_increment_);
      }
    }

    // Calculate accelerations for ankles with the same circular boundary.
    for (int i = 0; i < 100; i++) {
      const int im = (i + 99) % 100;
      const int ip = (i + 1) % 100;
      for (int j = 0; j < 3; j++) {
        ankle_path_xyz_ddot_[j][i] =
          (ankle_path_xyz_dot_[j][ip] - ankle_path_xyz_dot_[j][im]) / (2 * t_increment_);
      }
    }
  }

  // Smooth the velocity and acceleration data with circular boundaries.
  std::vector<Eigen::Matrix<float, 3, 1>> ankle_dot_converted(
    ankle_path_xyz_dot_[0].size(), Eigen::Matrix<float, 3, 1>::Zero());
  for (size_t i = 0; i < ankle_dot_converted.size(); ++i) {
    ankle_dot_converted[i] << ankle_path_xyz_dot_[0][i],
                              ankle_path_xyz_dot_[1][i],
                              ankle_path_xyz_dot_[2][i];
  }
  smoothGaussian(ankle_dot_converted, 5);
  for (size_t i = 0; i < ankle_dot_converted.size(); ++i) {
    ankle_path_xyz_dot_[0][i] = ankle_dot_converted[i](0);
    ankle_path_xyz_dot_[1][i] = ankle_dot_converted[i](1);
    ankle_path_xyz_dot_[2][i] = ankle_dot_converted[i](2);
  }
  // Convert ankle_path_xyz_ddot_ to std::vector<Eigen::Matrix<float, 3, 1>> for smoothing
  std::vector<Eigen::Matrix<float, 3, 1>> ankle_ddot_converted(ankle_path_xyz_ddot_[0].size(), Eigen::Matrix<float, 3, 1>::Zero());
  for (size_t i = 0; i < ankle_path_xyz_ddot_[0].size(); ++i) {
    ankle_ddot_converted[i] << ankle_path_xyz_ddot_[0][i], ankle_path_xyz_ddot_[1][i], ankle_path_xyz_ddot_[2][i];
  }

  smoothGaussian(ankle_ddot_converted, 10);

  // Convert back to std::vector<std::vector<float>> after smoothing
  for (size_t i = 0; i < ankle_ddot_converted.size(); ++i) {
    ankle_path_xyz_ddot_[0][i] = ankle_ddot_converted[i](0);
    ankle_path_xyz_ddot_[1][i] = ankle_ddot_converted[i](1);
    ankle_path_xyz_ddot_[2][i] = ankle_ddot_converted[i](2);
  }
  

  // Generate the interpolated objects
  
  // Pelvis
  interp_pelvis_x_.setData(gait_percentage_vector_,pelvis_path_xyz_[0]);
  interp_pelvis_y_.setData(gait_percentage_vector_,pelvis_path_xyz_[1]);
  interp_pelvis_z_.setData(gait_percentage_vector_,pelvis_path_xyz_[2]);
  interp_pelvis_jointRot_roll_.setData(gait_percentage_vector_,pelvis_jointRot_roll_absolute_);
  interp_pelvis_jointRot_pitch_.setData(gait_percentage_vector_,pelvis_jointRot_pitch_absolute_);
  interp_pelvis_jointRot_yaw_.setData(gait_percentage_vector_,pelvis_jointRot_yaw_absolute_);

  // Hip
  interp_hip_x_.setData(gait_percentage_vector_,hip_path_xyz_[0]);
  interp_hip_y_.setData(gait_percentage_vector_,hip_path_xyz_[1]);
  interp_hip_z_.setData(gait_percentage_vector_,hip_path_xyz_[2]);
  interp_hip_jointRot_roll_.setData(gait_percentage_vector_,hip_jointRot_roll_absolute_);  
  interp_hip_jointRot_pitch_.setData(gait_percentage_vector_,hip_jointRot_pitch_absolute_);
  interp_hip_jointRot_yaw_.setData(gait_percentage_vector_,hip_jointRot_yaw_absolute_);
  interp_hip_dot_x_.setData(gait_percentage_vector_,hip_path_xyz_dot_[0]);
  interp_hip_dot_y_.setData(gait_percentage_vector_,hip_path_xyz_dot_[1]);
  interp_hip_dot_z_.setData(gait_percentage_vector_,hip_path_xyz_dot_[2]);

  // Knee
  interp_knee_x_.setData(gait_percentage_vector_,knee_path_xyz_[0]);
  interp_knee_y_.setData(gait_percentage_vector_,knee_path_xyz_[1]);
  interp_knee_z_.setData(gait_percentage_vector_,knee_path_xyz_[2]);
  interp_knee_jointRot_roll_.setData(gait_percentage_vector_,knee_jointRot_roll_absolute_);
  interp_knee_jointRot_pitch_.setData(gait_percentage_vector_,knee_jointRot_pitch_absolute_);
  interp_knee_jointRot_yaw_.setData(gait_percentage_vector_,knee_jointRot_yaw_absolute_);
  interp_knee_dot_x_.setData(gait_percentage_vector_,knee_path_xyz_dot_[0]);
  interp_knee_dot_y_.setData(gait_percentage_vector_,knee_path_xyz_dot_[1]);
  interp_knee_dot_z_.setData(gait_percentage_vector_,knee_path_xyz_dot_[2]);

  // Ankle
  interp_ankle_x_.setData(gait_percentage_vector_,ankle_path_xyz_[0]);
  interp_ankle_y_.setData(gait_percentage_vector_,ankle_path_xyz_[1]);
  interp_ankle_z_.setData(gait_percentage_vector_,ankle_path_xyz_[2]);
  interp_ankle_jointRot_roll_.setData(gait_percentage_vector_,ankle_jointRot_roll_absolute_);
  interp_ankle_jointRot_pitch_.setData(gait_percentage_vector_,ankle_jointRot_pitch_absolute_);
  interp_ankle_jointRot_yaw_.setData(gait_percentage_vector_,ankle_jointRot_yaw_absolute_);
  interp_ankle_dot_x_.setData(gait_percentage_vector_,ankle_path_xyz_dot_[0]);
  interp_ankle_dot_y_.setData(gait_percentage_vector_,ankle_path_xyz_dot_[1]);
  interp_ankle_dot_z_.setData(gait_percentage_vector_,ankle_path_xyz_dot_[2]);

  // Ankle accelerations
  interp_ankle_ddot_x_.setData(gait_percentage_vector_, ankle_path_xyz_ddot_[0]);
  interp_ankle_ddot_y_.setData(gait_percentage_vector_, ankle_path_xyz_ddot_[1]);
  interp_ankle_ddot_z_.setData(gait_percentage_vector_, ankle_path_xyz_ddot_[2]);
}
