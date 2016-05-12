#pragma once

#include <mc_rbdyn/robot.h>


#include <mc_control/generic_gripper.h>
#include <mc_control/msg/QPResult.h>
#include <mc_solver/QPSolver.h>
#include <mc_solver/CollisionsConstraint.h>
#include <mc_solver/ContactConstraint.h>
#include <mc_solver/DynamicsConstraint.h>
#include <mc_solver/KinematicsConstraint.h>

#include <Tasks/QPTasks.h>

#include <json/json.h>

namespace mc_rbdyn
{
  struct Contact;
}

#include <mc_control/api.h>

namespace mc_control
{

/** \class ControllerResetData
 * \brief Contains information allowing the controller to start smoothly from
 * the current state of the robot
 * \note
 * For now, this only contains the state of the robot (free flyer and joints state)
 */
struct MC_CONTROL_DLLAPI ControllerResetData
{
  /** Contains free flyer + joints state information */
  const std::vector< std::vector<double> > & q;
};

struct MCGlobalController;

/** \class MCController
 * \brief MCController is the base class to implement all controllers. It
 * assumes that at least two robots are provided. The first is considered as the
 * "main" robot. Some common constraints and a posture task are defined (but not
 * added to the solver) for this robot
 */
struct MC_CONTROL_DLLAPI MCController
{
  friend struct MCGlobalController;
public:
  virtual ~MCController();
  /** This function is called at each time step of the process driving the robot
   * (i.e. simulation or robot's controller). This function is the most likely
   * to be overriden for complex controller behaviours.
   * \return True if the solver succeeded, false otherwise
   * \note
   * This is meant to run in real-time hence some precaution should apply (e.g.
   * no i/o blocking calls, no thread instantiation and such)
   *
   * \note
   * The default implementation does the bare minimum (i.e. call run on QPSolver)
   * It is recommended to use it in your override.
   */
  virtual bool run();

  /** Gives access to the result of the QP execution
   * \param t Unused at the moment
   */
  virtual const QPResultMsg & send(const double & t);

  /** Reset the controller with data provided by ControllerResetData. This is
   * called at two possible points during a simulation/live execution:
   *   1. Actual start
   *   2. Switch from a previous (MCController-like) controller
   * In the first case, the data comes from the simulation/controller. In the
   * second case, the data comes from the previous MCController instance.
   * \param reset_data Contains information allowing to reset the controller
   * properly
   * \note
   * The default implementation reset the main robot's state to that provided by
   * reset_data (with a null speed/acceleration). It maintains the contacts as
   * they were set by the controller previously.
   */
  virtual void reset(const ControllerResetData & reset_data);

  /** Pass force sensors information to the controller if available at the
   * simulation/controller level
   * \param wrenches Force/Torque sensors information provided by the
   * simulation/controller
   */
  virtual void setWrenches(const std::map<std::string, sva::ForceVecd>& wrenches);

  /** Get the current wrenches information
   * \return A vector of sva::ForceVecd representing force/torque pairs.
   */
  const std::map<std::string, sva::ForceVecd> & getWrenches();

  /** Get the current encoder values provided by simulation/low-level controller
   * \return A vector of encoder values ordered by the controller/simulator
   */
  const std::vector<double> & getEncoderValues();

  /** Get the current joint torques provided by the low-level controller
   * \return A vector of joint torques ordered according to RobotModule::ref_joint_order()
   */
  const std::vector<double> & getJointTorques();

  /** Get the sensor orientation
   * \return The sensor orientation if provided, Eigen::Vector3d::Zero() otherwise
   */
  const Eigen::Vector3d & getSensorOrientation();

  /** Get the sensor angular velocity
   * \return The sensor angular velocity if provided, Eigen::Vector3d::Zero() otherwise
   */
  const Eigen::Vector3d & getSensorVelocity();

  /** Get the sensor acceleration
   * \return The sensor acceleration if provided, Eigen::Vector3d::Zero() otherwise
   */
  const Eigen::Vector3d & getSensorAcceleration();

  /** Return the main robot (first robot provided in the constructor
   * \anchor mc_controller_robot_const_doc
   */
  virtual const mc_rbdyn::Robot & robot() const;

  /** Return the env "robot"
   * \note
   * In multi-robot scenarios, the env robot is either:
   *   1. The first robot with zero dof
   *   2. The last robot provided at construction
   * \anchor mc_controller_env_const_doc
   */
  virtual const mc_rbdyn::Robot & env() const;

  /** Return the mc_rbdyn::Robots controlled by this controller
   * \anchor mc_controller_robots_const_doc
   */
  virtual const mc_rbdyn::Robots & robots() const;

  /** Non-const variant of \ref mc_controller_robots_const_doc "robots()" */
  virtual mc_rbdyn::Robots & robots();

  /** Non-const variant of \ref mc_controller_robot_const_doc "robot()" */
  virtual mc_rbdyn::Robot & robot();

  /** Non-const variant of \ref mc_controller_env_const_doc "env()" */
  virtual mc_rbdyn::Robot & env();

  /** Return the mc_solver::QPSolver instance attached to this controller
   * \anchor mc_controller_qpsolver_const_doc
   */
  const mc_solver::QPSolver & solver() const;

  /** Non-const variant of \ref mc_controller_qpsolver_const_doc "solver()" */
  mc_solver::QPSolver & solver();

  /** Set a joint position to the desired value
   * \param jname Name of the joint to control
   * \param pos Desired position (radians)
   * \return True if jname is valid, false otherwise
   * \note
   * No control is made on the value of pos to ensure it is within the joints'
   * limits of jname. It is assumed that this is done via the controller own
   * constraint set
   *
   * \note
   * The default implementation only works on the main robot.
   */
  virtual bool set_joint_pos(const std::string & jname, const double & pos);

  /** Get a joint position
   * \param jname Name of the desired joint
   * \param position position of jname joint (radians)
   * \return True if jname is valid, false otherwise
   * \note
   * Due to the overhead of checks on joint validity, it is not recommended
   * to use this function to repeatedly access a specific joint value.
   *
   * \note
   * The default implementation only works on the main robot.
   */
  bool get_joint_pos(const std::string & jname, double & pos);

  /** Change the currently controlled end-effector
   * \param name End of the name effector
   * \return False if the controller does not implement this kind of control or
   * if name is not a end-effector, true otherwise
   */
  virtual bool change_ef(const std::string & name);

  /** Move the currently controlled end-effector
   * \param t Translation amount
   * \param m Rotation applied to current orientation
   * \return False if the controller does not implement this control, true
   * otherwise
   */
  virtual bool move_ef(const Eigen::Vector3d & t, const Eigen::Matrix3d & m);

  /** Move the CoM
   * \param t CoM translation
   * \return False if the controller does not implement this control, true
   */
  virtual bool move_com(const Eigen::Vector3d & t);

  /** Trigger next step in a FSM controller
   * \return False if the controller does not implement this or if the switch
   * cannot happen, true otherwise
   */
  virtual bool play_next_stance();

  /** Driving service
   * \param wheel Wheel rotation angle
   * \param ankle Ankle angle (related to acceleration)
   * \param pan Head pan angle
   * \param tilt Head tilt angle
   * \return False if the controller does not implement this, true otherwise
   */
  virtual bool driving_service(double wheel, double ankle, double pan, double tilt);

  /** Generic message passing interface, write-only version
   * \param msg A message passed to the controller
   * \return True if the controller was able to do something out of msg, false
   * otherwise
   * \note
   * The default implementation does nothing and always return false.
   */
  virtual bool read_msg(std::string & msg);

  /** Generic message passing interface, read/write version
   * \param msg A message passed to the controller
   * \param out A message passed back to the caller
   * \return True if the controller was able to do something out of msg, false
   * otherwise
   * \note
   * The default implementation does nothing and always return false.
   */
  virtual bool read_write_msg(std::string & msg, std::string & out);

  /** Logging interface: add information to the csv header. The csv entries are
   * separated by a ";"
   * \param os Datastream containing the default header
   * \note
   * This is called once (when the controller is started)
   *
   * \note
   * The default implementation does nothing.
   */
  virtual std::ostream& log_header(std::ostream & os);

  /** Logging interface: add data to the csv file. The csv entries are
   * separated by a ";"
   * \param os Datastream containing the default header
   * \note
   * This is called at each timestep after run() invokation
   *
   * \note
   * The default implementation does nothing.
   */
  virtual std::ostream& log_data(std::ostream & os);

  /** Returns a list of robots supported by the controller.
   * \return Vector of supported robots designed by name (as returned by
   * RobotModule::name())
   * \note
   * Default implementation returns an empty list which indicates that all
   * robots are supported.
   */
  virtual std::vector<std::string> supported_robots() const;
protected:
  /** Builds a controller base with an empty environment
   * \param robot Pointer to the main RobotModule
   * \param dt Controller timestep
   * your controller
   */
  MCController(std::shared_ptr<mc_rbdyn::RobotModule> robot, double dt);

  /** Builds a multi-robot controller base
   * \param robots Collection of robot modules used by the controller
   * \param dt Timestep of the controller
   */
  MCController(const std::vector<std::shared_ptr<mc_rbdyn::RobotModule>> & robot_modules, double dt);
protected:
  /** Encoder values provided by the low-level controller */
  std::vector<double> encoderValues;
  /** Joint torques provided by the low-level controller */
  std::vector<double> jointTorques;
  /** Force/Torque sensors */
  std::map<std::string, sva::ForceVecd> wrenches;
  /** Robot orientation provided by sensors */
  Eigen::Vector3d sensorOri;
  /** Robot acceleration provided by sensors */
  Eigen::Vector3d sensorAcc;
  /** Robot angular velocity provided by sensors */
  Eigen::Vector3d sensorVel;
  /** QP solver */
  std::shared_ptr<mc_solver::QPSolver> qpsolver;
public:
  /** Controller timestep */
  const double timeStep;
  /** Reference joint order see mc_rbdyn::RobotModule */
  std::vector<std::string> ref_joint_order;
  /** Grippers */
  std::map<std::string, std::shared_ptr<mc_control::Gripper>> grippers;
  /** Contact constraint for the main robot */
  mc_solver::ContactConstraint contactConstraint;
  /** Dynamics constraints for the main robot */
  mc_solver::DynamicsConstraint dynamicsConstraint;
  /** Kinematics constraints for the main robot */
  mc_solver::KinematicsConstraint kinematicsConstraint;
  /** Self collisions constraint for the main robot */
  mc_solver::CollisionsConstraint selfCollisionConstraint;
  /** Posture task for the main robot */
  std::shared_ptr<tasks::qp::PostureTask> postureTask;
};

}

#ifdef WIN32
  #define CONTROLLER_MODULE_API __declspec(dllexport)
#else
  #define CONTROLLER_MODULE_API
#endif

/** Provides a handle to construct the controller with Json config */
#define CONTROLLER_CONSTRUCTOR(NAME, TYPE)\
extern "C"\
{\
  CONTROLLER_MODULE_API const char * CLASS_NAME() { return NAME; }\
  CONTROLLER_MODULE_API void destroy(mc_control::MCController * ptr) { delete ptr; }\
  CONTROLLER_MODULE_API mc_control::MCController * create(const std::shared_ptr<mc_rbdyn::RobotModule> & robot, const double & dt, const Json::Value & conf) { return new TYPE(robot, dt, conf); }\
}

/** Provides a handle to construct a generic controller */
#define SIMPLE_CONTROLLER_CONSTRUCTOR(NAME, TYPE)\
extern "C"\
{\
  CONTROLLER_MODULE_API const char * CLASS_NAME() { return NAME; }\
  CONTROLLER_MODULE_API void destroy(mc_control::MCController * ptr) { delete ptr; }\
  CONTROLLER_MODULE_API mc_control::MCController * create(const std::shared_ptr<mc_rbdyn::RobotModule> & robot, const double & dt, const Json::Value &) { return new TYPE(robot, dt); }\
}
