#pragma once

#include <mc_rbdyn/Collision.h>
#include <mc_rbdyn/robot.h>

#include <map>
#include <vector>

#include <mc_rbdyn/api.h>
#include <mc_rbdyn_urdf/urdf.h>

/* This is an interface designed to provide additionnal information about a robot */

namespace mc_rbdyn
{

struct Robot;

/* TODO Functions are declared const here but most implementations will likely not respect the constness */
struct MC_RBDYN_DLLAPI RobotModule
{
  /*! Holds necessary information to create a gripper */
  struct Gripper
  {
    /*! Gripper's name */
    std::string name;
    /*! Active joints in the gripper */
    std::vector<std::string> joints;
    /*! Whether the limits should be reversed, see mc_control::Gripper */
    bool reverse_limits;
  };

  RobotModule(const std::string & path, const std::string & name)
  : RobotModule(path, name, path + "/urdf/" + name + ".urdf")
  {}

  RobotModule(const std::string & path, const std::string & name, const std::string & urdf_path)
  : path(path), name(name), urdf_path(urdf_path),
    rsdf_dir(path + "/rsdf/" + name),
    calib_dir(path + "/calib/" + name)
  {}

  virtual ~RobotModule() {}

  /* If implemented, returns limits in this order:
      - joint limits (lower/upper)
      - velocity limits (lower/upper)
      - torque limits (lower/upper)
    Limits are maps of (jointName, limit)
  */
  virtual const std::vector< std::map<std::string, std::vector<double> > > & bounds() const { return _bounds; }

  /* return the initial configuration of the robot as a map (jointName, jointPos) */
  virtual const std::map<std::string, std::vector<double> > & stance() const { return _stance; }

  /* return a map (name, (bodyName, PolyhedronURL)) */
  virtual const std::map<std::string, std::pair<std::string, std::string> > & convexHull() const { return _convexHull; }

  /* return a map (name, (bodyName, STPBVURL)) */
  virtual const std::map<std::string, std::pair<std::string, std::string> > & stpbvHull() const { return _stpbvHull; }

  /* return a map (bodyName, sva::PTransformd) */
  virtual const std::map<std::string, sva::PTransformd> & collisionTransforms() const { return _collisionTransforms; }

  /* return flexibilities */
  virtual const std::vector<Flexibility> & flexibility() const { return _flexibility; }

  /* return force sensors */
  virtual const std::vector<ForceSensor> & forceSensors() const { return _forceSensors; }

  const std::string & accelerometerBody() const { return _accelerometerBody; }

  virtual const Springs & springs() const { return _springs; }

  /** Return a minimal self-collision set */
  virtual const std::vector<mc_rbdyn::Collision> & minimalSelfCollisions() { return _minimalSelfCollisions; }
  /** Return a broader set of the most common self-collisions */
  virtual const std::vector<mc_rbdyn::Collision> & commonSelfCollisions() { return _commonSelfCollisions; }

  /** Return a map of gripper. Keys represents the gripper name. Values indicate the active joints in the gripper. */
  virtual const std::vector<Gripper> & grippers() { return _grippers; }

  /** Return the reference (native controller) joint order of the robot */
  virtual const std::vector<std::string> & ref_joint_order() { return _ref_joint_order; }

  std::string path;
  std::string name;
  std::string urdf_path;
  std::string rsdf_dir;
  std::string calib_dir;
  rbd::MultiBody mb;
  rbd::MultiBodyConfig mbc;
  rbd::MultiBodyGraph mbg;
  std::vector< std::map<std::string, std::vector<double> > > _bounds;
  std::map<std::string, std::vector<double> > _stance;
  std::map<std::string, std::pair<std::string, std::string> > _convexHull;
  std::map<std::string, std::pair<std::string, std::string> > _stpbvHull;
  std::map<std::string, std::vector<mc_rbdyn_urdf::Visual>> _visual;
  std::map<std::string, sva::PTransformd> _collisionTransforms;
  std::vector<Flexibility> _flexibility;
  std::vector<ForceSensor> _forceSensors;
  std::string _accelerometerBody;
  Springs _springs;
  std::vector<mc_rbdyn::Collision> _minimalSelfCollisions;
  std::vector<mc_rbdyn::Collision> _commonSelfCollisions;
  std::vector<Gripper> _grippers;
  std::vector<std::string> _ref_joint_order;
};

} // namespace mc_rbdyn

/* Set of macros to assist with the writing of a RobotModule */

#ifdef WIN32
#define ROBOT_MODULE_API __declspec(dllexport)
#else
#define ROBOT_MODULE_API
#endif

/*! ROBOT_MODULE_COMMON
 * Declare a destroy symbol and CLASS_NAME symbol
 * Constructor should be declared by the user
*/
#define ROBOT_MODULE_COMMON(NAME)\
  ROBOT_MODULE_API const char * CLASS_NAME() { return NAME; }\
  ROBOT_MODULE_API void destroy(mc_rbdyn::RobotModule * ptr) { delete ptr; }

/*! ROBOT_MODULE_DEFAULT_CONSTRUCTOR
 * Declare an external symbol for creation using a default constructor
 * Also declare destruction symbol
 * Exclusive of ROBOT_MODULE_CANONIC_CONSTRUCTOR
*/
#define ROBOT_MODULE_DEFAULT_CONSTRUCTOR(NAME, TYPE)\
extern "C"\
{\
  ROBOT_MODULE_COMMON(NAME)\
  ROBOT_MODULE_API mc_rbdyn::RobotModule * create() { return new TYPE(); }\
}

/*! ROBOT_MODULE_CANONIC_CONSTRUCTOR
 * Declare an external symbol for creation using the cannonical constructor (const string &, const string &)
 * Also declare destruction symbol
 * Exclusive of ROBOT_MODULE_DEFAULT_CONSTRUCTOR
*/
#define ROBOT_MODULE_CANONIC_CONSTRUCTOR(NAME, TYPE)\
extern "C"\
{\
  ROBOT_MODULE_COMMON(NAME)\
  ROBOT_MODULE_API mc_rbdyn::RobotModule * create(const std::string & path, const std::string & name) { return new TYPE(path, name); }\
}
