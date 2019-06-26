/*
 * Copyright 2015-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once

#include <mc_rbdyn/RobotModule.h>

#include <mc_rbdyn_urdf/urdf.h>

#include "api.h"

namespace mc_robots
{

struct MC_ROBOTS_DLLAPI PolarisRangerEgressRobotModule : public mc_rbdyn::RobotModule
{
public:
  PolarisRangerEgressRobotModule();

protected:
  std::map<std::string, std::pair<std::string, std::string>> getConvexHull(
      const std::map<std::string, std::pair<std::string, std::string>> & files) const;

  void readUrdf(const std::vector<std::string> & filteredLinks);

  std::map<std::string, std::vector<double>> halfSittingPose(const rbd::MultiBody & mb) const;

  std::vector<std::map<std::string, std::vector<double>>> nominalBounds(const mc_rbdyn_urdf::Limits & limits) const;

  std::map<std::string, std::pair<std::string, std::string>> stdCollisionsFiles(const rbd::MultiBody & mb) const;

  virtual const std::map<std::string, std::pair<std::string, std::string>> & convexHull() const override;

  virtual const std::vector<std::map<std::string, std::vector<double>>> & bounds() const override;

  virtual const std::map<std::string, std::vector<double>> & stance() const override;

public:
  std::vector<std::string> virtualLinks;
  std::vector<std::string> gripperLinks;
  std::map<std::string, std::vector<double>> halfSitting;
  std::map<std::string, std::vector<sva::PTransformd>> visual_tfs;
  mc_rbdyn_urdf::Limits limits;
};

} // namespace mc_robots

ROBOT_MODULE_DEFAULT_CONSTRUCTOR("PolarisRangerEgress", mc_robots::PolarisRangerEgressRobotModule)
