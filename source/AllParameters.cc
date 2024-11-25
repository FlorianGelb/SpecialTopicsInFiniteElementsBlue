//
// Created by Florian on 25.11.2024.
//

#include "AllParameters.h"

using namespace dealii;
  AllParameters::AllParameters(const std::string &input_file)
  {
    ParameterHandler prm;
    declare_parameters(prm);
    prm.parse_input(input_file);
    parse_parameters(prm);
  }

  void AllParameters::declare_parameters(ParameterHandler &prm)
  {
    FESystem::declare_parameters(prm);
    Geometry::declare_parameters(prm);
    Materials::declare_parameters(prm);
    LinearSolver::declare_parameters(prm);
    NonlinearSolver::declare_parameters(prm);
    Time::declare_parameters(prm);
  }

  void AllParameters::parse_parameters(ParameterHandler &prm)
  {
    FESystem::parse_parameters(prm);
    Geometry::parse_parameters(prm);
    Materials::parse_parameters(prm);
    LinearSolver::parse_parameters(prm);
    NonlinearSolver::parse_parameters(prm);
    Time::parse_parameters(prm);
  }

