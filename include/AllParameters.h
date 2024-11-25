//
// Created by Florian on 25.11.2024.
//

#ifndef LAPLACIAN_ALLPARAMETERS_H
#define LAPLACIAN_ALLPARAMETERS_H
#include <deal.II/base/parameter_handler.h>

#include "FESystem.h"
#include "Geometry.h"
#include "LinearSolver.h"
#include "Materials.h"
#include "NonlinearSolver.h"
#include "time.h"

using namespace dealii;
  struct AllParameters : public FESystem,
                         public Geometry,
                         public Materials,
                         public LinearSolver,
                         public NonlinearSolver,
                         public Time

  {
    AllParameters(const std::string &input_file);

    static void declare_parameters(ParameterHandler &prm);

    void parse_parameters(ParameterHandler &prm);
  };



#endif // LAPLACIAN_ALLPARAMETERS_H
