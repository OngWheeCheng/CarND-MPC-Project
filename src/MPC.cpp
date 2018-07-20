#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

typedef CPPAD_TESTVECTOR(double) Dvector;

class FG_eval {
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  FG_eval(Eigen::VectorXd coeffs) { this->coeffs = coeffs; }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {
    // TODO: implement MPC
    // `fg` a vector of the cost constraints, `vars` is a vector of variable values (state & actuators)
    // NOTE: You'll probably go back and forth between this function and
    // the Solver function below.
  unsigned int t;
    // The cost is stored is the first element of `fg`.
    // Any additions to the cost should be added to `fg[0]`.
    fg[0] = 0.0;

     // The part of the cost based on the reference state.
    for(t = 0; t < N; t++) {
      fg[0] += WGT_CTE * CppAD::pow(vars[CTE_START + t], 2);   // cross track error
      fg[0] += WGT_EPSI * CppAD::pow(vars[EPSI_START + t], 2);  // heading error
      fg[0] += CppAD::pow(vars[V_START + t] - MAX_VELOCITY, 2);  // velocity error
    }

    // Minimize the use of actuators (minimize change rate)
    for (t = 0; t < N - 1; t++) {
      fg[0] += WGT_DELTA * CppAD::pow(vars[DELTA_START + t], 2);
      fg[0] += WGT_ACCEL * CppAD::pow(vars[A_START + t], 2);
    }

    // Minimize the value gap between sequential actuations.
    // (how smooth the actuations are)
    for (t = 0; t < N - 2; t++) {
      fg[0] += WGT_DDELTA * CppAD::pow(vars[DELTA_START + t + 1] - vars[DELTA_START + t], 2);
      fg[0] += WGT_DACCEL * CppAD::pow(vars[A_START + t + 1] - vars[A_START + t], 2);
    }

    //
    // Setup Constraints
    //
    // NOTE: In this section you'll setup the model constraints.

    // Initial constraints
    //
    // We add 1 to each of the starting indices due to cost being located at
    // index 0 of `fg`.
    // This bumps up the position of all the other values.
    fg[1 + X_START] = vars[X_START];
    fg[1 + Y_START] = vars[Y_START];
    fg[1 + PSI_START] = vars[PSI_START];
    fg[1 + V_START] = vars[V_START];
    fg[1 + CTE_START] = vars[CTE_START];
    fg[1 + EPSI_START] = vars[EPSI_START];

    // The rest of the constraints
    for (t = 1; t < N; t++) {
      // The state at time t+1 .
      AD<double> x1 = vars[X_START + t];
      AD<double> y1 = vars[Y_START + t];
      AD<double> psi1 = vars[PSI_START + t];
      AD<double> v1 = vars[V_START + t];
      AD<double> cte1 = vars[CTE_START + t];
      AD<double> epsi1 = vars[EPSI_START + t];

      // The state at time t.
      AD<double> x0 = vars[X_START + t - 1];
      AD<double> y0 = vars[Y_START + t - 1];
      AD<double> psi0 = vars[PSI_START + t - 1];
      AD<double> v0 = vars[V_START + t - 1];
      AD<double> cte0 = vars[CTE_START + t - 1];
      AD<double> epsi0 = vars[EPSI_START + t - 1];

      // Only consider the actuation at time t.
      AD<double> delta0 = vars[DELTA_START + t - 1];
      AD<double> a0 = vars[A_START + t - 1];

      AD<double> f0 = coeffs[0] + coeffs[1] * x0 + coeffs[2] * CppAD::pow(x0, 2) + coeffs[3] * CppAD::pow(x0, 3);
      AD<double> psides0 = CppAD::atan(coeffs[1] + 2 * coeffs[2] * x0 + 3 * coeffs[3] * CppAD::pow(x0, 2));

      // Equations for the model:
      // x_[t] = x[t-1] + v[t-1] * cos(psi[t-1]) * dt
      // y_[t] = y[t-1] + v[t-1] * sin(psi[t-1]) * dt
      // psi_[t] = psi[t-1] + v[t-1] / Lf * delta[t-1] * dt
      // v_[t] = v[t-1] + a[t-1] * dt
      // cte[t] = f(x[t-1]) - y[t-1] + v[t-1] * sin(epsi[t-1]) * dt
      // epsi[t] = psi[t] - psides[t-1] + v[t-1] * delta[t-1] / Lf * dt
      fg[1 + X_START + t] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
      fg[1 + Y_START + t] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
      fg[1 + PSI_START + t] = psi1 - (psi0 - v0 / Lf * delta0 * dt);
      fg[1 + V_START + t] = v1 - (v0 + a0 * dt);
      fg[1 + CTE_START + t] = cte1 - ((f0 - y0) + (v0 * CppAD::sin(epsi0) * dt));
      fg[1 + EPSI_START + t] = epsi1 - ((psi0 - psides0) - v0 / Lf * delta0 * dt);
    }
  }
};

//
// MPC class definition implementation.
//
MPC::MPC() {}
MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
  //bool ok = true;
  unsigned int i;

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(NUM_VARS);
  for (i = 0; i < NUM_VARS; i++) {
    vars[i] = 0;
  }

  Dvector vars_lowerbound(NUM_VARS);
  Dvector vars_upperbound(NUM_VARS);
  // TODO: Set lower and upper limits for variables.
  // Set the initial variable values

  // Set all non-actuators upper and lower limits
  // to the max negative and positive values.
  for (i = 0; i < DELTA_START; i++ ) {
    vars_lowerbound[i] = -1.0e19;
    vars_upperbound[i] = 1.0e19;
  }

  // Steering upper and lower limits
  for (i = DELTA_START; i < A_START; i++) {
    vars_lowerbound[i] = -0.436332;
    vars_upperbound[i] = 0.436332;
  }

  // Acceleration/decceleration upper and lower limits.
  for (i = A_START; i < NUM_VARS; i++) {
    vars_lowerbound[i] = -1.0;
    vars_upperbound[i] = 1.0;
  }

  // Lower and upper limits for the constraints
  // All of these should be 0 except the initial
  // state indices.
  Dvector constraints_lowerbound(NUM_CONSTRAINTS);
  Dvector constraints_upperbound(NUM_CONSTRAINTS);
  for (i = 0; i < NUM_CONSTRAINTS; i++) {
    constraints_lowerbound[i] = 0.0;
    constraints_upperbound[i] = 0.0;
  }

  const double x = state[0];
  const double y = state[1];
  const double psi = state[2];
  const double v = state[3];
  const double cte = state[4];
  const double epsi = state[5];

  constraints_lowerbound[X_START] = x;
  constraints_lowerbound[Y_START] = y;
  constraints_lowerbound[PSI_START] = psi;
  constraints_lowerbound[V_START] = v;
  constraints_lowerbound[CTE_START] = cte;
  constraints_lowerbound[EPSI_START] = epsi;

  constraints_upperbound[X_START] = x;
  constraints_upperbound[Y_START] = y;
  constraints_upperbound[PSI_START] = psi;
  constraints_upperbound[V_START] = v;
  constraints_upperbound[CTE_START] = cte;
  constraints_upperbound[EPSI_START] = epsi;

  // Object that computes objective and constraints
  FG_eval fg_eval(coeffs);

  //
  // NOTE: You don't have to worry about these options
  //
  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  //ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Cost
  //auto cost = solution.obj_value;
  //std::cout << "Cost " << cost << std::endl;

  // TODO: Return the first actuator values. The variables can be accessed with
  // `solution.x[i]`.
  //
  // {...} is shorthand for creating a vector, so auto x1 = {1.0,2.0}
  // creates a 2 element double vector.
  vector<double> actuatorValues;

  actuatorValues.push_back(solution.x[DELTA_START]);
  actuatorValues.push_back(solution.x[A_START]);

  for (i = 0; i < N - 2; i++) {
    actuatorValues.push_back(solution.x[X_START + i + 1]);
    actuatorValues.push_back(solution.x[Y_START + i + 1]);
  }
  return actuatorValues;
}