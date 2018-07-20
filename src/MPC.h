#ifndef MPC_H
#define MPC_H

#include <vector>
#include "Eigen-3.3/Eigen/Core"

using namespace std;

// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
const double Lf = 2.67;

// TODO: Set the timestep length and duration
const size_t N = 12;
const double dt = 0.1;
const double MAX_VELOCITY = 80;  // mph

// TODO: Set the number of model variables (includes both states and inputs).
// For example: If the state is a 4 element vector, the actuators is a 2
// element vector and there are 10 timesteps. The number of variables is:
//
// 4 * 10 + 2 * 9
const int NUM_OF_STATES = 6;
const size_t NUM_VARS = N * NUM_OF_STATES + (N - 1) * 2;
// TODO: Set the number of constraints
const size_t NUM_CONSTRAINTS = N * NUM_OF_STATES;

const size_t X_START = 0;
const size_t Y_START = X_START + N;
const size_t PSI_START = Y_START + N;
const size_t V_START = PSI_START + N;
const size_t CTE_START = V_START + N;
const size_t EPSI_START = CTE_START + N;
const size_t DELTA_START = EPSI_START + N;
const size_t A_START = DELTA_START + N - 1;

// weights for cost computations
const double WGT_CTE = 1100.0;
const double WGT_EPSI = 1100.0;
const double WGT_DELTA = 50.0;
const double WGT_ACCEL = 50.0;
const double WGT_DDELTA = 210000.0; // weight cost for difference between consecutive steering actuations
const double WGT_DACCEL = 1100.0;   // weight cost for difference between consecutive acceleration actuations

class MPC {
 public:
  MPC();

  virtual ~MPC();

  // Solve the model given an initial state and polynomial coefficients.
  // Return the first actuatotions.
  vector<double> Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs);
};

#endif /* MPC_H */
