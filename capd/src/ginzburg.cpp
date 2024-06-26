#include "capd/capdlib.h"

using namespace capd;
using namespace std;
using capd::autodiff::Node;

// Generic version of the vector field
void vectorField(Node xi, Node in[], int /*dimIn*/, Node out[], int /*dimOut*/, Node params[], int /*noParams*/)
{
  Node omega = params[0];
  Node sigma = params[1];
  Node delta = params[2];
  Node d = params[3];

  Node a = in[0];
  Node b = in[1];
  Node alpha = in[2];
  Node beta = in[3];
  Node kappa = in[4];
  Node epsilon = in[5];

  Node a2b2_sigma = exp(sigma * log((a^2) + (b^2)));

  Node F1 = -(d - 1) / xi * (alpha + epsilon * beta) +
    kappa * xi * beta +
    kappa / sigma * b +
    omega * a -
    a2b2_sigma * a +
    delta * a2b2_sigma * b;

  Node F2 = -(d - 1) / xi * (beta - epsilon * alpha) -
    kappa * xi * alpha -
    kappa / sigma * a +
    omega * b -
    a2b2_sigma * b -
    delta * a2b2_sigma * a;

  Node one_p_epsilon2 = 1 + (epsilon^2);

  out[0] = alpha;
  out[1] = beta;
  out[2] = (F1 - epsilon * F2) / one_p_epsilon2;
  out[3] = (epsilon * F1 + F2) / one_p_epsilon2;
  out[4] = 0 * kappa;
  out[5] = 0 * epsilon;
}

// Vector field handling the case d == 1
void vectorField_d1(Node xi, Node in[], int /*dimIn*/, Node out[], int /*dimOut*/, Node params[], int /*noParams*/)
{
  Node omega = params[0];
  Node sigma = params[1];
  Node delta = params[2];

  Node a = in[0];
  Node b = in[1];
  Node alpha = in[2];
  Node beta = in[3];
  Node kappa = in[4];
  Node epsilon = in[5];

  Node a2b2_sigma = exp(sigma * log((a^2) + (b^2)));

  Node F1 = kappa * xi * beta +
    kappa / sigma * b +
    omega * a -
    a2b2_sigma * a +
    delta * a2b2_sigma * b;

  Node F2 = -kappa * xi * alpha -
    kappa / sigma * a +
    omega * b -
    a2b2_sigma * b -
    delta * a2b2_sigma * a;

  Node one_p_epsilon2 = 1 + (epsilon^2);

  out[0] = alpha;
  out[1] = beta;
  out[2] = (F1 - epsilon * F2) / one_p_epsilon2;
  out[3] = (epsilon * F1 + F2) / one_p_epsilon2;
  out[4] = 0 * kappa;
  out[5] = 0 * epsilon;
}

// Vector field optimized for the case d == 1, omega == 1 and delta == 0
void vectorField_d1_optimized(Node xi, Node in[], int /*dimIn*/, Node out[], int /*dimOut*/, Node params[], int /*noParams*/)
{
  Node sigma = params[0];

  Node a = in[0];
  Node b = in[1];
  Node alpha = in[2];
  Node beta = in[3];
  Node kappa = in[4];
  Node epsilon = in[5];

  Node a2b2_sigma_m1 = exp(sigma * log((a^2) + (b^2))) - 1;
  Node kappa_xi = kappa * xi;
  Node kappa_div_sigma = kappa / sigma;

  Node F1 = kappa_xi * beta + kappa_div_sigma * b - a2b2_sigma_m1 * a;

  Node F2 = -kappa_xi * alpha - kappa_div_sigma * a - a2b2_sigma_m1 * b;

  Node one_p_epsilon2 = 1 + (epsilon^2);

  out[0] = alpha;
  out[1] = beta;
  out[2] = (F1 - epsilon * F2) / one_p_epsilon2;
  out[3] = (epsilon * F1 + F2) / one_p_epsilon2;
  out[4] = 0 * kappa;
  out[5] = 0 * epsilon;
}

// Vector field optimized for the case d == 3, omega == 1, sigma == 1
// and delta == 0
void vectorField_d3_optimized(Node xi, Node in[], int /*dimIn*/, Node out[], int /*dimOut*/, Node* /*params*/, int /*noParams*/)
{
  Node a = in[0];
  Node b = in[1];
  Node alpha = in[2];
  Node beta = in[3];
  Node kappa = in[4];
  Node epsilon = in[5];

  Node a2b2_m1 = (a^2) + (b^2) - 1;
  Node kappa_xi = kappa * xi;

  Node F1 = -2 * (alpha + epsilon * beta) / xi +
    kappa_xi * beta +
    kappa * b -
    a2b2_m1 * a;

  Node F2 = -2 * (beta - epsilon * alpha) / xi -
    kappa_xi * alpha -
    kappa * a -
    a2b2_m1 * b;

  Node one_p_epsilon2 = 1 + (epsilon^2);

  out[0] = alpha;
  out[1] = beta;
  out[2] = (F1 - epsilon * F2) / one_p_epsilon2;
  out[3] = (epsilon * F1 + F2) / one_p_epsilon2;
  out[4] = 0 * kappa;
  out[5] = 0 * epsilon;
}

// Vector field optimized for the case d == 3, omega == 1, sigma == 1,
// epsilon = 0 and delta == 0
void vectorField_d3_optimized_epsilon_0(Node xi, Node in[], int /*dimIn*/, Node out[], int /*dimOut*/, Node* /*params*/, int /*noParams*/)
{
  Node a = in[0];
  Node b = in[1];
  Node alpha = in[2];
  Node beta = in[3];
  Node kappa = in[4];
  Node epsilon = in[5];

  Node a2b2_m1 = (a^2) + (b^2) - 1;
  Node kappa_xi = kappa * xi;

  Node F1 = -2 * alpha / xi +
    kappa_xi * beta +
    kappa * b -
    a2b2_m1 * a;

  Node F2 = -2 * beta / xi -
    kappa_xi * alpha -
    kappa * a -
    a2b2_m1 * b;

  out[0] = alpha;
  out[1] = beta;
  out[2] = F1;
  out[3] = F2;
  out[4] = 0 * kappa;
  out[5] = 0 * epsilon;
}

int main()
{
  cout.precision(17); // Enough to exactly recover Float64 values
  cerr.precision(17); // Enough to exactly recover Float64 values

  // Read initial value
  IVector u0(6);
  cin >> u0[0] >> u0[1] >> u0[2] >> u0[3];

  // Read parameter values
  int d;
  interval kappa, epsilon, omega, sigma, delta;
  cin >> d;
  cin >> kappa;
  cin >> epsilon;
  cin >> omega;
  cin >> sigma;
  cin >> delta;

  u0[4] = kappa;
  u0[5] = epsilon;

  // Read time span
  interval T0, T1;
  cin >> T0 >> T1;

  // Read flags for if to output Jacobian and for if the Jacobian
  // should be with respect to epsilon instead of kappa.
  int output_jacobian;
  int jacobian_epsilon;
  cin >> output_jacobian;
  cin >> jacobian_epsilon;

  // Read tolerance to use
  double tol;
  cin >> tol;

  // Create the vector field and the parameters
  int dim = 6;
  IMap vf;

  if (d == 1 && omega == 1 && delta == 0) {
    vf = IMap(vectorField_d1_optimized, dim, dim, 1);
    vf.setParameter(0, sigma);
  } else if (d == 1) {
    vf = IMap(vectorField_d1, dim, dim, 3);
    vf.setParameter(0, omega);
    vf.setParameter(1, sigma);
    vf.setParameter(2, delta);
  } else if (d == 3 && omega == 1 && sigma == 1 && delta == 0 && epsilon == 0 && !jacobian_epsilon) {
      vf = IMap(vectorField_d3_optimized_epsilon_0, dim, dim, 0);
  } else if (d == 3 && omega == 1 && sigma == 1 && delta == 0) {
    vf = IMap(vectorField_d3_optimized, dim, dim, 0);
  } else {
    vf = IMap(vectorField, dim, dim, 4);
    vf.setParameter(0, omega);
    vf.setParameter(1, sigma);
    vf.setParameter(2, delta);
    vf.setParameter(3, interval(d));
  }

  // Create the solver and the time map
  IOdeSolver solver(vf, 20);

  solver.setAbsoluteTolerance(tol);
  solver.setRelativeTolerance(tol);

  ITimeMap timeMap(solver);

  try {
    if (output_jacobian) {
      // Define a representation of the initial value
      C1HORect2Set s(u0, T0);

      // Solve the system
      IVector result = timeMap(T1, s);

      IMatrix m = (IMatrix)(s);

      if (jacobian_epsilon) {
          for (int i = 0; i < 4; i++)
              // Only print derivatives of u[1], ..., u[4]
              for (int j = 0; j < 4; j++)
                  cout << m[j][i] << endl;

          // Derivative w.r.t. epsilon
          for (int j = 0; j < 4; j++)
              cout << m[j][5] << endl;
      } else {
          for (int i = 0; i < 5; i++)
              // Only print derivatives of u[1], ..., u[4]
              for (int j = 0; j < 4; j++)
                  cout << m[j][i] << endl;
      }
    } else {
      // Define a doubleton representation of the initial value
      C0HORect2Set s(u0, T0);

      // Solve the system
      IVector result = timeMap(T1, s);

      for (int i = 0; i < 4; i++)
        cout << result[i] << endl;

    }
  } catch(exception& e) {
    cout << "\n\nException caught!\n" << e.what() << endl << endl;
  }
} // END
