#include "capd/capdlib.h"

using namespace capd;
using namespace std;
using capd::autodiff::Node;

// Vector field in the case d == 1
void vectorField_d1(Node t, Node in[], int /*dimIn*/, Node out[], int /*dimOut*/, Node params[], int /*noParams*/)
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

  Node F1 = kappa * t * beta +
    kappa / sigma * b +
    omega * a -
    a2b2_sigma * a +
    delta * a2b2_sigma * b;

  Node F2 = -kappa * t * alpha -
    kappa / sigma * a +
    omega * b -
    a2b2_sigma * b -
    delta * a2b2_sigma * a;

  Node one_p_epsilon2 = 1 + (epsilon^2);

  out[0] = in[2];
  out[1] = in[3];
  out[2] = (F1 - epsilon * F2) / one_p_epsilon2;
  out[3] = (epsilon * F1 + F2) / one_p_epsilon2;
  out[4] = 0 * kappa;
  out[5] = 0 * epsilon;
}

// Vector field in the case d != 1. The only difference is the
// definition of F1 and F2
void vectorField(Node t, Node in[], int /*dimIn*/, Node out[], int /*dimOut*/, Node params[], int /*noParams*/)
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

  Node F1 = -(d - 1) / t * (alpha + epsilon * beta) +
    kappa * t * beta +
    kappa / sigma * b +
    omega * a -
    a2b2_sigma * a +
    delta * a2b2_sigma * b;

  Node F2 = -(d - 1) / t * (beta - epsilon * alpha) -
    kappa * t * alpha -
    kappa / sigma * a +
    omega * b -
    a2b2_sigma * b -
    delta * a2b2_sigma * a;

  Node one_p_epsilon2 = 1 + (epsilon^2);

  out[0] = in[2];
  out[1] = in[3];
  out[2] = (F1 - epsilon * F2) / one_p_epsilon2;
  out[3] = (epsilon * F1 + F2) / one_p_epsilon2;
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

  // Create the vector field and the parameters
  int dim = 6, noParam = 3;
  IMap vf;

  if (d == 1)
    vf = IMap(vectorField_d1, dim, dim, noParam);
  else
    vf = IMap(vectorField, dim, dim, noParam + 1);

  vf.setParameter(0, omega);
  vf.setParameter(1, sigma);
  vf.setParameter(2, delta);

  if (d != 1)
    vf.setParameter(3, interval(d));

  int output_jacobian;
  int jacobian_epsilon;
  cin >> output_jacobian;
  cin >> jacobian_epsilon;

  // Create the solver and the time map
  IOdeSolver solver(vf, 20);

  // IMPROVE: Consider choosing the tolerance depending on the input.
  //double tol;
  //tol = capd::max((u0[0].rightBound() - u0[0].leftBound()) * 1e-2, 1e-10);
  //solver.setAbsoluteTolerance(tol);
  //solver.setRelativeTolerance(tol);

  solver.setAbsoluteTolerance(1e-10);
  solver.setRelativeTolerance(1e-10);

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
