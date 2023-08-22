#include "capd/capdlib.h"

using namespace capd;
using namespace std;
using capd::autodiff::Node;

// Vector field in the case d == 1
void vectorField_d1(Node t, Node in[], int /*dimIn*/, Node out[], int /*dimOut*/, Node params[], int /*noParams*/)
{
  Node omega = params[0];
  Node sigma = params[1];
  Node epsilon = params[2];
  Node delta = params[3];

  Node a = in[0];
  Node b = in[1];
  Node alpha = in[2];
  Node beta = in[3];
  Node kappa = in[4];

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

  out[0] = in[2];
  out[1] = in[3];
  out[2] = F1 - epsilon * F2;
  out[3] = epsilon * F1 + F2;
  out[4] = 0 * kappa;
}

// Vector field in the case d != 1. The only difference is the
// definition of F1 and F2
void vectorField(Node t, Node in[], int /*dimIn*/, Node out[], int /*dimOut*/, Node params[], int /*noParams*/)
{
  Node omega = params[0];
  Node sigma = params[1];
  Node epsilon = params[2];
  Node delta = params[3];
  Node d = params[4];

  Node a = in[0];
  Node b = in[1];
  Node alpha = in[2];
  Node beta = in[3];
  Node kappa = in[4];

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

  out[0] = in[2];
  out[1] = in[3];
  out[2] = F1 - epsilon * F2;
  out[3] = epsilon * F1 + F2;
  out[4] = 0 * kappa;
}

int main()
{
  cout.precision(17); // Enough to exactly recover Float64 values

  IVector u0(5);

  // Read initial value
  cin >> u0[0] >> u0[1] >> u0[2] >> u0[3];

  // Read parameter values
  int d;
  interval omega, sigma, epsilon, delta, kappa;
  cin >> d;
  cin >> omega;
  cin >> sigma;
  cin >> epsilon;
  cin >> delta;
  cin >> kappa;

  u0[4] = kappa;

  // Read time span
  interval T0, T1;
  cin >> T0 >> T1;

  // Create the vector field and the parameters
  int dim = 5, noParam = 4;
  IMap vf;

  if (d == 1)
    vf = IMap(vectorField_d1, dim, dim, noParam);
  else
    vf = IMap(vectorField, dim, dim, noParam + 1);

  vf.setParameter(0, omega);
  vf.setParameter(1, sigma);
  vf.setParameter(2, epsilon);
  vf.setParameter(3, delta);

  if (d != 1)
    vf.setParameter(4, interval(d));

  // Create the solver and the time map
  IOdeSolver solver(vf, 20);
  solver.setAbsoluteTolerance(1e-10);
  solver.setRelativeTolerance(1e-10);

  ITimeMap timeMap(solver);

  // Define a doubleton representation of the initial value
  // TODO: Read initial condition for variational equation
  C1HORect2Set s(u0, T0);

  try{
    // Solve the system
    IVector result = timeMap(T1, s);

    for (int i = 0; i < 4; i++)
      cout << result[i] << endl;

    IMatrix m = (IMatrix)(s);

    for (int i = 0; i < 5; i++)
      for (int j = 0; j < 5; j++)
        cout << m[j][i] << endl;

  }catch(exception& e)
  {
    cout << "\n\nException caught!\n" << e.what() << endl << endl;
  }
} // END
