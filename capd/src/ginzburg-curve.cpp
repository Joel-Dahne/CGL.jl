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

// Print an enclosure of the curve on the given domain, also prints
// enclosure of the first and second derivative of the squared
// absolute value.
// In case either the first or second derivative of the squared
// absolute value contains zero then bisect the domain in two and call
// this recursively. Stops after a maximum of 5 bisections.
// On the first call (depth == 0) it checks if the first or second
// derivative of the squared absolute value contains zero when
// evaluating the curve on the right endpoint. If this is the case it
// doesn't try to bisect any further, since the requirement could
// never be fulfilled anyway.
void print_curve(const IOdeSolver::SolutionCurve &curve, interval domain, interval prevTime, int depth) {
    if (depth == 0) {
        // Evaluate at right endpoint to see if the requirement is
        // satisfied there.
        IVector v = curve(domain.rightBound());
        IVector dv = curve.timeDerivative(domain.rightBound());

        interval abs2_Q_derivative = 2 * (v[2] * v[0] + v[3] * v[1]);
        interval abs2_Q_derivative2 = 2 * (dv[2] * v[0] + v[2] * v[2] + dv[3] * v[1] + v[3] * v[3]);

        bool abs2_Q_derivative_non_zero = !interval(0).subset(abs2_Q_derivative);
        bool abs2_Q_derivative2_non_zero = !interval(0).subset(abs2_Q_derivative2);

        if (!abs2_Q_derivative_non_zero && !abs2_Q_derivative2_non_zero) {
            print_curve(curve, domain, prevTime, 5);
            return;
        }
    }

    // Here we evaluated curve at the interval domain. v will
    // contain rigorous bound for the trajectory for this time
    // interval.
    IVector v = curve(domain);
    IVector dv = curve.timeDerivative(domain);

    interval abs2_Q_derivative = 2 * (v[2] * v[0] + v[3] * v[1]);
    interval abs2_Q_derivative2 = 2 * (dv[2] * v[0] + v[2] * v[2] + dv[3] * v[1] + v[3] * v[3]);

    bool abs2_Q_derivative_non_zero = !interval(0).subset(abs2_Q_derivative);
    bool abs2_Q_derivative2_non_zero = !interval(0).subset(abs2_Q_derivative2);

    if (abs2_Q_derivative_non_zero || abs2_Q_derivative2_non_zero || depth >= 5) {
        cout << prevTime + domain; // xi value
        cout << ";" << v[0]; // a
        cout << ";" << v[1]; // b
        cout << ";" << v[2]; // alpha
        cout << ";" << v[3]; // beta
        cout << ";" << dv[2]; // second derivative of a
        cout << ";" << dv[3]; // second derivative of b
        cout << ";" << abs2_Q_derivative; // Derivative of abs(Q)^2
        cout << ";" << abs2_Q_derivative2 << endl; // Second derivative of abs(Q)^2 / 2
    } else {
        interval domain_left =
            interval(domain.leftBound(), (domain.leftBound() + domain.rightBound()) / 2);
        interval domain_right =
            interval((domain.leftBound() + domain.rightBound()) / 2, domain.rightBound());

        print_curve(curve, domain_left, prevTime, depth + 1);
        print_curve(curve, domain_right, prevTime, depth + 1);
    }
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
  } else if (d == 3 && omega == 1 && sigma == 1 && delta == 0 && epsilon == 0) {
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

  timeMap.stopAfterStep(true);

  interval prevTime(T0);

  // Define a representation of the initial value
  C0HORect2Set s(u0, T0);

  try {
      do {
	  timeMap(T1, s);

	  interval stepMade = solver.getStep();

	  // This is how we can extract an information about the
	  // trajectory between time steps. The type CurveType is a
	  // function defined on the interval [0,stepMade]. It can be
	  // evaluated at a point (or interval). The curve can be also
	  // differentiated wrt to time. We can also extract from it the
	  // 1-st order derivatives wrt.
	  const IOdeSolver::SolutionCurve& curve = solver.getCurve();
	  interval domain = interval(0,1) * stepMade;

          print_curve(curve, domain, prevTime, 0);

	  prevTime = timeMap.getCurrentTime();
      } while (!timeMap.completed());
  } catch(exception& e) {
    cout << "\n\nException caught!\n" << e.what() << endl << endl;
  }
} // END
