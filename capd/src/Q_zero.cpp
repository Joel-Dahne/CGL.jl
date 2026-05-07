#include "capd/capdlib.h"
#include <vector>

using namespace capd;
using namespace std;
using capd::autodiff::Node;

// ===========================================================================
// VECTOR FIELD DEFINITIONS
// ===========================================================================

// Generic version of the CGL vector field as a first-order real system.
// State: in = [a, b, alpha, beta, kappa, epsilon], where a + i*b = Q and
// alpha + i*beta = Q'. Params: [omega, sigma, delta, d].
// Outputs [alpha, beta, a'', b'', 0, 0] (kappa and epsilon are kept constant).
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

// Specialization for d == 1: the -(d-1)/xi * (...) singular term vanishes.
// Params: [omega, sigma, delta].
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

// Further specialization for d == 1, omega == 1, delta == 0: the omega*a - |Q|^{2sigma}*a
// term is rewritten using (|Q|^{2sigma} - 1), and the delta terms vanish. Params: [sigma].
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

// Specialization for d == 3, omega == 1, sigma == 1, delta == 0.
// The -2/xi terms in F1 and F2 cancel the (1+epsilon^2) denominator after
// forming (F1 - epsilon*F2) and (epsilon*F1 + F2); see the inline comment
// for details. No params required.
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
  Node one_p_epsilon2 = 1 + (epsilon^2);

  // The -2*(alpha + epsilon*beta)/xi term in F1 and -2*(beta - epsilon*alpha)/xi
  // term in F2, after forming (F1 - epsilon*F2) and (epsilon*F1 + F2), combine
  // to -2*alpha*(1 + epsilon^2)/xi and -2*beta*(1 + epsilon^2)/xi respectively,
  // cancelling the (1 + epsilon^2) denominator. We factor these out to avoid
  // multiplying wide intervals F1, F2 by epsilon.
  Node a_m_eps_b = a - epsilon * b;
  Node b_p_eps_a = b + epsilon * a;

  Node G1 =  kappa_xi * (beta  + epsilon * alpha) + kappa * b_p_eps_a - a2b2_m1 * a_m_eps_b;
  Node G2 = -kappa_xi * (alpha - epsilon * beta)  - kappa * a_m_eps_b - a2b2_m1 * b_p_eps_a;

  out[0] = alpha;
  out[1] = beta;
  out[2] = -2 * alpha / xi + G1 / one_p_epsilon2;
  out[3] = -2 * beta  / xi + G2 / one_p_epsilon2;
  out[4] = 0 * kappa;
  out[5] = 0 * epsilon;
}

// Specialization for d == 3, omega == 1, sigma == 1, epsilon == 0, delta == 0.
// With epsilon == 0 there are no cross terms and no denominator, giving a
// simpler and tighter enclosure than vectorField_d3_optimized. No params required.
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

// ===========================================================================
// SHARED LOGIC
// ===========================================================================

// Construct the most specialized IMap for the given parameters.
// Set computing_epsilon_deriv = true when computing derivatives w.r.t. epsilon to
// prevent selecting the epsilon_0 specialization, which assumes epsilon is fixed at 0.
IMap build_vector_field(int d, interval omega, interval sigma, interval delta, interval epsilon, bool computing_epsilon_deriv) {
    int dim = 6;
    IMap vf;

    bool d1_opt = (d == 1) && (omega == 1) && (delta == 0);
    bool d3_opt_eps0 = (d == 3) && (omega == 1) && (sigma == 1) && (delta == 0) && (epsilon == 0) && (!computing_epsilon_deriv);
    bool d3_opt = (d == 3) && (omega == 1) && (sigma == 1) && (delta == 0);

    if (d1_opt) {
	// Specialized for Case I in paper
        vf = IMap(vectorField_d1_optimized, dim, dim, 1);
        vf.setParameter(0, sigma);
    } else if (d == 1) {
	// Generic version for d = 1
        vf = IMap(vectorField_d1, dim, dim, 3);
        vf.setParameter(0, omega);
	vf.setParameter(1, sigma);
	vf.setParameter(2, delta);
    } else if (d3_opt_eps0) {
	// Specialized for Case II in paper when epsilon = 0
        vf = IMap(vectorField_d3_optimized_epsilon_0, dim, dim, 0);
    } else if (d3_opt) {
	// Specialized for Case II in paper
        vf = IMap(vectorField_d3_optimized, dim, dim, 0);
    } else {
	// Generic version
        vf = IMap(vectorField, dim, dim, 4);
        vf.setParameter(0, omega);
	vf.setParameter(1, sigma);
	vf.setParameter(2, delta);
	vf.setParameter(3, interval(d));
    }

    return vf;
}

// Integrate the CGL ODE from xi_0 to xi_1 and print the result to stdout.
// Outputs [a, b, alpha, beta] at xi_1 as intervals, one per line.
// Returns 0 on success, 1 on solver error (outputs NaN intervals).
int Q_zero(
    /** Initial value Q(xi_0) **/
    IVector Q_xi_0,
    /** Parameters **/
    int d,
    interval kappa,
    interval epsilon,
    interval omega,
    interval sigma,
    interval delta,
    /** Integration interval **/
    interval xi_0,
    interval xi_1,
    /** Settings **/
    double tol
) {
    // To get better enclosures for wide kappa and epsilon we treat
    // them as variables in the ODE.
    IVector Q_xi_0_with_parameters(6);
    for (int i = 0; i < 4; i++)
	Q_xi_0_with_parameters[i] = Q_xi_0[i];
    Q_xi_0_with_parameters[4] = kappa;
    Q_xi_0_with_parameters[5] = epsilon;

    // Choose optimized vector field based on parameters
    IMap vf = build_vector_field(d, omega, sigma, delta, epsilon, false);

    // Create the solver and the time map
    IOdeSolver solver(vf, 20);

    solver.setAbsoluteTolerance(tol);
    solver.setRelativeTolerance(tol);

    ITimeMap timeMap(solver);

    try {
	// Initial value for solver
	C0HORect2Set s(Q_xi_0_with_parameters, xi_0);

	// Solve the system
	IVector result = timeMap(xi_1, s);

	for (int i = 0; i < 4; i++) {
	    cout << result[i] << endl;
	}

	return 0; // Success
    } catch(...) {
	for (int i = 0; i < 4; i++) {
	    cout << interval(NAN) << endl;
	}

	return 1; // Error
    }
}

// Integrate the CGL ODE and compute its Jacobian from xi_0 to xi_1.
// Outputs 16 values (4x4 Jacobian w.r.t. initial conditions, column-major)
// followed by 4 values (derivative w.r.t. kappa or epsilon per wrt_epsilon).
// All values are intervals, one per line.
// Returns 0 on success, 1 on solver error (outputs 20 NaN intervals).
int Q_zero_jacobian(
    /** Initial value Q(xi_0) **/
    IVector Q_xi_0,
    /** Parameters **/
    int d,
    interval kappa,
    interval epsilon,
    interval omega,
    interval sigma,
    interval delta,
    /** Integration interval **/
    interval xi_0,
    interval xi_1,
    /** Settings **/
    bool wrt_epsilon,
    double tol
) {
    // To get better enclosures for wide kappa and epsilon we treat
    // them as variables in the ODE.
    IVector Q_xi_0_with_parameters(6);
    for (int i = 0; i < 4; i++)
	Q_xi_0_with_parameters[i] = Q_xi_0[i];
    Q_xi_0_with_parameters[4] = kappa;
    Q_xi_0_with_parameters[5] = epsilon;

    // Choose optimized vector field based on parameters
    IMap vf = build_vector_field(d, omega, sigma, delta, epsilon, wrt_epsilon);

    // Create the solver and the time map
    IOdeSolver solver(vf, 20);

    solver.setAbsoluteTolerance(tol);
    solver.setRelativeTolerance(tol);

    ITimeMap timeMap(solver);

    timeMap.stopAfterStep(true);

    try {
	// Define a representation of the initial value
	C1HORect2Set s(Q_xi_0_with_parameters, xi_0);

        // CAPD does not seem to support any way to limit the maximum
        // number of steps. This is a way to work limit the number of
        // steps by manually taking steps and aborting after too many
        // steps. The reason this is needed is that for too wide
        // input, it sometimes take an excessive number of steps only
        // to fail in the end anyway. Limiting the number of steps
        // avoids very large computational times in these extreme
        // cases.
        int steps = 0;
        int max_steps;
        // We take the maximum number of steps depending on xi_1.
        if (xi_1 < 50)
            max_steps = 1000;
        else if (xi_1 < 100)
            max_steps = 2000;
        else {
            max_steps = 10000;
        }
        do {
            timeMap(xi_1, s);
            steps++;
            if (steps >= max_steps) {
                for (int j = 0; j < 5; j++)
                    for (int i = 0; i < 4; i++) {
                        cout << interval(NAN) << endl;
                    }

                return 1; // Error
            }
        } while (!timeMap.completed());

	IMatrix m = (IMatrix)(s);

	for (int j = 0; j < 4; j++)
	    // Jacobian w.r.t. initial conditions
	    for (int i = 0; i < 4; i++) {
		cout << m[i][j] << endl;
	    }

	if (wrt_epsilon) {
	    // Derivative w.r.t. epsilon
	    for (int i = 0; i < 4; i++) {
		cout << m[i][5] << endl;
	    }
	} else {
	    // Derivative w.r.t. kappa
	    for (int i = 0; i < 4; i++) {
		cout << m[i][4] << endl;
	    }
	}

	return 0; // Success
    } catch(...) {
	for (int j = 0; j < 5; j++)
	    for (int i = 0; i < 4; i++) {
		cout << interval(NAN) << endl;
	    }

	return 1; // Error
    }
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
        cout << ";" << abs2_Q_derivative2 << endl; // Second derivative of abs(Q)^2
    } else {
        interval domain_left =
            interval(domain.leftBound(), (domain.leftBound() + domain.rightBound()) / 2);
        interval domain_right =
            interval((domain.leftBound() + domain.rightBound()) / 2, domain.rightBound());

        print_curve(curve, domain_left, prevTime, depth + 1);
        print_curve(curve, domain_right, prevTime, depth + 1);
    }
}

// Integrate the CGL ODE step-by-step and print an enclosure of the solution
// curve to stdout. Each row covers one sub-interval of [xi_0, xi_1] and
// contains: xi; a; b; alpha; beta; a''; b''; d|Q|^2/dxi; d^2|Q|^2/dxi^2,
// separated by semicolons. Sub-intervals where neither monotonicity witness
// is conclusive are bisected up to 5 levels; see print_curve for details.
// Returns 0 on success, 1 on solver error.
int Q_zero_curve(
    /** Initial value Q(xi_0) **/
    IVector Q_xi_0,
    /** Parameters **/
    int d,
    interval kappa,
    interval epsilon,
    interval omega,
    interval sigma,
    interval delta,
    /** Integration interval **/
    interval xi_0,
    interval xi_1,
    /** Settings **/
    double tol
) {
    // To get better enclosures for wide kappa and epsilon we treat
    // them as variables in the ODE.
    IVector Q_xi_0_with_parameters(6);
    for (int i = 0; i < 4; i++)
	Q_xi_0_with_parameters[i] = Q_xi_0[i];
    Q_xi_0_with_parameters[4] = kappa;
    Q_xi_0_with_parameters[5] = epsilon;

    // Choose optimized vector field based on parameters
    IMap vf = build_vector_field(d, omega, sigma, delta, epsilon, false);

    // Create the solver and the time map
    IOdeSolver solver(vf, 20);

    solver.setAbsoluteTolerance(tol);
    solver.setRelativeTolerance(tol);

    ITimeMap timeMap(solver);

    timeMap.stopAfterStep(true);

    interval prevTime(xi_0);

    try {
	// Initial value for solver
	C0HORect2Set s(Q_xi_0_with_parameters, xi_0);

	do {
	    timeMap(xi_1, s);

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

	return 0; // Success
    } catch(...) {
	// IMPROVE: Should we print something here?

	return 1; // Error
    }
}

int main()
{
  cout.precision(17); // Enough to exactly recover Float64 values
  cerr.precision(17); // Enough to exactly recover Float64 values

  // Read input from stdin

  // Read initial value
  IVector Q_xi_0(4);
  cin >> Q_xi_0[0] >> Q_xi_0[1] >> Q_xi_0[2] >> Q_xi_0[3];

  // Read parameter values
  int d;
  interval kappa, epsilon, omega, sigma, delta;
  cin >> d;
  cin >> kappa;
  cin >> epsilon;
  cin >> omega;
  cin >> sigma;
  cin >> delta;

  // Read time span
  interval xi_0, xi_1;
  cin >> xi_0 >> xi_1;

  // Read flags for if to output Jacobian and for if the Jacobian
  // should be with respect to epsilon instead of kappa.
  int output_jacobian;
  int wrt_epsilon;
  int output_curve;
  cin >> output_jacobian;
  cin >> wrt_epsilon;
  cin >> output_curve;

  // Read tolerance to use
  double tol;
  cin >> tol;

  if (output_jacobian) {
      return Q_zero_jacobian(Q_xi_0, d, kappa, epsilon, omega, sigma, delta, xi_0, xi_1, wrt_epsilon, tol);
  } else if (output_curve) {
      return Q_zero_curve(Q_xi_0, d, kappa, epsilon, omega, sigma, delta, xi_0, xi_1, tol);
  } else {
      return Q_zero(Q_xi_0, d, kappa, epsilon, omega, sigma, delta, xi_0, xi_1, tol);
  }
}
