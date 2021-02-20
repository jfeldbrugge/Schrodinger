// Parameters

const double ProperTime = 9.;
const double Leng = 20 + 0.00001;

const double mu = Leng / 2.;
const double sigma = 1;
const double p = 1.;

const double m = 1.;
const double hbar = 1.;

const double tau = 0.0001;
const double delta = 0.02;

const int L = int(Leng/delta) + 1;
int n = int(ProperTime/tau) + 1;

const int output = int( n / ( 8. * ProperTime ) );

// Initial Wave Function
void Initialize(vector1d &Phi, const double mu, const double sigma, const double p) {
#pragma omp parallel for
    for( int i = 0; i < L; i++) {
        Phi[i] =  1. / sqrt(sqrt(2 * pi * sigma * sigma))
        * exp(- (delta * i - mu) * (delta * i - mu) / (4 * sigma * sigma) + I * p * (delta * double(i) - mu) / hbar);
    }
}

// Potential
double V ( double x ) {
    return pow(x - Leng / 2., 2.);
}

// Suzuki-Trotter method to Fourth order in the space-time manifold and proper time.
const double a1 = hbar * hbar / (24. * m * delta * delta);
const double a2 = -2. * hbar * hbar / (3. * m * delta * delta);
double diag(double x) {return 5. * hbar * hbar / (4. * m * delta * delta) + V(x);}
