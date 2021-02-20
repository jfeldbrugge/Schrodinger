//Parameters
const double hbar = 0.75;

const double muX = 25.;
const double muY = 5.;
const double sigmaX = 0.5;
const double sigmaY = 0.5;
const double pX = 0.;
const double pY = 5.;

const double m = 5.;
const double e = 0.;

const double ProperTime = 20.;
const double LengthY = 50.;
const double LengthX = 50.;

const double tau = 0.025;
const double delta = 0.05;

const int LY = int(LengthY/delta) + 1;
const int LX = int(LengthX/delta) + 1;

int n = int(ProperTime/tau) + 1;

const int output = int( n / ( 2. * ProperTime ) );

// Initialize
void Initialize(vector2d &Phi, const double muY, const double muX, const double sigmaY, const double sigmaX, const double pY, const double pX) {
#pragma omp parallel for
    for( int i = 0; i < LY; i++) {
        for( int j = 0; j < LX; j++) {
            Phi[i][j] =  1/sqrt(2 * pi * sigmaX * sigmaY )
            * exp(-(delta * j - muX) * (delta * j - muX) / (4 * sigmaX * sigmaX)
                  -(delta * i - muY) * (delta * i - muY) / (4 * sigmaY * sigmaY)
                  - I * (pY - e * Ax(i,j)) * (delta * double(i) - muY) / hbar
                  - I * (pX - e * Ay(i,j)) * (delta * double(j) - muX) / hbar);
        }
    }
}

// Vector potential
double Ax(const int ii, const int jj) {
    return 0.;
}

double Ay(const int ii, const int jj) {
    return 0.;
}

// Potential
double V(const double x, const double y) {
    double pot = 0.;
    const double yscreen = 25.;
    const double screenwidth = 1.;
    const double xslit1 = 23.;
    const double xslit2 = 27.;
    const double slitwidth = 1.;
    
    if(y > yscreen - screenwidth / 2. && y < yscreen + screenwidth / 2.) {
        if((x < xslit1 - slitwidth / 2. || (x > xslit1 + slitwidth / 2. && x < xslit2 - slitwidth / 2.) || x > xslit2 + slitwidth / 2.)) {
            pot = 100.;
        }
    }
    
    return pot;
}

inline std::complex<double> a1x(int i1, int j1, int i2, int j2) {
    return +1. / (24. * m) * ( hbar * hbar / (delta * delta) + I * e * hbar / delta * ( Ax(i1,j1) + Ax(i2,j2)));
}

inline std::complex<double> a1y(int i1, int j1, int i2, int j2) {
    return +1. / (24. * m) * ( hbar * hbar / (delta * delta) + I * e * hbar / delta * ( Ay(i1,j1) + Ay(i2,j2)));
}

inline std::complex<double> a2x(int i1, int j1, int i2, int j2) {
    return -1. / (3. * m) * ( 2. * hbar * hbar / (delta * delta) + I * e * hbar / delta * ( Ax(i1,j1) + Ax(i2,j2)));
}

inline std::complex<double> a2y(int i1, int j1, int i2, int j2) {
    return -1. / (3. * m) * ( 2. * hbar * hbar / (delta * delta) + I * e * hbar / delta * ( Ay(i1,j1) + Ay(i2,j2)));
}

inline std::complex<double> diag(int i1, int j1) {
    return 1. / (2. * m) * ( 5. * hbar * hbar / (delta * delta) + e * e * (Ax(i1,j1) * Ax(i1,j1) + Ay(i1,j1) * Ay(i1,j1))) + V(delta * j1, delta * i1);
}
