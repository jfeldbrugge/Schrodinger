//Parameters
const double hbar = 0.75;

const double muX = 25.;
const double muY = 5.;
const double sigmaX = 0.5;
const double sigmaY = 0.5;
const double pX = 0.;
const double pY = 5.;

const double muTr = 5.;
const double muXr = 5.;
const double sigmaTr = 0.5;
const double sigmaXr = 0.5;
const double pTr = 1.;
const double pXr = 0.;

const double m = 5.;

const double ProperTime = 30.;
const double LengthY = 50.;
const double LengthX = 50.;

const double tau = 0.1;
const double delta = 0.1;

const int LY = int(LengthY/delta) + 1;
const int LX = int(LengthX/delta) + 1;

int n = int(ProperTime/tau) + 1;

const int output = int( n / ( 2. * ProperTime ) );

// Initialize
void Initialize(vector2d &Phi, const double muY, const double muX, const double sigmaY, const double sigmaX, const double pY, const double pX)
{
#pragma omp parallel for
    for( int i = 0; i < LY; i++) {
        for( int j = 0; j < LX; j++) {
            Phi[i][j] =  1/sqrt(2 * pi * sigmaX * sigmaY )
            * exp(
                  -(delta * j - muX) * (delta * j - muX) / (4 * sigmaX * sigmaX)
                  -(delta * i - muY) * (delta * i - muY) / (4 * sigmaY * sigmaY)
                  - I * pY * (delta * double(i) - muY) / hbar
                  - I * pX * (delta * double(j) - muX) / hbar);
        }
    }
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

// The finite difference expansion of the Hamiltonian
inline std::complex<double> a1(int i1, int j1, int i2, int j2){return -2. * hbar * hbar / (3.  * delta * delta * m);}
inline std::complex<double> a2(int i1, int j1, int i2, int j2){return 1. * hbar * hbar / (24. * delta * delta * m);}
inline std::complex<double> diag(int i, int j){return 5. * hbar * hbar/ (2. * delta * delta * m) + V(delta * j, delta * i);}
