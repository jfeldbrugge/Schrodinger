// Implementation
void U1(vector1d &Phi, const double tau)
{
    ExpA1(Phi, tau);
    ExpA2(Phi, tau);
    ExpA3(Phi, tau);
    ExpA4(Phi, tau);
    ExpA5(Phi, tau);
}

void U2(vector1d &Phi, const double tau)
{
    ExpA1(Phi, tau / 2.);
    ExpA2(Phi, tau / 2.);
    ExpA3(Phi, tau / 2.);
    ExpA4(Phi, tau / 2.);
    ExpA5(Phi, tau);
    ExpA4(Phi, tau / 2);
    ExpA3(Phi, tau / 2);
    ExpA2(Phi, tau / 2.);
    ExpA1(Phi, tau / 2.);
}

void U4(vector1d &Phi, const double tau)
{
    const double p = 1. / (4. - cbrt(4.));
    U2(Phi, p * tau);
    U2(Phi, p * tau);
    U2(Phi, (1. - 4. * p ) * tau);
    U2(Phi, p * tau);
    U2(Phi, p * tau);
}

// Exponentiation
void ExpA1(vector1d &Phi, const double tau)
{
#pragma omp parallel for
    for(int i1 = 0; i1 < L - 2; i1 = i1 + 4)
    {
        Exp(Phi, tau / hbar, -a1, i1, i1 + 2);
    }
#pragma omp parallel for
    for(int i1 = 1; i1 < L - 2; i1 = i1 + 4)
    {
        Exp(Phi, tau / hbar, -a1, i1, i1 + 2);
    }
}

void ExpA2(vector1d &Phi, const double tau)
{
#pragma omp parallel for
    for(int i1 = 2; i1 < L - 2; i1 = i1 + 4)
    {
        Exp(Phi, tau / hbar, -a1, i1, i1 + 2);
    }
#pragma omp parallel for
    for(int i1 = 3; i1 < L - 2; i1 = i1 + 4)
    {
        Exp(Phi, tau / hbar, -a1, i1, i1 + 2);
    }
}

void ExpA3(vector1d &Phi, const double tau)
{
#pragma omp parallel for
    for(int i1 = 0; i1 < L - 1; i1 = i1 + 2)
    {
        Exp(Phi, tau / hbar, -a2, i1, i1 + 1);
    }
}

void ExpA4(vector1d &Phi, const double tau)
{
#pragma omp parallel for
    for(int i1 = 1; i1 < L - 1; i1 = i1 + 2)
    {
        Exp(Phi, tau / hbar, -a2, i1, i1 + 1);
    }
}

void ExpA5(vector1d &Phi, const double tau)
{
#pragma omp parallel for
    for( int i = 0; i < L; i++)
    {
        Phi[i] = exp(-I * tau / hbar * diag(delta * i)) * Phi[i];
    }
}

void Exp(vector1d &Phi, const double tau, const std::complex<double> a, const int i1, const int i2)
{
    std::complex<double> temp = Phi[i1];
    const double norm = std::abs(a);
    
    Phi[i1] = cos(tau * norm) * Phi[i1]   + I *      a  / norm * sin(tau * norm) * Phi[i2];
    Phi[i2] = cos(tau * norm) * Phi[i2]   + I * std::conj(a) / norm * sin(tau * norm) * temp;
}
