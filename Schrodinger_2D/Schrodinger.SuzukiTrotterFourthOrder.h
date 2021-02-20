// Implementation of Suzuki-Trotter method fo fourth order
void U1(vector2d &Phi, const double tau) {
    ExpA1(Phi, tau);
    ExpA2(Phi, tau);
    ExpA3(Phi, tau);
    ExpA4(Phi, tau);
    ExpA5(Phi, tau);
    ExpA6(Phi, tau);
    ExpA7(Phi, tau);
    ExpA8(Phi, tau);
    ExpA9(Phi, tau);
}

void U2(vector2d &Phi, const double tau) {
    ExpA1(Phi, tau / 2.);
    ExpA2(Phi, tau / 2.);
    ExpA3(Phi, tau / 2.);
    ExpA4(Phi, tau / 2.);
    ExpA5(Phi, tau / 2.);
    ExpA6(Phi, tau / 2.);
    ExpA7(Phi, tau / 2.);
    ExpA8(Phi, tau / 2.);
    ExpA9(Phi, tau);
    ExpA8(Phi, tau / 2.);
    ExpA7(Phi, tau / 2.);
    ExpA6(Phi, tau / 2.);
    ExpA5(Phi, tau / 2.);
    ExpA4(Phi, tau / 2.);
    ExpA3(Phi, tau / 2.);
    ExpA2(Phi, tau / 2.);
    ExpA1(Phi, tau / 2.);
}

void U4(vector2d &Phi, const double tau) {
    const double p = 1. / (4. - cbrt(4.));
    U2(Phi, p * tau);
    U2(Phi, p * tau);
    U2(Phi, (1. - 4. * p ) * tau);
    U2(Phi, p * tau);
    U2(Phi, p * tau);
}

// The partition of the Hamiltonian in block diagonal form
void ExpA1(vector2d &Phi, const double tau) {
#pragma omp parallel for
    for( int j1 = 0; j1 < LX; j1++)
    {
        for(int i1 = 0; i1 < LY - 2; i1 = i1 + 4){Exp(Phi, tau / hbar, a1x(i1, j1, i1 + 2, j1), i1, j1, i1 + 2, j1);}
        for(int i1 = 1; i1 < LY - 2; i1 = i1 + 4){Exp(Phi, tau / hbar, a1x(i1, j1, i1 + 2, j1), i1, j1, i1 + 2, j1);}
    }
}

void ExpA2(vector2d &Phi, const double tau) {
#pragma omp parallel for
    for( int j1 = 0; j1 < LX; j1++)
    {
        for(int i1 = 2; i1 < LY - 2; i1 = i1 + 4){Exp(Phi, tau / hbar, a1x(i1, j1, i1 + 2, j1), i1, j1, i1 + 2, j1);}
        for(int i1 = 3; i1 < LY - 2; i1 = i1 + 4){Exp(Phi, tau / hbar, a1x(i1, j1, i1 + 2, j1), i1, j1, i1 + 2, j1);}
    }
}

void ExpA3(vector2d &Phi, const double tau) {
#pragma omp parallel for
    for(int j1 = 0; j1 < LX; j1++)
    {
        for(int i1 = 0; i1 < LY - 1; i1 = i1 + 2){Exp(Phi, tau / hbar, a2x(i1, j1, i1 + 1, j1), i1, j1, i1 + 1, j1);}
    }
}

void ExpA4(vector2d &Phi, const double tau) {
#pragma omp parallel for
    for(int j1 = 0; j1 < LX; j1++)
    {
        for(int i1 = 1; i1 < LY - 1; i1 = i1 + 2){Exp(Phi, tau / hbar, a2x(i1, j1, i1 + 1, j1), i1, j1, i1 + 1, j1);}
    }
}

void ExpA5(vector2d &Phi, const double tau) {
#pragma omp parallel for
    for(int i1 = 0; i1 < LY; i1++)
    {
        for(int j1 = 0; j1 < LX - 2; j1 = j1 + 4){Exp(Phi, tau / hbar, a1y(i1, j1, i1, j1 + 2), i1, j1, i1, j1 + 2);}
        for(int j1 = 1; j1 < LX - 2; j1 = j1 + 4){Exp(Phi, tau / hbar, a1y(i1, j1, i1, j1 + 2), i1, j1, i1, j1 + 2);}
    }
}

void ExpA6(vector2d &Phi, const double tau) {
#pragma omp parallel for
    for(int i1 = 0; i1 < LY; i1++)
    {
        for(int j1 = 2; j1 < LX - 2; j1 = j1 + 4){Exp(Phi, tau / hbar, a1y(i1, j1, i1, j1 + 2), i1, j1, i1, j1 + 2);}
        for(int j1 = 3; j1 < LX - 2; j1 = j1 + 4){Exp(Phi, tau / hbar, a1y(i1, j1, i1, j1 + 2), i1, j1, i1, j1 + 2);}
    }
}

void ExpA7(vector2d &Phi, const double tau) {
#pragma omp parallel for
    for(int i1 = 0; i1 < LY; i1++)
    {
        for( int j1 = 0; j1 < LX - 1; j1 = j1 + 2){Exp(Phi, tau / hbar, a2y(i1, j1, i1, j1 + 1), i1, j1, i1, j1 + 1);}
    }
}

void ExpA8(vector2d &Phi, const double tau) {
#pragma omp parallel for
    for(int i1 = 0; i1 < LY; i1++)
    {
        for( int j1 = 1; j1 < LX - 1; j1 = j1 + 2){Exp(Phi, tau / hbar, a2y(i1, j1, i1, j1 + 1), i1, j1, i1, j1 + 1);}
    }
}

void ExpA9(vector2d &Phi, const double tau) {
#pragma omp parallel for
    for( int i1 = 0; i1 < LY; i1++)
    {
        for(int j1 = 0; j1 < LX; j1++){Phi[i1][j1] = exp( -I * tau / hbar * diag(i1,j1)) * Phi[i1][j1];}
    }
}

// Exponentiation
void Exp(vector2d &Phi, const double tau, const std::complex<double> a, const int i1, const int j1, const int i2, const int j2) {
    std::complex<double> temp = Phi[i1][j1];
    const double norm = abs(a);
    
    Phi[i1][j1] = cos(tau * norm) * Phi[i1][j1]   + I *      a  / norm * sin(tau * norm) * Phi[i2][j2];
    Phi[i2][j2] = cos(tau * norm) * Phi[i2][j2]   + I * conj(a) / norm * sin(tau * norm) * temp;
}
