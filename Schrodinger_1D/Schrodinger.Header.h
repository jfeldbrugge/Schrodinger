#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <complex>
#include <cstdlib>
#if defined(_OPENMP)
#include <omp.h>
extern const bool parallelism_enabled = true;
#else
extern const bool parallelism_enabled = false;
#endif

const double pi = 3.14159265359;
const std::complex<double> I (0,1.0);

//Data Types
typedef std::vector<std::complex<double> > vector1d;

//Prototypes
void Initialize(vector1d &Phi, const double mu, const double sigma, const double p);
double V(double x);

void Exp(vector1d &Phi, const double tau, const std::complex<double> a, const int i1, const int i2);

void ExpA1(vector1d &Phi, const double tau);
void ExpA2(vector1d &Phi, const double tau);
void ExpA3(vector1d &Phi, const double tau);
void ExpA4(vector1d &Phi, const double tau);
void ExpA5(vector1d &Phi, const double tau);

void U1(vector1d &Phi, const double tau);
void U2(vector1d &Phi, const double tau);
void U4(vector1d &Phi, const double tau);

void Integrate(const vector1d &Phi, vector1d &PhiT);
double Norm(const vector1d &Phi);

void WriteSlice(const vector1d &Phi, std::ofstream &writer);
void ReadSlice(vector1d &Phi, std::ifstream &readR, std::ifstream &readI);
