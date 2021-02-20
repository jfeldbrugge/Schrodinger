#include <iostream>
#include <fstream>
//#include <omp.h>
#include <vector>
#include <cmath>
#include <complex>
#include <cstdlib>

const double pi = 3.14159265359;
const std::complex<double> I (0,1.0);

//Data Types
typedef std::vector<std::vector<std::complex<double> > > vector2d;

//Prototypes
void Initialize(vector2d &Phi, const double muT, const double muX, const double sigmaT, const double sigmaX, const double pT, const double pX);

void Exp(vector2d &Phi, const double tau, const std::complex<double> a, const int i1, const int j1, const int i2, const int j2);

void ExpA1(vector2d &Phi, const double tau);  // T Derivative
void ExpA2(vector2d &Phi, const double tau);  // T Derivative
void ExpA3(vector2d &Phi, const double tau);  // T Derivative
void ExpA4(vector2d &Phi, const double tau);  // T Derivative
void ExpA5(vector2d &Phi, const double tau);  // X-Derivative
void ExpA6(vector2d &Phi, const double tau);  // X-Derivative
void ExpA7(vector2d &Phi, const double tau);  // X-Derivative
void ExpA8(vector2d &Phi, const double tau);  // X-Derivative
void ExpA9(vector2d &Phi, const double tau);  // Diagonal term

void U1(vector2d &Phi, const double tau);
void U2(vector2d &Phi, const double tau);
void U4(vector2d &Phi, const double tau);

double Ax(const int ii, const int jj);
double Ay(const int ii, const int jj);
double V(const double x, const double y);

void Integrate(const vector2d &Phi, vector2d &PhiT);

double Norm(const vector2d &Phi);

void WriteSlice(const vector2d &Phi, std::ofstream &writer);
