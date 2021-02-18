// Integrate
void Integrate(const vector2d &Phi, vector2d &PhiT)
{
#pragma omp parallel for
    for( int ii = 0; ii < LY; ii++ )
    {
        for( int jj = 0; jj < LX; jj++)
        {
            PhiT[ii][jj] += Phi[ii][jj];
        }
    }
}

// Check normalization
double Norm(const vector2d &Phi)
{
    double sum = 0;
    for( int ii = 0; ii < LY; ii++)
    {
        for( int jj = 0; jj < LX; jj++)
        {
            sum = sum + abs(Phi[ii][jj]) * abs(Phi[ii][jj]);
        }
    }
    return sum * delta * delta;
}

// Output
void WriteSlice(const vector2d &Phi, std::ofstream &writer)
{
    for(int ii = 0; ii < LY; ii++)
    {
        for(int jj = 0; jj < LX; jj++)
        {
            double gg = real(Phi[ii][jj]), hh = imag(Phi[ii][jj]);
            writer.write((char*) &gg, sizeof(double));
            writer.write((char*) &hh, sizeof(double));
        }
    }
}
