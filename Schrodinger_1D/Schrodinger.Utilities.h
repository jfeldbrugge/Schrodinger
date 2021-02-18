// Integrate
void Integrate ( const vector1d &Phi, vector1d &PhiT )
{
#pragma omp parallel for
    for( int i = 0; i < L; i++ )
    {
        PhiT[i] += Phi[i];
    }
}

// Compute Norm
double Norm(const vector1d &Phi)
{
    double sum = 0;
    for( int i = 0; i < L; i++)
    {
        sum = sum + abs(Phi[i]) * abs(Phi[i]);
    }
    return sum * delta;
}

// Output
void WriteSlice(const vector1d &Phi, std::ofstream &writer)
{
    for(int i = 0; i < L; i++)
    {
        double gg = real(Phi[i]), hh = imag(Phi[i]);
        writer.write((char*) &gg, sizeof(double));
        writer.write((char*) &hh, sizeof(double));
    }
}

// Input
void ReadSlice(vector1d &Phi, std::ifstream &readR, std::ifstream &readI)
{
    if(readR.is_open() && readI.is_open())
    {
        for(int i = 0; i < L; i++)
        {
            double gg, hh;
            readR.read((char*) &gg, sizeof(double));
            readI.read((char*) &hh, sizeof(double));
            Phi[i] = gg + I * hh;
        }
        
        std::cout << Phi[0] << std::endl;
        std::cout << Phi[1] << std::endl;
    } else
    {
        std::cout << "Error: could not open initial wave function files." << std::endl;
    }
}
