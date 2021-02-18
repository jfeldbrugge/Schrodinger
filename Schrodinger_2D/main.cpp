#include "Schrodinger.Header.h"
#include "Schrodinger.Physics.h"
#include "Schrodinger.SuzukiTrotterFourthOrder.h"
#include "Schrodinger.Utilities.h"

//Main Program
int main(int argc, char* argv[])
{
    //Arguments
    if(argc == 2) {n = atoi(argv[1]);} else { }
    
    //Initialize
    std::cout << "Schrodinger solver" << std::endl;
    
    std::ofstream myfileB; myfileB.open("output.bin", std::ios::binary);
    
    vector2d Phi(LY, std::vector<std::complex<double> > (LX));
    
    Initialize(Phi, muY, muX, sigmaY, sigmaX, pY, pX);
    WriteSlice(Phi, myfileB);
    std::cout << "t = " << 0 << " |Phi|^2 = " << Norm(Phi) << std::endl;
    
    //Integrate, forward in time
    for(int iteration = 1; iteration < n; iteration++) {
        U4(Phi, tau / hbar);
        if(iteration % output == 0) {
            std::cout << "t = " << tau * iteration << " |Phi|^2 = " << Norm(Phi) << std::endl;
            WriteSlice(Phi, myfileB);
        }
    }
    myfileB.close();
    return 0;
}
