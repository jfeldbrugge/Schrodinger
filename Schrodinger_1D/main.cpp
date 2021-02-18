#include "Schrodinger.Header.h"
#include "Schrodinger.Physics.h"
//#include "Schrodinger.SuzukiTrotterFirstOrder.h"
#include "Schrodinger.SuzukiTrotterFourthOrder.h"
#include "Schrodinger.Utilities.h"

//Main Program
int main(int argc, char* argv[]) {
    //Arguments
    if(argc == 2) {n = atoi(argv[1]);} else { }
    
    //Initialize
    std::cout << "Schrodinger solver"<< std::endl;
    
    std::ofstream myfileB; myfileB.open("output.bin", std::ios::binary);
    
    vector1d Phi(L);
    Initialize(Phi, mu, sigma, p);
    WriteSlice(Phi, myfileB);
    std::cout << "n = " << 0 << " |Phi|^2 = " << Norm(Phi) << std::endl;
    
    //Integrate
    for(int iteration = 1; iteration < n; iteration++) {
        U4(Phi, tau);
        if(iteration % output == 0) {
            std::cout << "t = " << tau * iteration << " |Phi|^2 = " << Norm(Phi) << std::endl;
            WriteSlice(Phi, myfileB);
        }
    }
    myfileB.close();
    
    return 0;
}
