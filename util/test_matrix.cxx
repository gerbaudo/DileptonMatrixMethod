#include <cassert>
#include <cstdlib>
#include <iostream>
#include <string>

#include "Cintex/Cintex.h"

#include "SameSignMatrixMethod/DiLeptonMatrixMethod.h"

/*

*/


using std::cout;
using std::endl;
using std::string;


std::string getRootCoreDir()
{
    using namespace std;
    string dir;
    char* rootcoredir = getenv("ROOTCOREDIR");
    bool envvarDefined(rootcoredir!=0);
    if (envvarDefined) { dir = rootcoredir; }
    else               { cout<<"getRootCoreDir: ROOTCOREDIR not defined"<<endl; }
    return dir;
}

int main(int argc, char** argv)
{
    ROOT::Cintex::Cintex::Enable();
    string matrixFilename = getRootCoreDir()+"/../SameSignMatrixMethod/data/FinalFakeHist_Nov_10.root";
    SameSignMatrixMethod::DiLeptonMatrixMethod dmm;
    dmm.configure(matrixFilename,
                  SameSignMatrixMethod::PT,     // Electron Real
                  SameSignMatrixMethod::PT,     // Electron Fake
                  SameSignMatrixMethod::PT,     // Muon Real
                  SameSignMatrixMethod::PT      // Muon Fake
        );



  return 0;
}
