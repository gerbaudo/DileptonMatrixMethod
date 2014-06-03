# include "SameSignMatrixMethod/DiLeptonMatrixMethod.h"

#include <iostream>
#include <string>
#include <cmath>

/*
  Test the matrix method with the recommended, most up-to-date matrix.

  davide.gerbaudo@gmail.com
  April 2014
*/

using namespace std;
namespace smm = SameSignMatrixMethod;
namespace sf = susy::fake;

//----------------------------------------------------------
bool testParametrization(string filename, smm::RATE_PARAM rp)
{
    bool success = false;
    smm::DiLeptonMatrixMethod matrix;
    if(matrix.configure(filename, rp, rp, rp, rp)) {
        float gev2mev(1.0e3);
        bool l0IsSig(true), l0IsEle(false);
        bool l1IsSig(false), l1IsEle(true);
        float l0Pt(30.0), l0Eta(+0.5);
        float l1Pt(25.0), l1Eta(-0.5);
        float metRel(20.0);
        sf::Region region = sf::CR_SSInc1j;
        smm::SYSTEMATIC sys = smm::SYS_NOM;
        float weight = matrix.getTotalFake(l0IsSig, l0IsEle, l0Pt*gev2mev, l0Eta,
                                           l1IsSig, l1IsEle, l1Pt*gev2mev, l1Eta,
                                           region, metRel*gev2mev, sys);
        float weightUp = matrix.getTotalFake(l0IsSig, l0IsEle, l0Pt*gev2mev, l0Eta,
                                             l1IsSig, l1IsEle, l1Pt*gev2mev, l1Eta,
                                             region, metRel*gev2mev, smm::SYS_EL_FR_UP);
        float weightDo = matrix.getTotalFake(l0IsSig, l0IsEle, l0Pt*gev2mev, l0Eta,
                                             l1IsSig, l1IsEle, l1Pt*gev2mev, l1Eta,
                                             region, metRel*gev2mev, smm::SYS_EL_FR_DOWN);
        cout<<"weight for"
            <<" "<<sf::region2str(region)<<", "
            <<" "<<smm::systematic_names[sys]
            <<" : "<<weight
            <<" (SYS_EL_FR_UP "<<weightUp
            <<", "
            <<" SYS_EL_FR_DOWN "<<weightDo
            <<")"
            <<endl;
        success = true;
    }
    return success;
}
//----------------------------------------------------------
int main(int argc, char **argv)
{

  string inputFilename="data/FinalFakeHist_Apr_10.root";
  inputFilename = "data/FinalFakeHist_May_20.root";
  size_t nFail=0;

  cout<<endl<<" --- test 1D parametrization (pt)     ---"<<endl;
  if(!testParametrization(inputFilename, smm::PT))
      nFail++;
  cout<<endl<<" --- test 2D parametrization (pt eta) ---"<<endl;
  if(!testParametrization(inputFilename, smm::PT_ETA))
      nFail++;
  cout<<endl<<" --- Number of failures: "<<nFail<<" --- "<<endl;
  return 0;
}
