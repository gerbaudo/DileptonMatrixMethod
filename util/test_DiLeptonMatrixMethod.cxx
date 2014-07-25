#include "DileptonMatrixMethod/DileptonMatrixMethod.h"

#include <iostream>
#include <string>
#include <cmath>

/*
  Test the matrix method with the recommended, most up-to-date matrix.

  davide.gerbaudo@gmail.com
  April 2014
*/

using namespace std;
namespace sf = susy::fake;
using sf::Parametrization;
using sf::Systematic;
using sf::Lepton;

//----------------------------------------------------------
bool testParametrization(string filename, Parametrization::Value rp)
{
    bool success = false;
    sf::DileptonMatrixMethod matrix;
    std::vector<std::string> regions;
    const std::string regionName = "CR_SSInc1j";
    regions.push_back(regionName);
    if(matrix.configure(filename, regions, rp, rp, rp, rp)) {
        float gev2mev(1.0e3);
        bool l0IsSig(true), l0IsEle(false);
        bool l1IsSig(false), l1IsEle(true);
        float l0Pt(30.0), l0Eta(+0.5);
        float l1Pt(25.0), l1Eta(-0.5);
        float metRel(20.0);
        size_t iRegion = matrix.getIndexRegion(regionName);
        Systematic::Value sys = Systematic::SYS_NOM;
        Systematic::Value sysUp = Systematic::SYS_EL_FR_UP;
        Systematic::Value sysDo = Systematic::SYS_EL_FR_DOWN;
        Lepton l0(l0IsSig, l0IsEle, l0Pt*gev2mev, l0Eta);
        Lepton l1(l1IsSig, l1IsEle, l1Pt*gev2mev, l1Eta);
        float weight   = matrix.getTotalFake(l0, l1, iRegion, metRel*gev2mev, sys);
        float weightUp = matrix.getTotalFake(l0, l1, iRegion, metRel*gev2mev, sysUp);
        float weightDo = matrix.getTotalFake(l0, l1, iRegion, metRel*gev2mev, sysDo);
        cout<<"weight for"
            <<" "<<regionName<<", "
            <<" "<<Systematic::str(sys)
            <<" : "<<weight
            <<" ("<<Systematic::str(sysUp)<<" "<<weightUp
            <<", "
            <<" "<<Systematic::str(sysDo)<<" "<<weightDo
            <<")"
            <<endl;
        success = true;
    } else {
        cout<<"Failed to configure, cannot test"<<endl;
    }
    return success;
}
//----------------------------------------------------------
int main(int argc, char **argv)
{

  string inputFilename="data/FinalFakeHist_Apr_10.root";
  inputFilename = "data/FinalFakeHist_May_16.root";
  size_t nFail=0;

  cout<<endl<<" --- test 1D parametrization (pt)     ---"<<endl;
  if(!testParametrization(inputFilename, Parametrization::PT))
      nFail++;
  cout<<endl<<" --- test 2D parametrization (pt eta) ---"<<endl;
  if(!testParametrization(inputFilename, Parametrization::PT_ETA))
      nFail++;
  cout<<endl<<" --- Number of failures: "<<nFail<<" --- "<<endl;
  return 0;
}
