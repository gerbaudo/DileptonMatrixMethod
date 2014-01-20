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
namespace ssmm = SameSignMatrixMethod;
namespace sf = susy::fake;

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

    struct MmPars {
        bool l0IsSig, l1IsSig, l0IsEle, l1IsEle;
        float l0Pt, l1Pt, l0Eta, l1Eta, metRel;
        susy::fake::Region region;
        ssmm::SYSTEMATIC sys;
        MmPars(bool s0, bool e0, float pt0, float eta0,
               bool s1, bool e1, float pt1, float eta1,
               float mrel,
               sf::Region r, ssmm::SYSTEMATIC s) :
            l0IsSig(s0), l1IsSig(s1), l0IsEle(e0), l1IsEle(e1),
            l0Pt(pt0), l1Pt(pt1), l0Eta(eta0), l1Eta(eta1), metRel(mrel),
            region(r), sys(s) {};
    };

int main(int argc, char** argv)
{
    ROOT::Cintex::Cintex::Enable();
    MmPars ps[] = {
        MmPars(true, true, 20.0, 1.2, true, true, 20.0, 1.2, 40.0, sf::CR_CR8ee, ssmm::SYS_NONE)
    };
    float gev2mev=  1.0e3;
    string matrixFilename = getRootCoreDir()+"/../SameSignMatrixMethod/data/FinalFakeHist_Nov_10.root";
    ssmm::DiLeptonMatrixMethod dmm;
    ssmm::RATE_PARAM parVar = ssmm::PT; // parametrization variable
    dmm.configure(matrixFilename, parVar, parVar, parVar, parVar);
    /*const*/ MmPars &p = ps[0];
    dmm.getTotalFake(p.l0IsSig, p.l0IsEle, p.l0Pt*gev2mev, p.l0Eta,
                     p.l1IsSig, p.l1IsEle, p.l1Pt*gev2mev, p.l1Eta,
                     p.region,
                     p.metRel*gev2mev,
                     p.sys);



  return 0;
}
