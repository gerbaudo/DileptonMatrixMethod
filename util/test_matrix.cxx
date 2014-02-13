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
    float expected;
    MmPars(bool s0, bool e0, float pt0, float eta0,
           bool s1, bool e1, float pt1, float eta1,
           float mrel,
           sf::Region r, ssmm::SYSTEMATIC s,
           float exp) :
        l0IsSig(s0), l1IsSig(s1), l0IsEle(e0), l1IsEle(e1),
        l0Pt(pt0), l1Pt(pt1), l0Eta(eta0), l1Eta(eta1), metRel(mrel),
        region(r), sys(s),
        expected(exp) {};
};

bool equal(float v1, float v2, float precision=1e-6)
{
    float pDiff((v1+v2) ? fabs(v1-v2)/(v1+v2) : (v1+v2));
    return pDiff < precision;
}

int main(int argc, char** argv)
{
    ROOT::Cintex::Cintex::Enable();
    sf::Region reg = sf::CR_CR8ee;
    ssmm::SYSTEMATIC sys = ssmm::SYS_NONE;
    MmPars ps[] = {
        MmPars(true,  true,  32.0, 1.2,  true,  true, 23.0, 1.2, 40.0, reg, sys, 0.5), // ee
        MmPars(true,  false, 32.0, 1.2,  true,  true, 23.0, 1.2, 40.0, reg, sys, 0.5),
        MmPars(false, true,  32.0, 1.2,  true,  true, 23.0, 1.2, 40.0, reg, sys, 0.5),
        MmPars(false, false, 32.0, 1.2,  true,  true, 23.0, 1.2, 40.0, reg, sys, 0.5),
        MmPars(true,  true,  32.0, 1.2,  true, false, 23.0, 1.2, 40.0, reg, sys, 0.5), // em
        MmPars(true,  false, 32.0, 1.2,  true, false, 23.0, 1.2, 40.0, reg, sys, 0.5),
        MmPars(false, true,  32.0, 1.2,  true, false, 23.0, 1.2, 40.0, reg, sys, 0.5),
        MmPars(false, false, 32.0, 1.2,  true, false, 23.0, 1.2, 40.0, reg, sys, 0.5),
        MmPars(true,  true,  32.0, 1.2, false,  true, 23.0, 1.2, 40.0, reg, sys, 0.5), // me
        MmPars(true,  false, 32.0, 1.2, false,  true, 23.0, 1.2, 40.0, reg, sys, 0.5),
        MmPars(false, true,  32.0, 1.2, false,  true, 23.0, 1.2, 40.0, reg, sys, 0.5),
        MmPars(false, false, 32.0, 1.2, false,  true, 23.0, 1.2, 40.0, reg, sys, 0.5),
        MmPars(true,  true,  32.0, 1.2, false, false, 23.0, 1.2, 40.0, reg, sys, 0.5), // mm
        MmPars(true,  false, 32.0, 1.2, false, false, 23.0, 1.2, 40.0, reg, sys, 0.5),
        MmPars(false, true,  32.0, 1.2, false, false, 23.0, 1.2, 40.0, reg, sys, 0.5),
        MmPars(false, false, 32.0, 1.2, false, false, 23.0, 1.2, 40.0, reg, sys, 0.5)
    };
    float gev2mev=  1.0e3;
    string matrixFilename = getRootCoreDir()+"/../SameSignMatrixMethod/data/FinalFakeHist_Nov_10.root";
    ssmm::DiLeptonMatrixMethod dmm;
    ssmm::RATE_PARAM parVar = ssmm::PT; // parametrization variable
    dmm.configure(matrixFilename, parVar, parVar, parVar, parVar);
    const size_t nPars = sizeof(ps)/sizeof(ps[0]);
    size_t nFail = 0;
    for(size_t i=0; i<nPars; ++i) {
        /*const*/ MmPars &p = ps[i];
        float weight(dmm.getTotalFake(p.l0IsSig, p.l0IsEle, p.l0Pt*gev2mev, p.l0Eta,
                                      p.l1IsSig, p.l1IsEle, p.l1Pt*gev2mev, p.l1Eta,
                                      p.region,
                                      p.metRel*gev2mev,
                                      p.sys));
        if(!equal(p.expected, weight)) {
            cout<<"["<<i<<"] : expected "<<p.expected<<" != "<<weight<<endl;
            nFail++;
        }
    }
    return 0;
}
