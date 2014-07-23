// -*- c++ -*-
#ifndef SUSY_FAKE_STATERRORTOOL_H
#define SUSY_FAKE_STATERRORTOOL_H

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//
// This tool aims to handle the cases where we are lacking sufficient  //
// pairs in a given region.  It will take as input the number of events//
// for various pairs, the current error, and fake region               //
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//


#include "DileptonMatrixMethod/DiLeptonMatrixMethod.h"
#include "DileptonMatrixMethod/FakeRegions.h"
#include "TRandom3.h"

namespace susy{
namespace fake{

  enum DileptonType {DT_ee = 0, DT_mm, DT_em, DT_N};

  class StatErrorTool : public susy::fake::DiLeptonMatrixMethod
  {
    
  public:
    StatErrorTool(std::string fakeFile, 
		  susy::fake::RATE_PARAM rate_param_real_el,
		  susy::fake::RATE_PARAM rate_param_fake_el,
		  susy::fake::RATE_PARAM rate_param_real_mu,
		  susy::fake::RATE_PARAM rate_param_fake_mu
		  );
    ~StatErrorTool();
    
    // Calculate the stat error based on Region
    // and the composition of the events.
    float getUpdatedError(float currentError,                    // current statistical error
			  int n_tt,                              // # of TT events
			  int n_tl,                              // # of TL events
			  int n_lt,                              // # of LT events
			  int n_ll,                              // # of LL events
              susy::fake::Region fr,       // Fake Region to get average rates
			  DileptonType dt);                      // Dilepton Type
    

  protected:
    
    // Set the average rates
    void loadRates();

    // Calculate the avg rate from rate histogram
    float getAvgRate(TH1* hist);

    // Get the rate given a signal region
    void getRealEff(DileptonType dt, susy::fake::Region fr, float &r1, float &r2);
    void getFakeRate(DileptonType dt, susy::fake::Region fr, float &f1, float &f2);

    // Get Number of Fakes
    float getWeight(float r1, float f1, float r2, float f2, float ntt, float ntl, float nlt, float nll){
      return getRF(r1,f1,r2,f2,ntt,ntl,nlt,nll) 
	   + getFR(r1,f1,r2,f2,ntt,ntl,nlt,nll) 
	   + getFF(r1,f1,r2,f2,ntt,ntl,nlt,nll);
    };
	
    // Number of Real-Real
    float getRR(float r1, float f1, float r2, float f2, float ntt, float ntl, float nlt, float nll){
      return r1*r2/(r1-f1)/(r2-f2)*( (1.-f1)*(1.-f2)*ntt
				     + (f1-1.)*f2*ntl
				     + f1*(f2-1.)*nlt
				     + f1*f2*nll);
    };
    // Number of Real-Fake
    float getRF(float r1, float f1, float r2, float f2, float ntt, float ntl, float nlt, float nll){
      return r1*f2/(r1-f1)/(r2-f2)*( (f1-1.)*(1.-r2)*ntt
				     + (1.-f1)*r2*ntl
				     + f1*(1.-r2)*nlt
				     - f1*r2*nll
				     );
    };
    // Number of Fake-Real
    float getFR(float r1, float f1, float r2, float f2, float ntt, float ntl, float nlt, float nll){
      return f1*r2/(r1-f1)/(r2-f2)*( (r1-1.)*(1.-f2)*ntt
				     + (1.-r1)*f2*ntl
				     + r1*(1.-f2)*nlt
				     - r1*f2*nll
				     );
    };
    // Number of Fake-Fake
    float getFF(float r1, float f1, float r2, float f2, float ntt, float ntl, float nlt, float nll){
      return f1*f2/(r1-f1)/(r2-f2)*( (1.-r1)*(1.-r2)*ntt
				     + (r1-1.)*r2*ntl
				     + r1*(r2-1.)*nlt
				     + r1*r2*nll
				     );
    };

    // Average rates
    float m_avgElReal[susy::fake::NumberOfSignalRegions];
    float m_avgMuReal[susy::fake::NumberOfSignalRegions];
    float m_avgElFake[susy::fake::NumberOfSignalRegions];
    float m_avgMuFake[susy::fake::NumberOfSignalRegions];

    TRandom3* rand;
    float m_upper;

  };    
} // fake
} // susy

#endif
