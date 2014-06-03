
#include "SameSignMatrixMethod/StatErrorTool.h"

#include "TH1.h"

using namespace FakeStatTool;

//---------------------------------------------------------------//
// Constructor
//---------------------------------------------------------------//
StatErrorTool::StatErrorTool(std::string fakeFile, 
			     SameSignMatrixMethod::RATE_PARAM rate_param_real_el,
			     SameSignMatrixMethod::RATE_PARAM rate_param_fake_el,
			     SameSignMatrixMethod::RATE_PARAM rate_param_real_mu,
			     SameSignMatrixMethod::RATE_PARAM rate_param_fake_mu
			     ) : 
  DiLeptonMatrixMethod(), // Could decouple these, but this is easiest for now
  rand(NULL),
  m_upper(1.29)
{

  // Range for error on prediction of 0 for 68%
  // 0-1.29 and this will be the range we sample 
  // from.  NOTE: I am doing this random.. should
  // it be weighted?
  rand = new TRandom3();

  // Load the rates
  configure(fakeFile, 
	    rate_param_real_el,
	    rate_param_fake_el,
	    rate_param_real_mu,
	    rate_param_fake_mu
	    ); 
  loadRates();
  
}

//---------------------------------------------------------------//  
// Destructor
//---------------------------------------------------------------//
StatErrorTool::~StatErrorTool() 
{

  if(rand) delete rand;

}

//---------------------------------------------------------------//
// Get the Updated Error
//---------------------------------------------------------------//
float StatErrorTool::getUpdatedError(float currentError, 
				    int n_tt,
				    int n_tl,
				    int n_lt,
				    int n_ll,
				    susy::fake::Region fr,
				    DileptonType dt)
{

  // Set the Seed
  int seed = (int) 1000* currentError * (3 + n_tt + n_tl + n_lt + n_ll);
  rand->SetSeed( seed );

  // Decide which combos are missing
  float ntt = n_tt == 0 ? rand->Rndm() * m_upper : 0;
  float ntl = n_tl == 0 ? rand->Rndm() * m_upper : 0;
  float nlt = n_lt == 0 ? rand->Rndm() * m_upper : 0;
  float nll = n_ll == 0 ? rand->Rndm() * m_upper : 0;

  // If we have all combos, just return the previous error
  if( !(ntt || ntl || nlt || nll) ) return currentError;
  
  // Get the avg rates for this event
  float r1 = 0;
  float r2 = 0;
  float f1 = 0;
  float f2 = 0;
  
  getRealEff(dt, fr, r1, r2);
  getFakeRate(dt, fr, f1, f2);

  // Now load the weights for the 'missing' events
  float wtt = getWeight(r1, f1, r2, f2, ntt, 0,   0,   0);
  float wtl = getWeight(r1, f1, r2, f2, 0,   ntl, 0,   0);
  float wlt = getWeight(r1, f1, r2, f2, 0,   0,   nlt, 0);
  float wll = getWeight(r1, f1, r2, f2, 0,   0,   0,   nll);

  // We treat the stat error the same as we do for mc weights
  // which is the Sqrt( Sum( w^2 ) )
  return sqrt( currentError*currentError 
	       + wtt*wtt 
	       + wtl*wtl 
	       + wlt*wlt 
	       + wll*wll );
  
}

//---------------------------------------------------------------//
// Load the average rates
//---------------------------------------------------------------//
void StatErrorTool::loadRates()
{

  // Load the rates
    for(int fr=0; fr<susy::fake::NumberOfSignalRegions; ++fr){
    m_avgElReal[fr] = getAvgRate(m_el_real_eff[fr]);
    m_avgMuReal[fr] = getAvgRate(m_mu_real_eff[fr]);
    m_avgElFake[fr] = getAvgRate(m_el_fake_rate[fr]);
    m_avgMuFake[fr] = getAvgRate(m_mu_fake_rate[fr]);
  }

}

//---------------------------------------------------------------//
// Get the average rate from histo
//---------------------------------------------------------------//
float StatErrorTool::getAvgRate(TH1* hist)
{

  // This will return the average fake rate. Not doing
  // anything fancy like taking weighted rates.

  float total = 0;
  int nBins   = hist->GetNbinsX();

  for(int bin = 1; bin<=nBins; ++bin) 
    total += hist->GetBinContent(bin);

  return total/nBins;

}

//---------------------------------------------------------------//
// Get the real effeciency
//---------------------------------------------------------------//
void StatErrorTool::getRealEff(DileptonType dt, susy::fake::Region fr, 
			       float &r1, float &r2)
{

  if(dt == DT_ee){
    r1 = m_avgElReal[fr];
    r2 = m_avgElReal[fr];
  }
  else if(dt == DT_mm){
    r1 = m_avgMuReal[fr];
    r2 = m_avgMuReal[fr];
  }
  else if(dt == DT_em){
    r1 = m_avgElReal[fr];  // Use electron as leading which
    r2 = m_avgMuReal[fr];  // and muon as subleading.. Better way??
  }
  else{
    r1 = 0; 
    r2 = 0;
    std::cout<<"StatErrorTool::getRealEff -- Unknown DileptonType!!\n\tr1=0\n\tr2=0"<<std::endl;
  }

}
//---------------------------------------------------------------//
// Get the fake rate
//---------------------------------------------------------------//
void StatErrorTool::getFakeRate(DileptonType dt, susy::fake::Region fr, 
				float &f1, float &f2)
{

  if(dt == DT_ee){
    f1 = m_avgElFake[fr];
    f2 = m_avgElFake[fr];
  }
  else if(dt == DT_mm){
    f1 = m_avgMuFake[fr];
    f2 = m_avgMuFake[fr];
  }
  else if(dt == DT_em){
    f1 = m_avgElFake[fr];  // Use electron as leading which
    f2 = m_avgMuFake[fr];  // and muon as subleading.. Better way??
  }
  else{
    f1 = 0; 
    f2 = 0;
    std::cout<<"StatErrorTool::getFakeEff -- Unknown DileptonType!!\n\tf1=0\n\tf2=0"<<std::endl;
  }

}
