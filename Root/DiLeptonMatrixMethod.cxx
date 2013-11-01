#include "SusyMatrixMethod/DiLeptonMatrixMethod.h"

#include <algorithm> // copy
#include <cassert>   // assert
#include <iterator>  // ostream_iterator

// -----------------------------------------------------------------------------
SusyMatrixMethod::DiLeptonMatrixMethod::DiLeptonMatrixMethod():
    m_hist_file(NULL),
    m_el_real_up(0.), m_el_real_down(0.),
    m_mu_real_up(0.), m_mu_real_down(0.),
    m_el_datamc(0.), m_mu_datamc(0.),
    m_el_region(0.), m_mu_region(0.),
    m_el_metrel(NULL), m_mu_metrel(NULL),
    m_el_eta(NULL), m_mu_eta(NULL),
    m_el_HFLFerr(0.)

{
  // do nothing
}

// -----------------------------------------------------------------------------
SusyMatrixMethod::DiLeptonMatrixMethod::~DiLeptonMatrixMethod()
{
  // Clean up pointers to member histograms
  std::cout << "checking m_hist_file\n";
  if (m_hist_file != NULL) {
    std::cout << "closing m_hist_file\n";
    m_hist_file->Close();
    std::cout << "deleting m_hist_file\n";
    // delete m_hist_file;
  }
  std::cout << "done with ~DiLeptonMatrixMethod\n";
}

// -----------------------------------------------------------------------------
bool SusyMatrixMethod::DiLeptonMatrixMethod::configure( std::string file_name
                                                      , RATE_PARAM rate_param_real_el
                                                      , RATE_PARAM rate_param_fake_el
                                                      , RATE_PARAM rate_param_real_mu
                                                      , RATE_PARAM rate_param_fake_mu
                                                      )
{
  // store the rate parameterization and systematic
  m_rate_param_real_el = rate_param_real_el;
  m_rate_param_fake_el = rate_param_fake_el;
  m_rate_param_real_mu = rate_param_real_mu;
  m_rate_param_fake_mu = rate_param_fake_mu;
  m_hist_file = TFile::Open(file_name.c_str(), "OPEN");
  if (!m_hist_file || !m_hist_file->IsOpen()) {
    std::cout << "ERROR: Failed to open fake rate file: "<< file_name << "\n";
    return false;
  }
  struct RetrieveOrNull {
    TFile *f_;  // TFile::Get breaks const-ness
    vector<string> miss_;
    RetrieveOrNull(TFile *f) : f_(f) {}
    bool anythingMissing() const { return miss_.size()>0; }
    TH1* operator() (const string &histoname) {
      cout<<"Retrieving "<<histoname<<endl;
      if(!f_) miss_.push_back(histoname);
      TH1* h = static_cast<TH1*>(f_->Get(histoname.c_str()));
      if(!h) miss_.push_back(histoname);
      return h;
    }
  } ron(m_hist_file);
  std::string s_el_real("el_real_eff_"), s_el_fake("el_fake_rate_");
  std::string s_mu_real("mu_real_eff_"), s_mu_fake("mu_fake_rate_");
  for(int r=0; r<FR_N; ++r) { // Get the Real eff and Fake rate for electrons and muons
    m_el_real_eff [r]  = ron(s_el_real + FRNames[r]);
    m_el_fake_rate[r]  = ron(s_el_fake + FRNames[r]);
    m_mu_real_eff [r]  = ron(s_mu_real + FRNames[r]);
    m_mu_fake_rate[r]  = ron(s_mu_fake + FRNames[r]);
  }
  if(ron.anythingMissing()) {
    cout<<"DiLeptonMatrixMethod::configure(): missing the following systematics:"<<endl;
    copy(ron.miss_.begin(), ron.miss_.end(), ostream_iterator<string>(cout,"\n\t"));
    cout<<"Exiting"<<endl;
    return false;
  }
  loadSysFromFile();
  return true;
}
// -----------------------------------------------------------------------------
float SusyMatrixMethod::DiLeptonMatrixMethod::getTotalFake(
    bool isTight1, bool isElectron1, float pt1, float eta1,
    bool isTight2, bool isElectron2, float pt2, float eta2,
    FAKE_REGION region,
    float MetRel,
    SYSTEMATIC syst) const
{
  SusyMatrixMethod::MatrixLepton lep1(isTight1, isElectron1, pt1, eta1);
  SusyMatrixMethod::MatrixLepton lep2(isTight2, isElectron2, pt2, eta2);

  return getTotalFake(lep1, lep2, region, MetRel, syst);
}

// -----------------------------------------------------------------------------
float SusyMatrixMethod::DiLeptonMatrixMethod::getRR(
    bool isTight1, bool isElectron1, float pt1, float eta1,
    bool isTight2, bool isElectron2, float pt2, float eta2,
    FAKE_REGION region,
    float MetRel,
    SYSTEMATIC syst) const
{
  SusyMatrixMethod::MatrixLepton lep1(isTight1, isElectron1, pt1, eta1);
  SusyMatrixMethod::MatrixLepton lep2(isTight2, isElectron2, pt2, eta2);

  return getRR(lep1, lep2, region, MetRel, syst);
}

// -----------------------------------------------------------------------------
float SusyMatrixMethod::DiLeptonMatrixMethod::getRF(
    bool isTight1, bool isElectron1, float pt1, float eta1,
    bool isTight2, bool isElectron2, float pt2, float eta2,
    FAKE_REGION region,
    float MetRel,
    SYSTEMATIC syst) const
{
  SusyMatrixMethod::MatrixLepton lep1(isTight1, isElectron1, pt1, eta1);
  SusyMatrixMethod::MatrixLepton lep2(isTight2, isElectron2, pt2, eta2);

  return getRF(lep1, lep2, region, MetRel, syst);
}

// -----------------------------------------------------------------------------
float SusyMatrixMethod::DiLeptonMatrixMethod::getFR(
    bool isTight1, bool isElectron1, float pt1, float eta1,
    bool isTight2, bool isElectron2, float pt2, float eta2,
    FAKE_REGION region,
    float MetRel,
    SYSTEMATIC syst) const
{
  SusyMatrixMethod::MatrixLepton lep1(isTight1, isElectron1, pt1, eta1);
  SusyMatrixMethod::MatrixLepton lep2(isTight2, isElectron2, pt2, eta2);

  return getFR(lep1, lep2, region, MetRel, syst);
}

// -----------------------------------------------------------------------------
float SusyMatrixMethod::DiLeptonMatrixMethod::getFF(
    bool isTight1, bool isElectron1, float pt1, float eta1,
    bool isTight2, bool isElectron2, float pt2, float eta2,
    FAKE_REGION region,
    float MetRel,
    SYSTEMATIC syst) const
{
  SusyMatrixMethod::MatrixLepton lep1(isTight1, isElectron1, pt1, eta1);
  SusyMatrixMethod::MatrixLepton lep2(isTight2, isElectron2, pt2, eta2);

  return getFF(lep1, lep2, region, MetRel, syst);
}

// -----------------------------------------------------------------------------
float SusyMatrixMethod::DiLeptonMatrixMethod::getTotalFake(
    const MatrixLepton& lep1,
    const MatrixLepton& lep2,
    FAKE_REGION region,
    float MetRel,
    SYSTEMATIC syst) const

{
  float fake_contribution = getRF(lep1, lep2, region, MetRel, syst);
  fake_contribution += getFR(lep1, lep2, region, MetRel, syst);
  fake_contribution += getFF(lep1, lep2, region, MetRel, syst);
  return fake_contribution;
}

// -----------------------------------------------------------------------------
float SusyMatrixMethod::DiLeptonMatrixMethod::getRR(
    const MatrixLepton& lep1,
    const MatrixLepton& lep2,
    FAKE_REGION region,
    float MetRel,
    SYSTEMATIC syst) const
{
  int n_tt = getTT(lep1, lep2);
  int n_tl = getTL(lep1, lep2);
  int n_lt = getLT(lep1, lep2);
  int n_ll = getLL(lep1, lep2);

  float r1 = getRate(lep1, REAL, region, MetRel, syst);
  float r2 = getRate(lep2, REAL, region, MetRel, syst);

  float f1 = getRate(lep1, FAKE, region, MetRel, syst);
  float f2 = getRate(lep2, FAKE, region, MetRel, syst);

  float rr = r1*r2/(r1-f1)/(r2-f2)*( (1.-f1)*(1.-f2)*n_tt
                                   + (f1-1.)*f2*n_tl
                                   + f1*(f2-1.)*n_lt
                                   + f1*f2*n_ll
                                   );

  int num_electrons = 0;
  if (lep1.isElectron()) ++num_electrons;
  if (lep2.isElectron()) ++num_electrons;
  std::string event_type = "";
  if (num_electrons == 2) event_type = "ee";
  else if (num_electrons == 1) event_type = "emu";
  else if (num_electrons == 0) event_type = "mumu";

  return rr;
}

// -----------------------------------------------------------------------------
float SusyMatrixMethod::DiLeptonMatrixMethod::getRF(
    const MatrixLepton& lep1,
    const MatrixLepton& lep2,
    FAKE_REGION region,
    float MetRel,
    SYSTEMATIC syst) const
{
  int n_tt = getTT(lep1, lep2);
  int n_tl = getTL(lep1, lep2);
  int n_lt = getLT(lep1, lep2);
  int n_ll = getLL(lep1, lep2);

  float r1 = getRate(lep1, REAL, region, MetRel, syst);
  float r2 = getRate(lep2, REAL, region, MetRel, syst);

  float f1 = getRate(lep1, FAKE, region, MetRel, syst);
  float f2 = getRate(lep2, FAKE, region, MetRel, syst);

  float rf = r1*f2/(r1-f1)/(r2-f2)*( (f1-1.)*(1.-r2)*n_tt
                                   + (1.-f1)*r2*n_tl
                                   + f1*(1.-r2)*n_lt
                                   - f1*r2*n_ll
                                   );
  return rf;
}

// -----------------------------------------------------------------------------
float SusyMatrixMethod::DiLeptonMatrixMethod::getFR(
    const MatrixLepton& lep1,
    const MatrixLepton& lep2,
    FAKE_REGION region,
    float MetRel,
    SYSTEMATIC syst) const
{
  int n_tt = getTT(lep1, lep2);
  int n_tl = getTL(lep1, lep2);
  int n_lt = getLT(lep1, lep2);
  int n_ll = getLL(lep1, lep2);

  float r1 = getRate(lep1, REAL, region, MetRel, syst);
  float r2 = getRate(lep2, REAL, region, MetRel, syst);

  float f1 = getRate(lep1, FAKE, region, MetRel, syst);
  float f2 = getRate(lep2, FAKE, region, MetRel, syst);

  float fr = f1*r2/(r1-f1)/(r2-f2)*( (r1-1.)*(1.-f2)*n_tt
                                   + (1.-r1)*f2*n_tl
                                   + r1*(1.-f2)*n_lt
                                   - r1*f2*n_ll
                                   );
  return fr;
}

// -----------------------------------------------------------------------------
float SusyMatrixMethod::DiLeptonMatrixMethod::getFF(
    const MatrixLepton& lep1,
    const MatrixLepton& lep2,
    FAKE_REGION region,
    float MetRel,
    SYSTEMATIC syst) const
{
  int n_tt = getTT(lep1, lep2);
  int n_tl = getTL(lep1, lep2);
  int n_lt = getLT(lep1, lep2);
  int n_ll = getLL(lep1, lep2);

  float r1 = getRate(lep1, REAL, region, MetRel, syst);
  float r2 = getRate(lep2, REAL, region, MetRel, syst);

  float f1 = getRate(lep1, FAKE, region, MetRel, syst);
  float f2 = getRate(lep2, FAKE, region, MetRel, syst);

  float ff = f1*f2/(r1-f1)/(r2-f2)*( (1.-r1)*(1.-r2)*n_tt
                                   + (r1-1.)*r2*n_tl
                                   + r1*(r2-1.)*n_lt
                                   + r1*r2*n_ll
                                   );

  return ff;
}

// -----------------------------------------------------------------------------
float SusyMatrixMethod::DiLeptonMatrixMethod::getRate(
    bool isTight, bool isElectron, float pt, float eta,
    RATE_TYPE rate_type, FAKE_REGION region, float MetRel,
    SYSTEMATIC syst) const
{
  SusyMatrixMethod::MatrixLepton lep(isTight, isElectron, pt, eta);
  return getRate(lep, rate_type, region, MetRel, syst);
}


// -----------------------------------------------------------------------------
float SusyMatrixMethod::DiLeptonMatrixMethod::getRate(
    const MatrixLepton& lep,
    RATE_TYPE rate_type,
    FAKE_REGION region,
    float MetRel,
    SYSTEMATIC syst) const
{
  // get histogram from member objects
  TH1* h_rate = NULL;
  RATE_PARAM rate_param = SusyMatrixMethod::PT;

  if (rate_type == REAL) {
    if (lep.isElectron()) {
      h_rate = m_el_real_eff[region];
      rate_param = m_rate_param_real_el;
    }
    else {
      h_rate = m_mu_real_eff[region];
      rate_param = m_rate_param_real_mu;
    }
  }
  else if (rate_type == FAKE) {
    if (lep.isElectron()) {
      h_rate = m_el_fake_rate[region];
      rate_param = m_rate_param_fake_el;
    }
    else {
      h_rate = m_mu_fake_rate[region];
      rate_param = m_rate_param_fake_mu;
    }
  }
  else {
    return 0;
  }

  // Get the rate bin
  int rate_bin = getRateBin(lep, h_rate, rate_param);

  /*
  if (m_rate_param == PT_ETA) {
    int maxbin = h_rate->GetYaxis()->GetNbins();
    float max  = h_rate->GetYaxis()->GetBinCenter(maxbin) +
                 h_rate->GetYaxis()->GetBinWidth(maxbin)/2.;
    // float pt   = lep.pt() > max ? max - 1e-4 : lep.pt();
    float pt   = lep.pt()/1000.;
    // TODO Clean up conversion
    pt = pt > max ? max - 1e-4 : pt;
    rate_bin   = h_rate->FindBin(pt, fabs(lep.eta()));
  }
  else if (m_rate_param == PT) {
    int maxbin = h_rate->GetXaxis()->GetNbins();
    float max  = h_rate->GetXaxis()->GetBinCenter(maxbin) +
                 h_rate->GetXaxis()->GetBinWidth(maxbin)/2.;
    // float pt   = lep.pt() > max ? max - 1e-4 : lep.pt();
    float pt   = lep.pt()/1000.;
    // TODO Clean up conversion
    pt = pt > max ? max - 1e-4 : pt;
    rate_bin   = h_rate->FindBin(pt);
  }
  */

  float rate = h_rate->GetBinContent(rate_bin);
  rate *= (1+getRateSyst(lep, rate_type, region, MetRel, syst));

  // Don't let rate go above 1
  if( rate > 1 ) rate = 1.;
  // Don't let the rate go below 0
  if( rate < 0 ) rate = 1e-5;

  return rate;
}

// -----------------------------------------------------------------------------
int SusyMatrixMethod::DiLeptonMatrixMethod::getRateBin( const MatrixLepton& lep,
							TH1* h_rate,
							RATE_PARAM rate_param) const
{

  // Handle 2-D param
  if( rate_param == PT_ETA ){
    int max_bin_x = h_rate->GetXaxis()->GetNbins();
    float max_x   = h_rate->GetXaxis()->GetBinCenter(max_bin_x);
    int max_bin_y = h_rate->GetYaxis()->GetNbins();
    float max_y   = h_rate->GetYaxis()->GetBinCenter(max_bin_y);

    // Reset x-var
    float var_x = lep.pt() / 1000.;
    var_x = var_x > max_x ? max_x : var_x;
    
    // Reset y-var
    float var_y = fabs(lep.eta());
    var_y = var_y > max_y ? max_y : var_y;

    return h_rate->FindBin(var_x,var_y);

  }
  
  // Handle 1-D param
  if( rate_param == PT ){
    int max_bin_x = h_rate->GetXaxis()->GetNbins();
    float max_x   = h_rate->GetXaxis()->GetBinCenter(max_bin_x);

    // Reset x-var
    float var_x = lep.pt() / 1000.;
    var_x = var_x > max_x ? max_x : var_x;

    return h_rate->FindBin(var_x);

  }

  // If we are here, then we don't have a supported
  // parameterization set. Print error message
  std::cout<<"*** Error: Did not set a supported rate-parameterization"<<std::endl;
  return 0;
  
}

// -----------------------------------------------------------------------------
float SusyMatrixMethod::DiLeptonMatrixMethod::getRateSyst(
    const MatrixLepton& lep,
    RATE_TYPE rate_type,
    FAKE_REGION region,
    float MetRel,
    SYSTEMATIC syst) const
{
  // No Systematic
  if( syst == SYS_NONE )
    return 0.;

  if ( syst == SYS_N || syst == SYS_N_USER) {
    std::cout << "WARNING: invalid SYSTEMATIC type\n";
  }

  // Real Eff Sys
  if( rate_type == REAL ){

    float statSys = getStatError(lep, rate_type, region);

    if( lep.isElectron() ){
      if( syst == SYS_EL_RE_UP )   return sqrt( m_el_real_up*m_el_real_up + statSys*statSys);
      if( syst == SYS_EL_RE_DOWN ) return -sqrt( m_el_real_down*m_el_real_down + statSys*statSys);
    }
    else{
      if( syst == SYS_MU_RE_UP )   return sqrt( m_mu_real_up*m_mu_real_up + statSys*statSys);
      if( syst == SYS_MU_RE_DOWN ) return -sqrt( m_mu_real_down*m_mu_real_down + statSys*statSys);
    }
    return 0.;
  }

  // Fake Rate Sys
  if( rate_type == FAKE ){
    if( lep.isElectron() ){
      // Fake Rate Up
      if( syst == SYS_EL_FR_UP ){
        //float metSys   = m_el_metrel->GetBinContent( m_el_metrel->FindBin(MetRel/1000.) );
        //metSys = metSys > 0 ? metSys : 0.;
	float etaSys   = m_el_eta->GetBinContent( m_el_eta->FindBin(fabs(lep.eta())) );
	etaSys = etaSys > 0 ? etaSys : 0.;
        float statSys  = getStatError(lep, rate_type, region);

        // treat systematics as independent, and add in quadrature
        return sqrt( m_el_HFLFerr*m_el_HFLFerr
		     + statSys*statSys
		     + etaSys*etaSys
		     + m_el_datamc*m_el_datamc
		     + m_el_region*m_el_region
		     );
      }

      // Fake Rate Down
      if( syst == SYS_EL_FR_DOWN){
        //float metSys   = m_el_metrel->GetBinContent( m_el_metrel->FindBin(MetRel/1000.) );
        //metSys = metSys < 0 ? metSys : 0.;
	float etaSys   = m_el_eta->GetBinContent( m_el_eta->FindBin(fabs(lep.eta())) );
	etaSys = etaSys < 0 ? etaSys : 0.;
        float statSys  = getStatError(lep, rate_type, region);

        // treat systematics as independent, and add in quadrature
        return (-1) * sqrt(m_el_HFLFerr*m_el_HFLFerr 
			   + statSys*statSys
			   + etaSys*etaSys
			   + m_el_datamc*m_el_datamc
			   + m_el_region*m_el_region
                          );
      }

      // HF/LF error
      if( syst == SYS_EL_HFLF_UP )   return m_el_HFLFerr;
      if( syst == SYS_EL_HFLF_DOWN ) return -m_el_HFLFerr;

      // Met Rel
      //if( syst == SYS_EL_METREL ) 
      //return m_el_metrel->GetBinContent( m_el_metrel->FindBin(MetRel/1000.) );	

      // Eta
      if( syst == SYS_EL_ETA )
	return m_el_eta->GetBinContent( m_el_eta->FindBin(fabs(lep.eta())) );

      // Statistical error from maps
      if( syst == SYS_EL_FR_STAT_UP )   return getStatError(lep, rate_type, region);
      if( syst == SYS_EL_FR_STAT_DOWN ) return (-1) * getStatError(lep, rate_type, region);

      // Data/MC error
      if( syst == SYS_EL_DATAMC_UP )   return m_el_datamc;
      if( syst == SYS_EL_DATAMC_DOWN ) return -m_el_datamc;

      // Error from combination
      if( syst == SYS_EL_REG_UP )   return m_el_region;
      if( syst == SYS_EL_REG_DOWN ) return -m_el_region;
      
      // Neither shift
      return 0.;

    }// end if Electron
    else{
      // Fake Rate Up
      if( syst == SYS_MU_FR_UP ){
        //float metSys   = m_mu_metrel->GetBinContent( m_mu_metrel->FindBin(MetRel/1000.) );
        //metSys = metSys > 0 ? metSys : 0.;
	float etaSys   = m_mu_eta->GetBinContent( m_mu_eta->FindBin(fabs(lep.eta())) );
	etaSys = etaSys > 0 ? etaSys : 0.;
        float statSys  = getStatError(lep, rate_type, region);

        // treat systematics as independent, and add in quadrature
        return sqrt( statSys*statSys
		     + etaSys*etaSys
		     + m_mu_datamc*m_mu_datamc
		     + m_mu_region*m_mu_region
		     );
      }

      // Fake Rate Down
      if( syst == SYS_MU_FR_DOWN){
        //float metSys   = m_mu_metrel->GetBinContent( m_mu_metrel->FindBin(MetRel/1000.) );
        //metSys = metSys < 0 ? metSys : 0.;
	float etaSys   = m_mu_eta->GetBinContent( m_mu_eta->FindBin(fabs(lep.eta())) );
	etaSys = etaSys < 0 ? etaSys : 0.;
        float statSys  = getStatError(lep, rate_type, region);

        // treat systematics as independent, and add in quadrature
        return (-1) * sqrt( statSys*statSys
			    + etaSys*etaSys
			    + m_mu_datamc*m_mu_datamc
			    + m_mu_region*m_mu_region
			    );
      }

      // Met Rel
      //if( syst == SYS_MU_METREL )
      //return m_mu_metrel->GetBinContent( m_mu_metrel->FindBin(MetRel/1000.) );

      // Eta
      if( syst == SYS_MU_ETA )
	return m_mu_eta->GetBinContent( m_mu_eta->FindBin(fabs(lep.eta())) );

      // Statistical error from maps
      if( syst == SYS_MU_FR_STAT_UP )   return getStatError(lep, rate_type, region);
      if( syst == SYS_MU_FR_STAT_DOWN ) return (-1) * getStatError(lep, rate_type, region);

      // Data/MC error
      if( syst == SYS_MU_DATAMC_UP )   return m_mu_datamc;
      if( syst == SYS_MU_DATAMC_DOWN ) return -m_mu_datamc;
      
      // Error from combination
      if( syst == SYS_MU_REG_UP )   return m_mu_region;
      if( syst == SYS_MU_REG_DOWN ) return -m_mu_region;

      // Neither shift
      return 0.;

    }// end if Muon
  }// end if FAKE

  return 0.;
}

// -----------------------------------------------------------------------------
int SusyMatrixMethod::DiLeptonMatrixMethod::getTT( const MatrixLepton& lep1
                                                 , const MatrixLepton& lep2
                                                 ) const
{
  if (lep1.isTight() && lep2.isTight()) return 1;
  return 0;
}

// -----------------------------------------------------------------------------
int SusyMatrixMethod::DiLeptonMatrixMethod::getTL( const MatrixLepton& lep1
                                                 , const MatrixLepton& lep2
                                                 ) const
{
  if (lep1.isTight() && !lep2.isTight()) return 1;
  return 0;
}

// -----------------------------------------------------------------------------
int SusyMatrixMethod::DiLeptonMatrixMethod::getLT( const MatrixLepton& lep1
                                                 , const MatrixLepton& lep2
                                                 ) const
{
  if (!lep1.isTight() && lep2.isTight()) return 1;
  return 0;
}

// -----------------------------------------------------------------------------
int SusyMatrixMethod::DiLeptonMatrixMethod::getLL( const MatrixLepton& lep1
                                                 , const MatrixLepton& lep2
                                                 ) const
{
  if (!lep1.isTight() && !lep2.isTight()) return 1;
  return 0;
}

// -----------------------------------------------------------------------------
void SusyMatrixMethod::DiLeptonMatrixMethod::loadSysFromFile()
{
  // I wanted to keep the systematics available as being separate in case
  // we want to investigate the impact of a specific systematic.  Also
  // I want to parameterize the metrel variation, as it kind of bounces around
  // and it isn't fair to assign +/- the max and min bounce.

  struct RetrieveOrZero {
    typedef TParameter<double>* tdp_t;
    TFile *f_;  // TFile::Get breaks const-ness
    vector<string> miss_;
    RetrieveOrZero(TFile *f) : f_(f) {}
    bool anythingMissing() const { return miss_.size()>0; }
    double operator() (const char* varname) {
      double val(0.0);
      if(!f_) miss_.push_back(varname);
      else if(tdp_t tdp = static_cast<tdp_t>(f_->Get(varname))) val = tdp->GetVal();
      else miss_.push_back(varname);
      return val;
    }
  } roz(m_hist_file);
  m_el_real_up   = roz("el_real_up");
  m_el_real_down = roz("el_real_down");
  m_mu_real_up   = roz("mu_real_up");
  m_mu_real_down = roz("mu_real_down");
  m_el_HFLFerr   = roz("el_HFLFerr");
  m_el_datamc    = roz("el_datamc");
  m_mu_datamc    = roz("mu_datamc");
  m_el_region    = roz("el_region");
  m_mu_region    = roz("mu_region");
  if(roz.anythingMissing()) {
    cout<<"DiLeptonMatrixMethod::loadSysFromFile(): missing the following systematics:"<<endl;
    copy(roz.miss_.begin(), roz.miss_.end(), ostream_iterator<string>(cout,"\n\t"));
    cout<<"Will use 0.0 for these"<<endl;
    assert(false); // should provide a bool return value and exit gracefully
  }
  // Met Rel -- Parameterized vs Met Rel
  m_el_metrel = static_cast<TH1F*>(m_hist_file->Get("el_metrel_sys"));
  m_mu_metrel = static_cast<TH1F*>(m_hist_file->Get("mu_metrel_sys"));
  // Eta -- Parameterized vs Eta
  m_el_eta = static_cast<TH1F*>(m_hist_file->Get("el_eta_sys"));
  m_mu_eta = static_cast<TH1F*>(m_hist_file->Get("mu_eta_sys"));

}
// -----------------------------------------------------------------------------
float SusyMatrixMethod::DiLeptonMatrixMethod::getStatError(const MatrixLepton& lep
                 , RATE_TYPE rate_type
                 , FAKE_REGION region) const

{

  // get histogram from member objects
  TH1* h_rate = NULL;
  RATE_PARAM rate_param = PT;

  if (rate_type == REAL) {
    if (lep.isElectron()){
      h_rate = m_el_real_eff[region];
      rate_param = m_rate_param_real_el;
    }
    else{
      h_rate = m_mu_real_eff[region];
      rate_param = m_rate_param_real_mu;
    }
  }
  else if (rate_type == FAKE) {
    if (lep.isElectron()){
      h_rate = m_el_fake_rate[region];
      rate_param = m_rate_param_fake_el;
    }
    else{
      h_rate = m_mu_fake_rate[region];
      rate_param = m_rate_param_fake_mu;
    }
  }
  else {
    std::cout << "WARNING: invalid RATE_TYPE\n";
    return 0;
  }

  int rate_bin = getRateBin(lep, h_rate, rate_param);

  // Need to add in protection in the event that the
  // lepton Pt goes beyond the rate histograms boundaries
  /*
  if (m_rate_param == PT_ETA) {
    int maxbin = h_rate->GetYaxis()->GetNbins();
    float max  = h_rate->GetYaxis()->GetBinCenter(maxbin) +
                 h_rate->GetYaxis()->GetBinWidth(maxbin)/2.;
    // float pt   = lep.pt() > max ? max - 1e-4 : lep.pt();
    float pt   = lep.pt()/1000.;
    // TODO Clean up conversion
    pt = pt > max ? max - 1e-4 : pt;
    rate_bin   = h_rate->FindBin(lep.eta(), pt);
  }
  else if (m_rate_param == PT) {
    int maxbin = h_rate->GetXaxis()->GetNbins();
    float max  = h_rate->GetXaxis()->GetBinCenter(maxbin) +
                 h_rate->GetXaxis()->GetBinWidth(maxbin)/2.;
    // float pt   = lep.pt() > max ? max - 1e-4 : lep.pt();
    float pt   = lep.pt()/1000.;
    // TODO Clean up conversion
    pt = pt > max ? max - 1e-4 : pt;
    rate_bin   = h_rate->FindBin(pt);
  }
  */

  float error = h_rate->GetBinError(rate_bin);

  return error;

}
// -----------------------------------------------------------------------------
void SusyMatrixMethod::DiLeptonMatrixMethod::printInfo(
    const MatrixLepton& lep1,
    const MatrixLepton& lep2,
    FAKE_REGION region, float MetRel,
    SYSTEMATIC syst) const
{
  float r1 = getRate(lep1, REAL, region, MetRel, syst);
  float r2 = getRate(lep2, REAL, region, MetRel, syst);

  float f1 = getRate(lep1, FAKE, region, MetRel, syst);
  float f2 = getRate(lep2, FAKE, region, MetRel, syst);

  int num_electrons = 0;
  if (lep1.isElectron()) ++num_electrons;
  if (lep2.isElectron()) ++num_electrons;

  std::string event_type = "";
  if (num_electrons == 2) event_type = "ee";
  else if (num_electrons == 1) event_type = "emu";
  else if (num_electrons == 0) event_type = "mumu";

  std::cout << "-----------------------------------------------\n";
  std::cout << "region: " << region << "\t-\tsyst: " << syst << "\n";
  std::cout << "lep 1:"
            << "\telectron: " << lep1.isElectron()
            << "\ttight: "    << lep1.isTight()
            << "\tpt: "       << lep1.pt()
            << "\tr: "        << r1
            << "\tf: "        << f1
            << "\n";
  std::cout << "lep 2:"
            << "\telectron: " << lep2.isElectron()
            << "\ttight: "    << lep2.isTight()
            << "\tpt: "       << lep2.pt()
            << "\tr: "        << r2
            << "\tf: "        << f2
            << "\n";
  std::cout << "Met Rel: "    << MetRel << "\n";
  std::cout << "event type: " << event_type << "\n";
  std::cout << "rr: " << getRR(lep1,lep2,region,syst) << "\n";
  std::cout << "rf: " << getRF(lep1,lep2,region,syst) << "\n";
  std::cout << "fr: " << getFR(lep1,lep2,region,syst) << "\n";
  std::cout << "ff: " << getFF(lep1,lep2,region,syst) << "\n";
  std::cout << "total fake: " << getTotalFake(lep1,lep2,region,syst) << "\n";
}
