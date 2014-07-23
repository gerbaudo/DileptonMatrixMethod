#include "DileptonMatrixMethod/DiLeptonMatrixMethod.h"

#include "TParameter.h"
#include "TVectorD.h"
#include "TArrayD.h"
#include "TAxis.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2.h"


#include <algorithm> // copy
#include <cassert>   // assert
#include <iterator>  // ostream_iterator, distance
#include <sstream>

using DileptonMatrixMethod::DiLeptonMatrixMethod;

// -----------------------------------------------------------------------------
DileptonMatrixMethod::DiLeptonMatrixMethod::DiLeptonMatrixMethod():
    m_hist_file(NULL),
    m_el_frac_up(NULL),
    m_el_frac_do(NULL),
    m_mu_frac_up(NULL),
    m_mu_frac_do(NULL),
    m_el_real_up(0.0), m_el_real_down(0.0),
    m_mu_real_up(0.0), m_mu_real_down(0.0),
    m_el_datamc(0.0), m_mu_datamc(0.0),
    m_el_region(0.0), m_mu_region(0.0),
    m_el_metrel(NULL), m_mu_metrel(NULL),
    m_el_eta(NULL), m_mu_eta(NULL),
    m_el_HFLFerr(0.0)

{
  // do nothing
}

// -----------------------------------------------------------------------------
DileptonMatrixMethod::DiLeptonMatrixMethod::~DiLeptonMatrixMethod()
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
bool DileptonMatrixMethod::DiLeptonMatrixMethod::configure( std::string file_name
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
  if(rate_param_real_el != rate_param_fake_el ||
     rate_param_fake_el != rate_param_real_mu ||
     rate_param_real_mu != rate_param_fake_mu) {
      cout<<"DiLeptonMatrixMethod::configure: mixed parametrization not implemented"<<endl;
      return false;
  }
  const RATE_PARAM rp = rate_param_real_el;
  std::string s_el_real("el_real_eff_"), s_el_fake("el_fake_rate_");
  std::string s_mu_real("mu_real_eff_"), s_mu_fake("mu_fake_rate_");
  if(rp==PT_ETA) { // todo: use string replacement (eff->eff2d, rate->rate2d)
      s_el_real = "el_real_eff2d_"; s_el_fake = "el_fake_rate2d_";
      s_mu_real = "mu_real_eff2d_"; s_mu_fake = "mu_fake_rate2d_";
  }
  for(int r=0; r<susy::fake::NumberOfSignalRegions; ++r) { // Get the Real eff and Fake rate for electrons and muons
    std::string regionName(susy::fake::region2str(susy::fake::SignalRegions[r]));
    m_el_real_eff [r]  = ron(s_el_real + regionName);
    m_el_fake_rate[r]  = ron(s_el_fake + regionName);
    m_mu_real_eff [r]  = ron(s_mu_real + regionName);
    m_mu_fake_rate[r]  = ron(s_mu_fake + regionName);
  }
  m_el_frac_up = ron("el_fake_rate2d_CR_SSInc1j_frac_up");
  m_el_frac_do = ron("el_fake_rate2d_CR_SSInc1j_frac_do");
  m_mu_frac_up = ron("mu_fake_rate2d_CR_SSInc1j_frac_up");
  m_mu_frac_do = ron("mu_fake_rate2d_CR_SSInc1j_frac_do");
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
float DileptonMatrixMethod::DiLeptonMatrixMethod::getTotalFake(
    bool isTight1, bool isElectron1, float pt1, float eta1,
    bool isTight2, bool isElectron2, float pt2, float eta2,
    susy::fake::Region region,
    float MetRel,
    SYSTEMATIC syst) const
{
  DileptonMatrixMethod::MatrixLepton lep1(isTight1, isElectron1, pt1, eta1);
  DileptonMatrixMethod::MatrixLepton lep2(isTight2, isElectron2, pt2, eta2);

  return getTotalFake(lep1, lep2, region, MetRel, syst);
}

// -----------------------------------------------------------------------------
float DileptonMatrixMethod::DiLeptonMatrixMethod::getRR(
    bool isTight1, bool isElectron1, float pt1, float eta1,
    bool isTight2, bool isElectron2, float pt2, float eta2,
    susy::fake::Region region,
    float MetRel,
    SYSTEMATIC syst) const
{
  DileptonMatrixMethod::MatrixLepton lep1(isTight1, isElectron1, pt1, eta1);
  DileptonMatrixMethod::MatrixLepton lep2(isTight2, isElectron2, pt2, eta2);

  return getRR(lep1, lep2, region, MetRel, syst);
}

// -----------------------------------------------------------------------------
float DileptonMatrixMethod::DiLeptonMatrixMethod::getRF(
    bool isTight1, bool isElectron1, float pt1, float eta1,
    bool isTight2, bool isElectron2, float pt2, float eta2,
    susy::fake::Region region,
    float MetRel,
    SYSTEMATIC syst) const
{
  DileptonMatrixMethod::MatrixLepton lep1(isTight1, isElectron1, pt1, eta1);
  DileptonMatrixMethod::MatrixLepton lep2(isTight2, isElectron2, pt2, eta2);

  return getRF(lep1, lep2, region, MetRel, syst);
}

// -----------------------------------------------------------------------------
float DileptonMatrixMethod::DiLeptonMatrixMethod::getFR(
    bool isTight1, bool isElectron1, float pt1, float eta1,
    bool isTight2, bool isElectron2, float pt2, float eta2,
    susy::fake::Region region,
    float MetRel,
    SYSTEMATIC syst) const
{
  DileptonMatrixMethod::MatrixLepton lep1(isTight1, isElectron1, pt1, eta1);
  DileptonMatrixMethod::MatrixLepton lep2(isTight2, isElectron2, pt2, eta2);

  return getFR(lep1, lep2, region, MetRel, syst);
}

// -----------------------------------------------------------------------------
float DileptonMatrixMethod::DiLeptonMatrixMethod::getFF(
    bool isTight1, bool isElectron1, float pt1, float eta1,
    bool isTight2, bool isElectron2, float pt2, float eta2,
    susy::fake::Region region,
    float MetRel,
    SYSTEMATIC syst) const
{
  DileptonMatrixMethod::MatrixLepton lep1(isTight1, isElectron1, pt1, eta1);
  DileptonMatrixMethod::MatrixLepton lep2(isTight2, isElectron2, pt2, eta2);

  return getFF(lep1, lep2, region, MetRel, syst);
}

// -----------------------------------------------------------------------------
float DileptonMatrixMethod::DiLeptonMatrixMethod::getTotalFake(
    const MatrixLepton& lep1,
    const MatrixLepton& lep2,
    susy::fake::Region region,
    float MetRel,
    SYSTEMATIC syst) const
{
  float fake_contribution = getRF(lep1, lep2, region, MetRel, syst);
  fake_contribution += getFR(lep1, lep2, region, MetRel, syst);
  fake_contribution += getFF(lep1, lep2, region, MetRel, syst);
  return fake_contribution;
}

// -----------------------------------------------------------------------------
float DileptonMatrixMethod::DiLeptonMatrixMethod::getRR(
    const MatrixLepton& lep1,
    const MatrixLepton& lep2,
    susy::fake::Region region,
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
float DileptonMatrixMethod::DiLeptonMatrixMethod::getRF(
    const MatrixLepton& lep1,
    const MatrixLepton& lep2,
    susy::fake::Region region,
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
float DileptonMatrixMethod::DiLeptonMatrixMethod::getFR(
    const MatrixLepton& lep1,
    const MatrixLepton& lep2,
    susy::fake::Region region,
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
float DileptonMatrixMethod::DiLeptonMatrixMethod::getFF(
    const MatrixLepton& lep1,
    const MatrixLepton& lep2,
    susy::fake::Region region,
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
float DileptonMatrixMethod::DiLeptonMatrixMethod::getRate(
    bool isTight, bool isElectron, float pt, float eta,
    RATE_TYPE rate_type, susy::fake::Region region, float MetRel,
    SYSTEMATIC syst) const
{
  DileptonMatrixMethod::MatrixLepton lep(isTight, isElectron, pt, eta);
  return getRate(lep, rate_type, region, MetRel, syst);
}

// --------------------------------------------------------
bool DileptonMatrixMethod::DiLeptonMatrixMethod::getHistoAndParametrization(const MatrixLepton &lep,
                                                                        const susy::fake::Region reg,
                                                                        const RATE_TYPE &rt,
                                                                        TH1* &h, RATE_PARAM &rp) const
{
    bool found=false;
    int iRegion(getIndexRegion(reg));
    if(rt==REAL) {
        if (lep.isElectron()) {
            h  = m_el_real_eff[iRegion];
            rp = m_rate_param_real_el;
            found = true;
        } else {
            h  = m_mu_real_eff[iRegion];
            rp = m_rate_param_real_mu;
            found = true;
        }
    } else if(rt==FAKE) {
        if (lep.isElectron()) {
            h  = m_el_fake_rate[iRegion];
            rp = m_rate_param_fake_el;
            found = true;
        } else {
            h  = m_mu_fake_rate[iRegion];
            rp = m_rate_param_fake_mu;
            found = true;
        }
    }
    return found;
}
// -----------------------------------------------------------------------------
float DileptonMatrixMethod::DiLeptonMatrixMethod::getRate(
    const MatrixLepton& lep,
    RATE_TYPE rate_type,
    susy::fake::Region region,
    float MetRel,
    SYSTEMATIC syst) const
{
    float rate = 0.0;
    TH1* h_rate = NULL;
    RATE_PARAM rate_param = m_rate_param_real_el; // all 4 must be the same; checked at configuration (DG todo improve this)
    if(getHistoAndParametrization(lep, region, rate_type, h_rate, rate_param)) {
        int rate_bin = getRateBin(lep, h_rate, rate_param);
        rate = h_rate->GetBinContent(rate_bin);
        rate *= (1.0+getRateSyst(lep, rate_type, region, MetRel, syst));
        if( rate > 1.0 ) rate = 1.0;
        if( rate < 0.0 ) rate = 1e-5;
    } else {
        cout<<"DiLeptonMatrixMethod::getRate: cannot get histo"<<endl;
    }
  return rate;
}

// -----------------------------------------------------------------------------
int DileptonMatrixMethod::DiLeptonMatrixMethod::getRateBin( const MatrixLepton& lep,
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
    float var_x = lep.pt() / 1000.0;
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
    float var_x = lep.pt() / 1000.0;
    var_x = var_x > max_x ? max_x : var_x;

    return h_rate->FindBin(var_x);

  }

  // If we are here, then we don't have a supported
  // parameterization set. Print error message
  std::cout<<"*** Error: Did not set a supported rate-parameterization"<<std::endl;
  return 0;

}

// -----------------------------------------------------------------------------
float DileptonMatrixMethod::DiLeptonMatrixMethod::getRateSyst(
    const MatrixLepton& lep,
    RATE_TYPE rate_type,
    susy::fake::Region region,
    float MetRel,
    SYSTEMATIC syst) const
{
  if( syst == SYS_NOM ) return 0.0;
  if ( syst == SYS_N || syst == SYS_N_USER) { std::cout << "WARNING: invalid SYSTEMATIC type\n"; }
  if( rate_type == REAL ){ // Real Eff Sys
      float statSys = getRelStatError(lep, rate_type, region);
      if( lep.isElectron() ){
          if( syst == SYS_EL_RE_UP )   return sqrt( m_el_real_up*m_el_real_up + statSys*statSys);
          if( syst == SYS_EL_RE_DOWN ) return -sqrt( m_el_real_down*m_el_real_down + statSys*statSys);
      } else {
          if( syst == SYS_MU_RE_UP )   return sqrt( m_mu_real_up*m_mu_real_up + statSys*statSys);
          if( syst == SYS_MU_RE_DOWN ) return -sqrt( m_mu_real_down*m_mu_real_down + statSys*statSys);
    }
      return 0.;
  }
  if( rate_type == FAKE ){ // Fake Rate Sys
      if(syst==SYS_EL_FRAC_UP || syst==SYS_EL_FRAC_DO || syst==SYS_MU_FRAC_UP || syst==SYS_MU_FRAC_DO)
          return getFracRelativeError(lep, rate_type, region, syst);
      if( lep.isElectron() ){
          if( syst == SYS_EL_FR_UP ){ // Fake Rate Up
              float etaSys   = m_el_eta->GetBinContent( m_el_eta->FindBin(fabs(lep.eta())) );
              etaSys = etaSys > 0 ? etaSys : 0.;
              float statSys = getRelStatError(lep, rate_type, region);
              return +1.0*sqrt(statSys*statSys);
//               return sqrt( m_el_HFLFerr*m_el_HFLFerr //add in quadrature
//                            + statSys*statSys
//                            + etaSys*etaSys
//                            + m_el_datamc*m_el_datamc
//                            + m_el_region*m_el_region);
          }
          if( syst == SYS_EL_FR_DOWN){ // Fake Rate Down
              float etaSys   = m_el_eta->GetBinContent( m_el_eta->FindBin(fabs(lep.eta())) );
              etaSys = etaSys < 0.0 ? etaSys : 0.0;
              float statSys  = getRelStatError(lep, rate_type, region);
              return -1.0*sqrt(statSys*statSys);
//               return (-1.0) * sqrt(m_el_HFLFerr*m_el_HFLFerr // add in quadrature
//                                    + statSys*statSys
//                                    + etaSys*etaSys
//                                    + m_el_datamc*m_el_datamc
//                                    + m_el_region*m_el_region);
          }
          if( syst == SYS_EL_HFLF_UP )   return +m_el_HFLFerr; // HF/LF error
          if( syst == SYS_EL_HFLF_DOWN ) return -m_el_HFLFerr;
          if( syst == SYS_EL_ETA ) return m_el_eta->GetBinContent( m_el_eta->FindBin(fabs(lep.eta())) );
          if( syst == SYS_EL_FR_STAT_UP )   return +1.0*getRelStatError(lep, rate_type, region); // Statistical error from maps
          if( syst == SYS_EL_FR_STAT_DOWN ) return -1.0*getRelStatError(lep, rate_type, region);
          if( syst == SYS_EL_DATAMC_UP )   return +m_el_datamc; // Data/MC error
          if( syst == SYS_EL_DATAMC_DOWN ) return -m_el_datamc;
          if( syst == SYS_EL_REG_UP )   return +m_el_region; // Error from combinatino
          if( syst == SYS_EL_REG_DOWN ) return -m_el_region;
          return 0.;
      } else { // end electron
          if( syst == SYS_MU_FR_UP ){ // Fake Rate Up
              float etaSys   = m_mu_eta->GetBinContent( m_mu_eta->FindBin(fabs(lep.eta())) );
              etaSys = etaSys > 0.0 ? etaSys : 0.0;
              float statSys  = getRelStatError(lep, rate_type, region);
              return +1.0*sqrt(statSys*statSys);
//               return sqrt( statSys*statSys
//                            + etaSys*etaSys
//                            + m_mu_datamc*m_mu_datamc
//                            + m_mu_region*m_mu_region);
          }
          if( syst == SYS_MU_FR_DOWN){ // Fake Rate Down
              float etaSys   = m_mu_eta->GetBinContent( m_mu_eta->FindBin(fabs(lep.eta())) );
              etaSys = etaSys < 0.0 ? etaSys : 0.0;
              float statSys  = getRelStatError(lep, rate_type, region);
              return -1.0*sqrt(statSys*statSys);
//               return (-1.0) * sqrt( statSys*statSys
//                                     + etaSys*etaSys
//                                     + m_mu_datamc*m_mu_datamc
//                                     + m_mu_region*m_mu_region);
          }
          if( syst == SYS_MU_ETA ) return m_mu_eta->GetBinContent( m_mu_eta->FindBin(fabs(lep.eta())) );
          if( syst == SYS_MU_FR_STAT_UP )   return getRelStatError(lep, rate_type, region); // Statistical error from maps
          if( syst == SYS_MU_FR_STAT_DOWN ) return (-1) * getRelStatError(lep, rate_type, region);
          if( syst == SYS_MU_DATAMC_UP )   return m_mu_datamc; // Data/MC error
          if( syst == SYS_MU_DATAMC_DOWN ) return -m_mu_datamc;
          if( syst == SYS_MU_REG_UP )   return m_mu_region; // Error from combination
          if( syst == SYS_MU_REG_DOWN ) return -m_mu_region;
          return 0.;
      }// end muon
  }// end fake
  return 0.;
}

// -----------------------------------------------------------------------------
int DileptonMatrixMethod::DiLeptonMatrixMethod::getTT( const MatrixLepton& lep1
                                                 , const MatrixLepton& lep2
                                                 ) const
{
  if (lep1.isTight() && lep2.isTight()) return 1;
  return 0;
}

// -----------------------------------------------------------------------------
int DileptonMatrixMethod::DiLeptonMatrixMethod::getTL( const MatrixLepton& lep1
                                                 , const MatrixLepton& lep2
                                                 ) const
{
  if (lep1.isTight() && !lep2.isTight()) return 1;
  return 0;
}

// -----------------------------------------------------------------------------
int DileptonMatrixMethod::DiLeptonMatrixMethod::getLT( const MatrixLepton& lep1
                                                 , const MatrixLepton& lep2
                                                 ) const
{
  if (!lep1.isTight() && lep2.isTight()) return 1;
  return 0;
}

// -----------------------------------------------------------------------------
int DileptonMatrixMethod::DiLeptonMatrixMethod::getLL( const MatrixLepton& lep1
                                                 , const MatrixLepton& lep2
                                                 ) const
{
  if (!lep1.isTight() && !lep2.isTight()) return 1;
  return 0;
}

// -----------------------------------------------------------------------------
void DileptonMatrixMethod::DiLeptonMatrixMethod::loadSysFromFile()
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
float DileptonMatrixMethod::DiLeptonMatrixMethod::getStatError(const MatrixLepton& lep
                 , RATE_TYPE rate_type
                 , susy::fake::Region region) const

{
    float error = 0.0;
    TH1* h_rate = NULL;
    RATE_PARAM rate_param = m_rate_param_real_el; // all 4 must be the same; checked at configuration (DG todo improve this)
    if(getHistoAndParametrization(lep, region, rate_type, h_rate, rate_param)) {
        int rate_bin = getRateBin(lep, h_rate, rate_param);
        error = h_rate->GetBinError(rate_bin);
    }
  return error;
}
// ---------------------------------------------------------
float DileptonMatrixMethod::DiLeptonMatrixMethod::getRelStatError(const MatrixLepton &lep, RATE_TYPE rt, susy::fake::Region region) const
{
    float rate(0.0), error(0.0), relativeError(0.0);
    TH1* h_rate = NULL;
    RATE_PARAM rate_param = m_rate_param_real_el; // all 4 must be the same; checked at configuration (DG todo improve this)
    if(getHistoAndParametrization(lep, region, rt, h_rate, rate_param)) {
        int rate_bin = getRateBin(lep, h_rate, rate_param);
        rate = h_rate->GetBinContent(rate_bin);
        error = h_rate->GetBinError(rate_bin);
        if(rate!=0.0) relativeError = error/rate;
        else cout<<"DiLeptonMatrixMethod::getRelStatError: rate "<<rate<<" returning "<<relativeError<<endl;
    }
    return relativeError;
}
// -----------------------------------------------------------------------------
void DileptonMatrixMethod::DiLeptonMatrixMethod::printInfo(
    const MatrixLepton& lep1,
    const MatrixLepton& lep2,
    susy::fake::Region region, float MetRel,
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
// -----------------------------------------------------------------------------
int DileptonMatrixMethod::DiLeptonMatrixMethod::getIndexRegion(susy::fake::Region region)
{
    const susy::fake::Region* begin = susy::fake::SignalRegions;
    const susy::fake::Region* end   = begin + susy::fake::NumberOfSignalRegions;
    const susy::fake::Region* it    = std::find(begin, end, region);
    if(it!=end) {
      return std::distance(begin, it);
    } else {
      std::cout<<"DileptonMatrixMethod::getIndexRegion :"
               <<" error, the region '"<<region2str(region)<<"'"
               <<" is not in the SignalRegions list"<<endl;
      assert(false);
    }
}
// -----------------------------------------------------------------------------
const TArrayD* DiLeptonMatrixMethod::getPtBins() const
{
    const TArrayD *bins = 0;
    if(const TAxis *ax = getPtAxis()) bins = ax->GetXbins();
    else cout<<"DiLeptonMatrixMethod::getPtBins(): error, invalid axis"<<endl;
    return bins;
}
//----------------------------------------------------------
const TArrayD* DiLeptonMatrixMethod::getEtaBins() const
{
    const TArrayD *bins = 0;
    if(const TAxis *ax = getEtaAxis()) bins = ax->GetXbins();
    else cout<<"DiLeptonMatrixMethod::getEtaBins(): error, invalid axis"<<endl;
    return bins;
}
//----------------------------------------------------------
std::string sys2str(const DileptonMatrixMethod::SYSTEMATIC s)
{
    return DileptonMatrixMethod::systematic_names[s];
}
//----------------------------------------------------------
std::string lep2str(const DileptonMatrixMethod::MatrixLepton l)
{
    std::ostringstream oss;
    oss<<(l.isElectron() ? "el" : "mu")
       <<" "<<(l.isTight() ? "tight" : "loose")
       <<" pt "<<l.pt()
       <<" eta "<<l.eta();
    return oss.str();
}
//----------------------------------------------------------
void DiLeptonMatrixMethod::printRateSystematics(const MatrixLepton &l, RATE_TYPE &rt, susy::fake::Region &reg) const
{
    // this is a useless parameter, it should be dropped everywhere (DG, 2014-05-18 TODO)
    float dummyMetRel(20.0);
    bool isEl(l.isElectron());
    float nomRate = getRate(l, rt, reg, dummyMetRel, SYS_NOM);
    float statSys = getStatError(l, rt, reg);
    vector<SYSTEMATIC> syss;
    cout<<lep2str(l)<<": "<<endl;
    cout<<"rate : "
        <<sys2str(SYS_NOM)<<" : "<<nomRate
        <<" stat sys : "<<statSys
        <<endl;
    if(rt==REAL){
        if(isEl){
            syss.push_back(SYS_EL_RE_UP   );
            syss.push_back(SYS_EL_RE_DOWN );
        } else {
            syss.push_back(SYS_MU_RE_UP   );
            syss.push_back(SYS_MU_RE_DOWN );
        }
    } else if(rt==FAKE){
        if(isEl){
            syss.push_back(SYS_EL_FR_UP       );
            syss.push_back(SYS_EL_FR_DOWN     );
            syss.push_back(SYS_EL_HFLF_UP     );
            syss.push_back(SYS_EL_HFLF_DOWN   );
            syss.push_back(SYS_EL_ETA         );
            syss.push_back(SYS_EL_FR_STAT_UP  );
            syss.push_back(SYS_EL_FR_STAT_DOWN);
            syss.push_back(SYS_EL_DATAMC_UP   );
            syss.push_back(SYS_EL_DATAMC_DOWN );
            syss.push_back(SYS_EL_REG_UP      );
            syss.push_back(SYS_EL_REG_DOWN    );
            syss.push_back(SYS_EL_FRAC_UP     );
            syss.push_back(SYS_EL_FRAC_DO     );
            syss.push_back(SYS_MU_FRAC_UP     );
            syss.push_back(SYS_MU_FRAC_DO     );
        } else { // isEl
            syss.push_back(SYS_MU_FR_UP        );
            syss.push_back(SYS_MU_FR_DOWN      );
            syss.push_back(SYS_MU_ETA          );
            syss.push_back(SYS_MU_FR_STAT_UP   );
            syss.push_back(SYS_MU_FR_STAT_DOWN );
            syss.push_back(SYS_MU_DATAMC_UP    );
            syss.push_back(SYS_MU_DATAMC_DOWN  );
            syss.push_back(SYS_MU_REG_UP       );
            syss.push_back(SYS_MU_REG_DOWN     );
            syss.push_back(SYS_EL_FRAC_UP      );
            syss.push_back(SYS_EL_FRAC_DO      );
            syss.push_back(SYS_MU_FRAC_UP      );
            syss.push_back(SYS_MU_FRAC_DO      );
        } // end isMu
    } // end isFake
    cout<<" fractional variations : ";
    for(size_t s=0; s<syss.size(); ++s)
        cout<<sys2str(syss[s])<<" : "<<getRateSyst(l, rt, reg, dummyMetRel, syss[s])<<" ";
    cout<<endl;
}
//----------------------------------------------------------
const TH1* DiLeptonMatrixMethod::getFirstPtEtaHisto() const
{
    const TH1 *first2dhisto=0;
    if(m_rate_param_fake_el==PT_ETA){
        if(susy::fake::NumberOfSignalRegions>0) first2dhisto = m_el_fake_rate[0];
        else cout<<"DiLeptonMatrixMethod::getFirstPtEtaHisto: error, need at least one signal region"<<endl;
    } else {
        cout<<"DiLeptonMatrixMethod::getFirstPtEtaHisto can only be called with the 2d parametrization."
            <<" Returning "<<first2dhisto
            <<endl;
    }
    return first2dhisto;
}
//----------------------------------------------------------
const TAxis* DiLeptonMatrixMethod::getPtAxis() const
{
    const TAxis* ax = 0;
    if(const TH2* histo2d = static_cast<const TH2*>(getFirstPtEtaHisto())) ax = histo2d->GetXaxis();
    else cout<<"DiLeptonMatrixMethod::getPtAxis() : error, invalid histo2d, returning "<<ax<<endl;
    return ax;
}
//----------------------------------------------------------
const TAxis* DiLeptonMatrixMethod::getEtaAxis() const
{
    const TAxis* ax = 0;
    if(const TH2* histo2d = static_cast<const TH2*>(getFirstPtEtaHisto())) ax = histo2d->GetYaxis();
    else cout<<"DiLeptonMatrixMethod::getEtaAxis() : error, invalid histo2d, returning "<<ax<<endl;
    return ax;
}
//----------------------------------------------------------
float DileptonMatrixMethod::DiLeptonMatrixMethod::getFracRelativeError(const MatrixLepton &lep,
                                                                   RATE_TYPE rt,
                                                                   susy::fake::Region region,
                                                                   SYSTEMATIC syst) const
{
    float relativeError=0.0;
    if(region==susy::fake::CR_SSInc1j){
        if(rt==FAKE &&
           (syst==SYS_EL_FRAC_UP || syst==SYS_EL_FRAC_DO ||
            syst==SYS_MU_FRAC_UP || syst==SYS_MU_FRAC_DO )){
            bool isEl(lep.isElectron()), isMu(lep.isMuon());
            bool isElSys(syst==SYS_EL_FRAC_UP || syst==SYS_EL_FRAC_DO);
            bool isMuSys(syst==SYS_MU_FRAC_UP || syst==SYS_MU_FRAC_DO);
            assert(isEl!=isMu);
            int iRegion(getIndexRegion(region));
            TH1* nomHisto = isEl ? m_el_fake_rate[iRegion] : m_mu_fake_rate[iRegion];
            TH1* sysHisto = NULL;
            if((isEl && isElSys) || (isMu && isMuSys)){
                if(isEl) sysHisto = syst==SYS_EL_FRAC_UP ? m_el_frac_up : m_el_frac_do;
                if(isMu) sysHisto = syst==SYS_MU_FRAC_UP ? m_mu_frac_up : m_mu_frac_do;
                if(!nomHisto) cout<<"cannot get nom histo"<<endl;
                if(!sysHisto) cout<<"cannot get sys histo"<<endl;
                if(nomHisto&&sysHisto){
                    int nomBin = getRateBin(lep, nomHisto, PT_ETA); // TODO : assert rate_param is pt_eta [DG 2014]
                    int sysBin = getRateBin(lep, sysHisto, PT_ETA);
                    float nom = nomHisto->GetBinContent(nomBin);
                    float sys = sysHisto->GetBinContent(sysBin);
                    if(nom) relativeError=(sys-nom)/nom;
                }
            }
        } // end if(FAKE); assume real is negligible
    } else {
        cout<<"getFracRelativeError: implemented only for CR_SSInc1j, not for "<<susy::fake::RegionNames[region]
            <<"; returning "<<relativeError<<endl;
    }
    return relativeError;
}
// -----------------------------------------------------------------------------
