#include "DileptonMatrixMethod/DileptonMatrixMethod.h"

#include "TParameter.h"
#include "TVectorD.h"
#include "TArrayD.h"
#include "TAxis.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2.h"


#include <iostream>
#include <math.h>
#include <algorithm> // copy
#include <cassert>   // assert
#include <iterator>  // ostream_iterator, distance
#include <sstream>

using susy::fake::DileptonMatrixMethod;
using susy::fake::Parametrization;
using susy::fake::Systematic;

// -----------------------------------------------------------------------------
DileptonMatrixMethod::DileptonMatrixMethod():
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
DileptonMatrixMethod::~DileptonMatrixMethod()
{
  // Clean up pointers to member histograms
  std::cout << "checking m_hist_file\n";
  if (m_hist_file != NULL) {
    std::cout << "closing m_hist_file\n";
    m_hist_file->Close();
    std::cout << "deleting m_hist_file\n";
    // delete m_hist_file;
  }
  std::cout << "done with ~DileptonMatrixMethod\n";
}

// -----------------------------------------------------------------------------
/// function object to retrieve parameters
struct RetrieveOrZero {
    typedef TParameter<double>* tdp_t;
    TFile *f_;  // TFile::Get breaks const-ness
    vector<string> miss_;
    RetrieveOrZero(TFile *f) : f_(f) {}
    bool anythingMissing() const { return miss_.size()>0; }
    double operator() (const char* varname) {
        cout<<"Retrieving parameter "<<varname<<endl;
        double val(0.0);
        if(!f_) miss_.push_back(varname);
        else if(tdp_t tdp = static_cast<tdp_t>(f_->Get(varname))) val = tdp->GetVal();
        else miss_.push_back(varname);
        return val;
    }
};
// -----------------------------------------------------------------------------
/// function object to retrieve the input histos
struct RetrieveOrNull {
    TFile *f_;  // TFile::Get breaks const-ness
    vector<string> miss_;
    RetrieveOrNull(TFile *f) : f_(f) {}
    bool anythingMissing() const { return miss_.size()>0; }
    TH1* operator() (const string &histoname) {
        cout<<"Retrieving histogram "<<histoname<<endl;
        if(!f_) miss_.push_back(histoname);
        TH1* h = static_cast<TH1*>(f_->Get(histoname.c_str()));
        if(!h) miss_.push_back(histoname);
        return h;
    }
};
// -----------------------------------------------------------------------------
bool DileptonMatrixMethod::configure(const std::string &file_name,
                                     const std::vector<std::string> &region_names,
                                     Parametrization::Value real_el,
                                     Parametrization::Value fake_el,
                                     Parametrization::Value real_mu,
                                     Parametrization::Value fake_mu)
{
    bool success=false;
    bool invalidParamChoice = (real_el != fake_el ||
                               fake_el != real_mu ||
                               real_mu != fake_mu );
    m_hist_file = TFile::Open(file_name.c_str(), "OPEN");
    bool invalidInputFile = (!m_hist_file || !m_hist_file->IsOpen());
    if(invalidParamChoice){
        cout<<"DileptonMatrixMethod::configure: mixed parametrization not implemented"<<endl;
    } else if(invalidInputFile) {
        cout<<"DileptonMatrixMethod::configure: failed to open fake rate file: "<<file_name<<endl;
    } else {
        m_rate_param_real_el = real_el;
        m_rate_param_fake_el = fake_el;
        m_rate_param_real_mu = real_mu;
        m_rate_param_fake_mu = fake_mu;
        success = loadNominalFromFile(region_names);
        if(!loadSysFromFile())
            cout<<"Warning: some systematic uncertainty will not be available"<<endl;
    }
    return success;
}
// -----------------------------------------------------------------------------
bool DileptonMatrixMethod::loadNominalFromFile(const std::vector<std::string> &region_names)
{
    RetrieveOrNull retrieve(m_hist_file);
    const Parametrization::Value rp = m_rate_param_real_el;
    std::string s_el_real("el_real_eff_"), s_el_fake("el_fake_rate_");
    std::string s_mu_real("mu_real_eff_"), s_mu_fake("mu_fake_rate_");
    if(rp==Parametrization::PT_ETA) { // todo: use string replacement (eff->eff2d, rate->rate2d)
        s_el_real = "el_real_eff2d_"; s_el_fake = "el_fake_rate2d_";
        s_mu_real = "mu_real_eff2d_"; s_mu_fake = "mu_fake_rate2d_";
    }
    for(size_t r=0; r<region_names.size(); ++r) { // Get the Real eff and Fake rate for electrons and muons
        const std::string &regionName = region_names[r];
        TH1* el_real_eff   = retrieve(s_el_real + regionName);
        TH1* el_fake_rate  = retrieve(s_el_fake + regionName);
        TH1* mu_real_eff   = retrieve(s_mu_real + regionName);
        TH1* mu_fake_rate  = retrieve(s_mu_fake + regionName);
        bool allHistosFound = (el_real_eff && el_fake_rate && mu_real_eff && mu_fake_rate);
        if(allHistosFound){
            m_signalRegions.push_back(regionName);
            m_el_real_eff .push_back(el_real_eff );
            m_el_fake_rate.push_back(el_fake_rate);
            m_mu_real_eff .push_back(mu_real_eff );
            m_mu_fake_rate.push_back(mu_fake_rate);
        } else {
            cout<<"DileptonMatrixMethod: missing nominal histos for '"<<regionName<<"': skipping this region"<<endl;
        }
    }
    if(retrieve.anythingMissing()){
        cout<<"DileptonMatrixMethod::loadNominalFromFile(): missing the following histograms:"<<endl;
        copy(retrieve.miss_.begin(), retrieve.miss_.end(), ostream_iterator<string>(cout,"\n\t"));
    }
    bool success = region_names.size()==m_signalRegions.size();
    return success;
}
//-----------------------------------------------------------------------------
float DileptonMatrixMethod::getTotalFake(
    bool isTight1, bool isElectron1, float pt1, float eta1,
    bool isTight2, bool isElectron2, float pt2, float eta2,
    size_t regionIndex,
    float MetRel,
    Systematic::Value syst) const
{
  susy::fake::MatrixLepton lep1(isTight1, isElectron1, pt1, eta1);
  susy::fake::MatrixLepton lep2(isTight2, isElectron2, pt2, eta2);

  return getTotalFake(lep1, lep2, regionIndex, MetRel, syst);
}

// -----------------------------------------------------------------------------
float DileptonMatrixMethod::getRR(
    bool isTight1, bool isElectron1, float pt1, float eta1,
    bool isTight2, bool isElectron2, float pt2, float eta2,
    size_t regionIndex,
    float MetRel,
    Systematic::Value syst) const
{
  susy::fake::MatrixLepton lep1(isTight1, isElectron1, pt1, eta1);
  susy::fake::MatrixLepton lep2(isTight2, isElectron2, pt2, eta2);

  return getRR(lep1, lep2, regionIndex, MetRel, syst);
}

// -----------------------------------------------------------------------------
float DileptonMatrixMethod::getRF(
    bool isTight1, bool isElectron1, float pt1, float eta1,
    bool isTight2, bool isElectron2, float pt2, float eta2,
    size_t regionIndex,
    float MetRel,
    Systematic::Value syst) const
{
  susy::fake::MatrixLepton lep1(isTight1, isElectron1, pt1, eta1);
  susy::fake::MatrixLepton lep2(isTight2, isElectron2, pt2, eta2);

  return getRF(lep1, lep2, regionIndex, MetRel, syst);
}

// -----------------------------------------------------------------------------
float DileptonMatrixMethod::getFR(
    bool isTight1, bool isElectron1, float pt1, float eta1,
    bool isTight2, bool isElectron2, float pt2, float eta2,
    size_t regionIndex,
    float MetRel,
    Systematic::Value syst) const
{
  susy::fake::MatrixLepton lep1(isTight1, isElectron1, pt1, eta1);
  susy::fake::MatrixLepton lep2(isTight2, isElectron2, pt2, eta2);

  return getFR(lep1, lep2, regionIndex, MetRel, syst);
}

// -----------------------------------------------------------------------------
float DileptonMatrixMethod::getFF(
    bool isTight1, bool isElectron1, float pt1, float eta1,
    bool isTight2, bool isElectron2, float pt2, float eta2,
    size_t regionIndex,
    float MetRel,
    Systematic::Value syst) const
{
  susy::fake::MatrixLepton lep1(isTight1, isElectron1, pt1, eta1);
  susy::fake::MatrixLepton lep2(isTight2, isElectron2, pt2, eta2);

  return getFF(lep1, lep2, regionIndex, MetRel, syst);
}

// -----------------------------------------------------------------------------
float DileptonMatrixMethod::getTotalFake(
    const MatrixLepton& lep1,
    const MatrixLepton& lep2,
    size_t regionIndex,
    float MetRel,
    Systematic::Value syst) const
{
  float fake_contribution = getRF(lep1, lep2, regionIndex, MetRel, syst);
  fake_contribution += getFR(lep1, lep2, regionIndex, MetRel, syst);
  fake_contribution += getFF(lep1, lep2, regionIndex, MetRel, syst);
  return fake_contribution;
}

// -----------------------------------------------------------------------------
float DileptonMatrixMethod::getRR(
    const MatrixLepton& lep1,
    const MatrixLepton& lep2,
    size_t regionIndex,
    float MetRel,
    Systematic::Value syst) const
{
  int n_tt = getTT(lep1, lep2);
  int n_tl = getTL(lep1, lep2);
  int n_lt = getLT(lep1, lep2);
  int n_ll = getLL(lep1, lep2);

  float r1 = getRate(lep1, REAL, regionIndex, MetRel, syst);
  float r2 = getRate(lep2, REAL, regionIndex, MetRel, syst);

  float f1 = getRate(lep1, FAKE, regionIndex, MetRel, syst);
  float f2 = getRate(lep2, FAKE, regionIndex, MetRel, syst);
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
float DileptonMatrixMethod::getRF(
    const MatrixLepton& lep1,
    const MatrixLepton& lep2,
    size_t regionIndex,
    float MetRel,
    Systematic::Value syst) const
{
  int n_tt = getTT(lep1, lep2);
  int n_tl = getTL(lep1, lep2);
  int n_lt = getLT(lep1, lep2);
  int n_ll = getLL(lep1, lep2);

  float r1 = getRate(lep1, REAL, regionIndex, MetRel, syst);
  float r2 = getRate(lep2, REAL, regionIndex, MetRel, syst);

  float f1 = getRate(lep1, FAKE, regionIndex, MetRel, syst);
  float f2 = getRate(lep2, FAKE, regionIndex, MetRel, syst);

  float rf = r1*f2/(r1-f1)/(r2-f2)*( (f1-1.)*(1.-r2)*n_tt
                                   + (1.-f1)*r2*n_tl
                                   + f1*(1.-r2)*n_lt
                                   - f1*r2*n_ll
                                   );
  return rf;
}

// -----------------------------------------------------------------------------
float DileptonMatrixMethod::getFR(
    const MatrixLepton& lep1,
    const MatrixLepton& lep2,
    size_t regionIndex,
    float MetRel,
    Systematic::Value syst) const
{
  int n_tt = getTT(lep1, lep2);
  int n_tl = getTL(lep1, lep2);
  int n_lt = getLT(lep1, lep2);
  int n_ll = getLL(lep1, lep2);

  float r1 = getRate(lep1, REAL, regionIndex, MetRel, syst);
  float r2 = getRate(lep2, REAL, regionIndex, MetRel, syst);

  float f1 = getRate(lep1, FAKE, regionIndex, MetRel, syst);
  float f2 = getRate(lep2, FAKE, regionIndex, MetRel, syst);

  float fr = f1*r2/(r1-f1)/(r2-f2)*( (r1-1.)*(1.-f2)*n_tt
                                   + (1.-r1)*f2*n_tl
                                   + r1*(1.-f2)*n_lt
                                   - r1*f2*n_ll
                                   );
  return fr;
}

// -----------------------------------------------------------------------------
float DileptonMatrixMethod::getFF(
    const MatrixLepton& lep1,
    const MatrixLepton& lep2,
    size_t regionIndex,
    float MetRel,
    Systematic::Value syst) const
{
  int n_tt = getTT(lep1, lep2);
  int n_tl = getTL(lep1, lep2);
  int n_lt = getLT(lep1, lep2);
  int n_ll = getLL(lep1, lep2);

  float r1 = getRate(lep1, REAL, regionIndex, MetRel, syst);
  float r2 = getRate(lep2, REAL, regionIndex, MetRel, syst);

  float f1 = getRate(lep1, FAKE, regionIndex, MetRel, syst);
  float f2 = getRate(lep2, FAKE, regionIndex, MetRel, syst);

  float ff = f1*f2/(r1-f1)/(r2-f2)*( (1.-r1)*(1.-r2)*n_tt
                                   + (r1-1.)*r2*n_tl
                                   + r1*(r2-1.)*n_lt
                                   + r1*r2*n_ll
                                   );

  return ff;
}

// -----------------------------------------------------------------------------
float DileptonMatrixMethod::getRate(
    bool isTight, bool isElectron, float pt, float eta,
    RATE_TYPE rate_type, size_t regionIndex, float MetRel,
    Systematic::Value syst) const
{
  susy::fake::MatrixLepton lep(isTight, isElectron, pt, eta);
  return getRate(lep, rate_type, regionIndex, MetRel, syst);
}

// --------------------------------------------------------
bool DileptonMatrixMethod::getHistoAndParametrization(const MatrixLepton &lep,
                                                      const size_t regionIndex,
                                                      const RATE_TYPE &rt,
                                                      TH1* &h, Parametrization::Value &rp) const
{
    bool found=false;
    if(regionIndex<m_signalRegions.size()){
        if(rt==REAL) {
            if (lep.isElectron()) {
                h  = m_el_real_eff[regionIndex];
                rp = m_rate_param_real_el;
                found = true;
            } else {
                h  = m_mu_real_eff[regionIndex];
                rp = m_rate_param_real_mu;
                found = true;
            }
        } else if(rt==FAKE) {
            if (lep.isElectron()) {
                h  = m_el_fake_rate[regionIndex];
                rp = m_rate_param_fake_el;
                found = true;
            } else {
                h  = m_mu_fake_rate[regionIndex];
                rp = m_rate_param_fake_mu;
                found = true;
            }
        }
    }
    return found;
}
// -----------------------------------------------------------------------------
float DileptonMatrixMethod::getRate(
    const MatrixLepton& lep,
    RATE_TYPE rate_type,
    const size_t regionIndex,
    float MetRel,
    Systematic::Value syst) const
{
    float rate = 0.0;
    TH1* h_rate = NULL;
    Parametrization::Value rate_param = m_rate_param_real_el; // all 4 must be the same; checked at configuration (DG todo improve this)
    if(getHistoAndParametrization(lep, regionIndex, rate_type, h_rate, rate_param)) {
        int rate_bin = getRateBin(lep, h_rate, rate_param);
        rate = h_rate->GetBinContent(rate_bin);
        rate *= (1.0+getRateSyst(lep, rate_type, regionIndex, MetRel, syst));
        if( rate > 1.0 ) rate = 1.0;
        if( rate < 0.0 ) rate = 1e-5;
    } else {
        cout<<"DileptonMatrixMethod::getRate: cannot get histo"<<endl;
    }
  return rate;
}

// -----------------------------------------------------------------------------
int DileptonMatrixMethod::getRateBin( const MatrixLepton& lep,
							TH1* h_rate,
							Parametrization::Value rate_param) const
{
  // Handle 2-D param
    if( rate_param == Parametrization::PT_ETA ){
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
    if( rate_param == Parametrization::PT ){
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
float DileptonMatrixMethod::getRateSyst(
    const MatrixLepton& lep,
    RATE_TYPE rate_type,
    const size_t regionIndex,
    float MetRel,
    Systematic::Value syst) const
{
  if( syst == Systematic::SYS_NOM ) return 0.0;
  if(!Systematic::isValid(syst)) { std::cout << "WARNING: invalid SYSTEMATIC type\n"; }
  if( rate_type == REAL ){ // Real Eff Sys
      float statSys = getRelStatError(lep, rate_type, regionIndex);
      if( lep.isElectron() ){
          if( syst == Systematic::SYS_EL_RE_UP )   return sqrt( m_el_real_up*m_el_real_up + statSys*statSys);
          if( syst == Systematic::SYS_EL_RE_DOWN ) return -sqrt( m_el_real_down*m_el_real_down + statSys*statSys);
      } else {
          if( syst == Systematic::SYS_MU_RE_UP )   return sqrt( m_mu_real_up*m_mu_real_up + statSys*statSys);
          if( syst == Systematic::SYS_MU_RE_DOWN ) return -sqrt( m_mu_real_down*m_mu_real_down + statSys*statSys);
    }
      return 0.;
  }
  if( rate_type == FAKE ){ // Fake Rate Sys
      if(syst==Systematic::SYS_EL_FRAC_UP || syst==Systematic::SYS_EL_FRAC_DO ||
         syst==Systematic::SYS_MU_FRAC_UP || syst==Systematic::SYS_MU_FRAC_DO )
          return getFracRelativeError(lep, rate_type, regionIndex, syst);
      if( lep.isElectron() ){
          if( syst == Systematic::SYS_EL_FR_UP ){ // Fake Rate Up
              float etaSys   = m_el_eta->GetBinContent( m_el_eta->FindBin(fabs(lep.eta())) );
              etaSys = etaSys > 0 ? etaSys : 0.;
              float statSys = getRelStatError(lep, rate_type, regionIndex);
              return +1.0*sqrt(statSys*statSys);
//               return sqrt( m_el_HFLFerr*m_el_HFLFerr //add in quadrature
//                            + statSys*statSys
//                            + etaSys*etaSys
//                            + m_el_datamc*m_el_datamc
//                            + m_el_region*m_el_region);
          }
          if( syst == Systematic::SYS_EL_FR_DOWN){ // Fake Rate Down
              float etaSys   = m_el_eta->GetBinContent( m_el_eta->FindBin(fabs(lep.eta())) );
              etaSys = etaSys < 0.0 ? etaSys : 0.0;
              float statSys  = getRelStatError(lep, rate_type, regionIndex);
              return -1.0*sqrt(statSys*statSys);
//               return (-1.0) * sqrt(m_el_HFLFerr*m_el_HFLFerr // add in quadrature
//                                    + statSys*statSys
//                                    + etaSys*etaSys
//                                    + m_el_datamc*m_el_datamc
//                                    + m_el_region*m_el_region);
          }
          if( syst == Systematic::SYS_EL_HFLF_UP )   return +m_el_HFLFerr; // HF/LF error
          if( syst == Systematic::SYS_EL_HFLF_DOWN ) return -m_el_HFLFerr;
          if( syst == Systematic::SYS_EL_ETA ) return m_el_eta->GetBinContent( m_el_eta->FindBin(fabs(lep.eta())) );
          if( syst == Systematic::SYS_EL_FR_STAT_UP )   return +1.0*getRelStatError(lep, rate_type, regionIndex); // Statistical error from maps
          if( syst == Systematic::SYS_EL_FR_STAT_DOWN ) return -1.0*getRelStatError(lep, rate_type, regionIndex);
          if( syst == Systematic::SYS_EL_DATAMC_UP )   return +m_el_datamc; // Data/MC error
          if( syst == Systematic::SYS_EL_DATAMC_DOWN ) return -m_el_datamc;
          if( syst == Systematic::SYS_EL_REG_UP )   return +m_el_region; // Error from combinatino
          if( syst == Systematic::SYS_EL_REG_DOWN ) return -m_el_region;
          return 0.;
      } else { // end electron
          if( syst == Systematic::SYS_MU_FR_UP ){ // Fake Rate Up
              float etaSys   = m_mu_eta->GetBinContent( m_mu_eta->FindBin(fabs(lep.eta())) );
              etaSys = etaSys > 0.0 ? etaSys : 0.0;
              float statSys  = getRelStatError(lep, rate_type, regionIndex);
              return +1.0*sqrt(statSys*statSys);
//               return sqrt( statSys*statSys
//                            + etaSys*etaSys
//                            + m_mu_datamc*m_mu_datamc
//                            + m_mu_region*m_mu_region);
          }
          if( syst == Systematic::SYS_MU_FR_DOWN){ // Fake Rate Down
              float etaSys   = m_mu_eta->GetBinContent( m_mu_eta->FindBin(fabs(lep.eta())) );
              etaSys = etaSys < 0.0 ? etaSys : 0.0;
              float statSys  = getRelStatError(lep, rate_type, regionIndex);
              return -1.0*sqrt(statSys*statSys);
//               return (-1.0) * sqrt( statSys*statSys
//                                     + etaSys*etaSys
//                                     + m_mu_datamc*m_mu_datamc
//                                     + m_mu_region*m_mu_region);
          }
          if( syst == Systematic::SYS_MU_ETA ) return m_mu_eta->GetBinContent( m_mu_eta->FindBin(fabs(lep.eta())) );
          if( syst == Systematic::SYS_MU_FR_STAT_UP )   return getRelStatError(lep, rate_type, regionIndex); // Statistical error from maps
          if( syst == Systematic::SYS_MU_FR_STAT_DOWN ) return (-1) * getRelStatError(lep, rate_type, regionIndex);
          if( syst == Systematic::SYS_MU_DATAMC_UP )   return m_mu_datamc; // Data/MC error
          if( syst == Systematic::SYS_MU_DATAMC_DOWN ) return -m_mu_datamc;
          if( syst == Systematic::SYS_MU_REG_UP )   return m_mu_region; // Error from combination
          if( syst == Systematic::SYS_MU_REG_DOWN ) return -m_mu_region;
          return 0.;
      }// end muon
  }// end fake
  return 0.;
}

// -----------------------------------------------------------------------------
int DileptonMatrixMethod::getTT( const MatrixLepton& lep1
                                                 , const MatrixLepton& lep2
                                                 ) const
{
  if (lep1.isTight() && lep2.isTight()) return 1;
  return 0;
}

// -----------------------------------------------------------------------------
int DileptonMatrixMethod::getTL( const MatrixLepton& lep1
                                                 , const MatrixLepton& lep2
                                                 ) const
{
  if (lep1.isTight() && !lep2.isTight()) return 1;
  return 0;
}

// -----------------------------------------------------------------------------
int DileptonMatrixMethod::getLT( const MatrixLepton& lep1
                                                 , const MatrixLepton& lep2
                                                 ) const
{
  if (!lep1.isTight() && lep2.isTight()) return 1;
  return 0;
}

// -----------------------------------------------------------------------------
int DileptonMatrixMethod::getLL( const MatrixLepton& lep1
                                                 , const MatrixLepton& lep2
                                                 ) const
{
  if (!lep1.isTight() && !lep2.isTight()) return 1;
  return 0;
}

// -----------------------------------------------------------------------------
bool DileptonMatrixMethod::loadSysFromFile()
{
  RetrieveOrZero roz(m_hist_file);
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
    cout<<"DileptonMatrixMethod::loadSysFromFile(): missing the following systematics:"<<endl;
    copy(roz.miss_.begin(), roz.miss_.end(), ostream_iterator<string>(cout,"\n\t"));
    cout<<"Will use 0.0 for these"<<endl;
  }
  RetrieveOrNull ron(m_hist_file);
  m_el_frac_up = ron("el_fake_rate2d_CR_SSInc1j_frac_up");
  m_el_frac_do = ron("el_fake_rate2d_CR_SSInc1j_frac_do");
  m_mu_frac_up = ron("mu_fake_rate2d_CR_SSInc1j_frac_up");
  m_mu_frac_do = ron("mu_fake_rate2d_CR_SSInc1j_frac_do");
  m_el_metrel  = ron("el_metrel_sys");
  m_mu_metrel  = ron("mu_metrel_sys");
  m_el_eta     = ron("el_eta_sys");
  m_mu_eta     = ron("mu_eta_sys");
  if(ron.anythingMissing()) {
      cout<<"DileptonMatrixMethod::loadSysFromFile(): missing the following systematic histograms:"<<endl;
      copy(ron.miss_.begin(), ron.miss_.end(), ostream_iterator<string>(cout,"\n\t"));
  }  
  return (!roz.anythingMissing() && !ron.anythingMissing());
}
// -----------------------------------------------------------------------------
float DileptonMatrixMethod::getStatError(const MatrixLepton& lep
                                         , RATE_TYPE rate_type
                                         , const size_t regionIndex) const

{
    float error = 0.0;
    TH1* h_rate = NULL;
    Parametrization::Value rate_param = m_rate_param_real_el; // all 4 must be the same; checked at configuration (DG todo improve this)
    if(getHistoAndParametrization(lep, regionIndex, rate_type, h_rate, rate_param)) {
        int rate_bin = getRateBin(lep, h_rate, rate_param);
        error = h_rate->GetBinError(rate_bin);
    }
  return error;
}
// ---------------------------------------------------------
float DileptonMatrixMethod::getRelStatError(const MatrixLepton &lep, RATE_TYPE rt, 
                                            const size_t regionIndex) const
{
    float rate(0.0), error(0.0), relativeError(0.0);
    TH1* h_rate = NULL;
    Parametrization::Value rate_param = m_rate_param_real_el; // all 4 must be the same; checked at configuration (DG todo improve this)
    if(getHistoAndParametrization(lep, regionIndex, rt, h_rate, rate_param)) {
        int rate_bin = getRateBin(lep, h_rate, rate_param);
        rate = h_rate->GetBinContent(rate_bin);
        error = h_rate->GetBinError(rate_bin);
        if(rate!=0.0) relativeError = error/rate;
        else cout<<"DileptonMatrixMethod::getRelStatError: rate "<<rate<<" returning "<<relativeError<<endl;
    }
    return relativeError;
}
// -----------------------------------------------------------------------------
void DileptonMatrixMethod::printInfo(
    const MatrixLepton& lep1,
    const MatrixLepton& lep2,
    const size_t regionIndex,
    float MetRel,
    Systematic::Value syst) const
{
  float r1 = getRate(lep1, REAL, regionIndex, MetRel, syst);
  float r2 = getRate(lep2, REAL, regionIndex, MetRel, syst);

  float f1 = getRate(lep1, FAKE, regionIndex, MetRel, syst);
  float f2 = getRate(lep2, FAKE, regionIndex, MetRel, syst);

  int num_electrons = 0;
  if (lep1.isElectron()) ++num_electrons;
  if (lep2.isElectron()) ++num_electrons;

  std::string event_type = "";
  if (num_electrons == 2) event_type = "ee";
  else if (num_electrons == 1) event_type = "emu";
  else if (num_electrons == 0) event_type = "mumu";

  std::cout << "-----------------------------------------------\n";
  std::cout << "region: " << m_signalRegions[regionIndex] << "\t-\tsyst: " << syst << "\n";
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
  std::cout << "rr: " << getRR(lep1,lep2,regionIndex,syst) << "\n";
  std::cout << "rf: " << getRF(lep1,lep2,regionIndex,syst) << "\n";
  std::cout << "fr: " << getFR(lep1,lep2,regionIndex,syst) << "\n";
  std::cout << "ff: " << getFF(lep1,lep2,regionIndex,syst) << "\n";
  std::cout << "total fake: " << getTotalFake(lep1,lep2,regionIndex,syst) << "\n";
}
// -----------------------------------------------------------------------------
size_t DileptonMatrixMethod::getIndexRegion(const std::string &regionName) const
{
    typedef std::vector<std::string>::const_iterator vsit_t;
    vsit_t begin = m_signalRegions.begin();
    vsit_t end   = m_signalRegions.end();
    vsit_t it    = std::find(begin, end, regionName);
    if(it!=end) {
      return std::distance(begin, it);
    } else {
      std::cout<<"DileptonMatrixMethod::getIndexRegion :"
               <<" error, the region '"<<regionName<<"'"
               <<" is not in the SignalRegions list"<<endl;
      assert(false);
    }
}
// -----------------------------------------------------------------------------
const TArrayD* DileptonMatrixMethod::getPtBins() const
{
    const TArrayD *bins = 0;
    if(const TAxis *ax = getPtAxis()) bins = ax->GetXbins();
    else cout<<"DileptonMatrixMethod::getPtBins(): error, invalid axis"<<endl;
    return bins;
}
//----------------------------------------------------------
const TArrayD* DileptonMatrixMethod::getEtaBins() const
{
    const TArrayD *bins = 0;
    if(const TAxis *ax = getEtaAxis()) bins = ax->GetXbins();
    else cout<<"DileptonMatrixMethod::getEtaBins(): error, invalid axis"<<endl;
    return bins;
}
//----------------------------------------------------------
std::string lep2str(const susy::fake::MatrixLepton l)
{
    std::ostringstream oss;
    oss<<(l.isElectron() ? "el" : "mu")
       <<" "<<(l.isTight() ? "tight" : "loose")
       <<" pt "<<l.pt()
       <<" eta "<<l.eta();
    return oss.str();
}
//----------------------------------------------------------
void DileptonMatrixMethod::printRateSystematics(const MatrixLepton &l, RATE_TYPE &rt, size_t regionIndex) const
{
    // this is a useless parameter, it should be dropped everywhere (DG, 2014-05-18 TODO)
    float dummyMetRel(20.0);
    bool isEl(l.isElectron());
    float nomRate = getRate(l, rt, regionIndex, dummyMetRel, Systematic::SYS_NOM);
    float statSys = getStatError(l, rt, regionIndex);
    vector<Systematic::Value> syss;
    cout<<lep2str(l)<<": "<<endl;
    cout<<"rate : "
        <<Systematic::str(Systematic::SYS_NOM)<<" : "<<nomRate
        <<" stat sys : "<<statSys
        <<endl;
    if(rt==REAL){
        if(isEl){
            syss.push_back(Systematic::SYS_EL_RE_UP   );
            syss.push_back(Systematic::SYS_EL_RE_DOWN );
        } else {
            syss.push_back(Systematic::SYS_MU_RE_UP   );
            syss.push_back(Systematic::SYS_MU_RE_DOWN );
        }
    } else if(rt==FAKE){
        if(isEl){
            syss.push_back(Systematic::SYS_EL_FR_UP       );
            syss.push_back(Systematic::SYS_EL_FR_DOWN     );
            syss.push_back(Systematic::SYS_EL_HFLF_UP     );
            syss.push_back(Systematic::SYS_EL_HFLF_DOWN   );
            syss.push_back(Systematic::SYS_EL_ETA         );
            syss.push_back(Systematic::SYS_EL_FR_STAT_UP  );
            syss.push_back(Systematic::SYS_EL_FR_STAT_DOWN);
            syss.push_back(Systematic::SYS_EL_DATAMC_UP   );
            syss.push_back(Systematic::SYS_EL_DATAMC_DOWN );
            syss.push_back(Systematic::SYS_EL_REG_UP      );
            syss.push_back(Systematic::SYS_EL_REG_DOWN    );
            syss.push_back(Systematic::SYS_EL_FRAC_UP     );
            syss.push_back(Systematic::SYS_EL_FRAC_DO     );
            syss.push_back(Systematic::SYS_MU_FRAC_UP     );
            syss.push_back(Systematic::SYS_MU_FRAC_DO     );
        } else { // isEl
            syss.push_back(Systematic::SYS_MU_FR_UP        );
            syss.push_back(Systematic::SYS_MU_FR_DOWN      );
            syss.push_back(Systematic::SYS_MU_ETA          );
            syss.push_back(Systematic::SYS_MU_FR_STAT_UP   );
            syss.push_back(Systematic::SYS_MU_FR_STAT_DOWN );
            syss.push_back(Systematic::SYS_MU_DATAMC_UP    );
            syss.push_back(Systematic::SYS_MU_DATAMC_DOWN  );
            syss.push_back(Systematic::SYS_MU_REG_UP       );
            syss.push_back(Systematic::SYS_MU_REG_DOWN     );
            syss.push_back(Systematic::SYS_EL_FRAC_UP      );
            syss.push_back(Systematic::SYS_EL_FRAC_DO      );
            syss.push_back(Systematic::SYS_MU_FRAC_UP      );
            syss.push_back(Systematic::SYS_MU_FRAC_DO      );
        } // end isMu
    } // end isFake
    cout<<" fractional variations : ";
    for(size_t s=0; s<syss.size(); ++s)
        cout<<Systematic::str(syss[s])<<" : "<<getRateSyst(l, rt, regionIndex, dummyMetRel, syss[s])<<" ";
    cout<<endl;
}
//----------------------------------------------------------
const TH1* DileptonMatrixMethod::getFirstPtEtaHisto() const
{
    const TH1 *first2dhisto=0;
    if(m_rate_param_fake_el==Parametrization::PT_ETA){
        if(m_el_fake_rate.size()>0 && m_el_fake_rate[0])
            first2dhisto = m_el_fake_rate[0];
        else
            cout<<"DileptonMatrixMethod::getFirstPtEtaHisto: error, need at least one signal region"<<endl;
    } else {
        cout<<"DileptonMatrixMethod::getFirstPtEtaHisto can only be called with the 2d parametrization."
            <<" Returning "<<first2dhisto
            <<endl;
    }
    return first2dhisto;
}
//----------------------------------------------------------
const TAxis* DileptonMatrixMethod::getPtAxis() const
{
    const TAxis* ax = 0;
    if(const TH2* histo2d = static_cast<const TH2*>(getFirstPtEtaHisto())) ax = histo2d->GetXaxis();
    else cout<<"DileptonMatrixMethod::getPtAxis() : error, invalid histo2d, returning "<<ax<<endl;
    return ax;
}
//----------------------------------------------------------
const TAxis* DileptonMatrixMethod::getEtaAxis() const
{
    const TAxis* ax = 0;
    if(const TH2* histo2d = static_cast<const TH2*>(getFirstPtEtaHisto())) ax = histo2d->GetYaxis();
    else cout<<"DileptonMatrixMethod::getEtaAxis() : error, invalid histo2d, returning "<<ax<<endl;
    return ax;
}
//----------------------------------------------------------
float DileptonMatrixMethod::getFracRelativeError(const MatrixLepton &lep,
                                                 RATE_TYPE rt,
                                                 size_t regionIndex,
                                                 susy::fake::Systematic::Value syst) const
{
    float relativeError=0.0;
    const std::string &region = m_signalRegions[regionIndex];
    bool validRegion = (region=="CR_SSInc1j");
    bool validParam = (m_rate_param_real_el==Parametrization::PT_ETA);
    if(validRegion && validParam){
        if(rt==FAKE &&
           (syst==Systematic::SYS_EL_FRAC_UP || syst==Systematic::SYS_EL_FRAC_DO ||
            syst==Systematic::SYS_MU_FRAC_UP || syst==Systematic::SYS_MU_FRAC_DO )){
            bool isEl(lep.isElectron()), isMu(lep.isMuon());
            bool isElSys(syst==Systematic::SYS_EL_FRAC_UP || syst==Systematic::SYS_EL_FRAC_DO);
            bool isMuSys(syst==Systematic::SYS_MU_FRAC_UP || syst==Systematic::SYS_MU_FRAC_DO);
            assert(isEl!=isMu);
            int iRegion(getIndexRegion(region));
            TH1* nomHisto = isEl ? m_el_fake_rate[iRegion] : m_mu_fake_rate[iRegion];
            TH1* sysHisto = NULL;
            if((isEl && isElSys) || (isMu && isMuSys)){
                if(isEl) sysHisto = syst==Systematic::SYS_EL_FRAC_UP ? m_el_frac_up : m_el_frac_do;
                if(isMu) sysHisto = syst==Systematic::SYS_MU_FRAC_UP ? m_mu_frac_up : m_mu_frac_do;
                if(!nomHisto) cout<<"cannot get nom histo"<<endl;
                if(!sysHisto) cout<<"cannot get sys histo"<<endl;
                if(nomHisto&&sysHisto){
                    int nomBin = getRateBin(lep, nomHisto, Parametrization::PT_ETA);
                    int sysBin = getRateBin(lep, sysHisto, Parametrization::PT_ETA);
                    float nom = nomHisto->GetBinContent(nomBin);
                    float sys = sysHisto->GetBinContent(sysBin);
                    if(nom) relativeError=(sys-nom)/nom;
                }
            }
        } // end if(FAKE); assume real is negligible
    } else {
        cout<<"getFracRelativeError: implemented only for CR_SSInc1j 2D param,"
            <<"not for '"<<region<<"'; returning "<<relativeError<<endl;
    }
    return relativeError;
}
// -----------------------------------------------------------------------------
