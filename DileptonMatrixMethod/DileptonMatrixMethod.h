// -*- c++ -*-
#ifndef SUSY_FAKE_DILEPTONMATRIXMETHOD_H
#define SUSY_FAKE_DILEPTONMATRIXMETHOD_H

#include "DileptonMatrixMethod/Lepton.h"
#include "DileptonMatrixMethod/Parametrization.h"
#include "DileptonMatrixMethod/Systematic.h"

#include <vector>
#include <string>

class TFile;
class TH1;
class TH1F;
class TArrayD;
class TAxis;

namespace susy{
namespace fake{
/**
   Class to perform the matrix method fake estimate.

   All physical quantities [E] are expressed in GeV.
   The format of the names expected for the input histograms is the
   one specified in generateNominalHistoname()

   Cleanup and rewrite, davide.gerbaudo@gmail.com, July 2014.

   @author Brett Jackson <Brett.David.Jackson@cern.ch>
   @data July 2012
*/


class DileptonMatrixMethod
  {
  public:
      /// internally we store p(T|L, Real) and p(T|L, Fake)
      enum RATE_TYPE { REAL, FAKE };

    public:
      DileptonMatrixMethod();
      ~DileptonMatrixMethod();
      /// configure the tool with a histogram file
      /**
         Find the necessary histograms for the requested
         regions, and store them in the internal pointers. If you are
         requesting signal regions for which the histograms are not
         there in the input file, the tool will warn you.
      */
      bool configure(const std::string &file_name,
                     const std::vector<std::string> &region_names,
                     Parametrization::Value real_el,
                     Parametrization::Value fake_el,
                     Parametrization::Value real_mu,
                     Parametrization::Value fake_mu);
      float getRR(const Lepton& lep1, const Lepton& lep2,
                  size_t regionIndex, float MetRel,
                  Systematic::Value syst = Systematic::SYS_NOM) const;
      float getRF(const Lepton& lep1, const Lepton& lep2,
                  size_t regionIndex, float MetRel,
                  Systematic::Value syst = Systematic::SYS_NOM) const;
      float getFR(const Lepton& lep1, const Lepton& lep2,
                  size_t regionIndex, float MetRel,
                  Systematic::Value syst = Systematic::SYS_NOM) const;
      float getFF(const Lepton& lep1, const Lepton& lep2,
                  size_t regionIndex, float MetRel,
                  Systematic::Value syst = Systematic::SYS_NOM) const;
      /// Get the total fake contribution for this event
      float getTotalFake(const Lepton& l1, const Lepton& l2, const size_t regionIndex,
                         float MetRel, Systematic::Value syst = Systematic::SYS_NOM) const;
      /// given a region, determine the internal index used to store its histograms; abort if invalid
      size_t getIndexRegion(const std::string &regionName) const;
      /// names of the available signal regions
      const std::vector<std::string>& signalRegions() const { return m_signalRegions; }

      const TArrayD* getPtBins() const;
      const TArrayD* getEtaBins() const;
      void printRateSystematics(const Lepton &l, RATE_TYPE &rt, size_t regionIndex) const;
      ///< name of the nominal histogram for a given region, parametrization, and lepton type
      /**
         Either el or mu; either real or fake.
         Histogram names for the systematic variations are slighly
         different. They are currently hardcoded in
         DileptonMatrixMethod::loadSysFromFile()
       */
      static std::string generateNominalHistoname(bool isEl, bool isReal,
                                                  const Parametrization::Value &p,
                                                  const std::string &region);
  protected:
      // Get the rate for this lepton -- real/fake for electron or muon
      // Specific rate depends on type of lepton supplied and the
      // RATE_TYPE parameter
      // for internal use
      float getRate(const Lepton&,
                    RATE_TYPE,
                    size_t regionIndex,
                    float MetRel,
                    Systematic::Value syst = Systematic::SYS_NOM) const;
      float getRateSyst(const Lepton&,
                        RATE_TYPE,
                        size_t regionIndex,
                        float MetRel,
                        Systematic::Value syst = Systematic::SYS_NOM) const;

      int getRateBin(const Lepton& lep,
                     TH1* h_rate,
                     Parametrization::Value rate_param) const;
      // Additional methods needed for systematics
      float getStatError(const Lepton&, RATE_TYPE, size_t regionIndex) const;
      float getRelStatError(const Lepton&, RATE_TYPE, size_t regionIndex) const;

      // Determine if the event is TT/TL/LT/LL -- for internal use
      int getTT(const Lepton& lep1, const Lepton& lep2) const;
      int getTL(const Lepton& lep1, const Lepton& lep2) const;
      int getLT(const Lepton& lep1, const Lepton& lep2) const;
      int getLL(const Lepton& lep1, const Lepton& lep2) const;

      // Verbose output
      void printInfo(const Lepton& lep1,
                     const Lepton& lep2,
                     size_t regionIndex,
                     float MetRel,
                     Systematic::Value syst = Systematic::SYS_NOM) const;
      /// retrieve nominal histos and paramerers
      bool loadNominalFromFile(const std::vector<std::string> &region_names);
      /// Get systematic value from file -- for internal use
      bool loadSysFromFile();

      const TH1* getFirstPtEtaHisto() const; /// get the first available pt_eta histo
      const TAxis* getPtAxis() const;  /// only consider pt_eta histos; assume all histos have the same binning
      const TAxis* getEtaAxis() const; /// only consider pt_eta histos; assume all histos have the same binning
      bool getHistoAndParametrization(const Lepton &lep,
                                      const size_t regionIndex,
                                      const RATE_TYPE &rt,
                                      TH1* &h,
                                      Parametrization::Value &rp) const;
      float getFracRelativeError(const Lepton &lep, RATE_TYPE rt, size_t regionIndex,
                                 Systematic::Value syst) const;

      /// input file holding the real efficiency and fake rates for leptons
      TFile* m_hist_file;
      /// names of the signal regions for which we can compute the fake weight
      std::vector<std::string> m_signalRegions;
      std::vector<TH1*> m_el_real_eff;
      std::vector<TH1*> m_el_fake_rate;
      std::vector<TH1*> m_mu_real_eff;
      std::vector<TH1*> m_mu_fake_rate;
      TH1* m_el_frac_up;
      TH1* m_el_frac_do;
      TH1* m_mu_frac_up;
      TH1* m_mu_frac_do;

      // Systematic uncertainties grabbed from the config file
      // Real errors
      double m_el_real_up;
      double m_el_real_down;
      double m_mu_real_up;
      double m_mu_real_down;
      // Fake errors
      double m_el_datamc;
      double m_mu_datamc;
      double m_el_region;
      double m_mu_region;
      TH1* m_el_metrel;   ///< Parameterized vs Met Rel for syst (obsolete?)
      TH1* m_mu_metrel;
      TH1* m_el_eta;     ///< Parameterized vs Eta for syst (obsolete?)
      TH1* m_mu_eta;
      double m_el_HFLFerr;

      // rate parameterization that describes how the rates are binned
      //Parametrization::Value m_rate_param;
      Parametrization::Value m_rate_param_real_el;
      Parametrization::Value m_rate_param_fake_el;
      Parametrization::Value m_rate_param_real_mu;
      Parametrization::Value m_rate_param_fake_mu;
  };
} // fake
} // susy

#endif
