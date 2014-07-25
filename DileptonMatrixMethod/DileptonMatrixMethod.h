// -*- c++ -*-
#ifndef SUSY_FAKE_DILEPTONMATRIXMETHOD_H
#define SUSY_FAKE_DILEPTONMATRIXMETHOD_H

/**
 * @author Brett Jackson <Brett.David.Jackson@cern.ch>
 * @data July 2012
 *
 * Class to perform the matrix method fake estimate
 */


#include "DileptonMatrixMethod/MatrixLepton.h"
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
class DileptonMatrixMethod
  {
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
      /// Get the total fake contribution for this event -- called by the user
      float getTotalFake(bool isTight1, bool isElectron1, float pt1, float eta1,
                         bool isTight2, bool isElectron2, float pt2, float eta2,
                         const size_t regionIndex,
                         float MetRel,
                         Systematic::Value syst = Systematic::SYS_NOM) const;

      /// Get the RR contribution for this event -- called by the user
      float getRR(bool isTight1, bool isElectron1, float pt1, float eta1,
                  bool isTight2, bool isElectron2, float pt2, float eta2,
                  const size_t regionIndex,
                  float MetRel,
                  Systematic::Value syst = Systematic::SYS_NOM) const;

      /// Get the RF contribution for this event -- called by the user
      float getRF(bool isTight1, bool isElectron1, float pt1, float eta1,
                  bool isTight2, bool isElectron2, float pt2, float eta2,
                  const size_t regionIndex,
                  float MetRel,
                  Systematic::Value syst = Systematic::SYS_NOM) const;

      /**
       * Get the FR contribution for this event -- called by the user
       */
      float getFR(bool isTight1, bool isElectron1, float pt1, float eta1,
                  bool isTight2, bool isElectron2, float pt2, float eta2,
                  const size_t regionIndex,
                  float MetRel,
                  Systematic::Value syst = Systematic::SYS_NOM) const;

      /// Get the FF contribution for this event -- called by the user
      float getFF(bool isTight1, bool isElectron1, float pt1, float eta1,
                  bool isTight2, bool isElectron2, float pt2, float eta2,
                  const size_t regionIndex,
                  float MetRel,
                  Systematic::Value syst = Systematic::SYS_NOM) const;
      /// given a region, determine the internal index used to store its histograms; abort if invalid
      size_t getIndexRegion(const std::string &regionName) const;

      enum RATE_TYPE { REAL
                     , FAKE
                     };

      /// Get the rate for this lepton -- real/fake for electron or muon
      /**
         Specific rate depends on type of lepton supplied and the RATE_TYPE parameter
      */
      float getRate(bool isTight,
                    bool isElectron,
                    float pt,
                    float eta,
                    RATE_TYPE rate_type,                     
                    size_t regionIndex,
                    float MetRel,
                    Systematic::Value syst = Systematic::SYS_NOM) const;

      const TArrayD* getPtBins() const;
      const TArrayD* getEtaBins() const;
      void printRateSystematics(const MatrixLepton &l, RATE_TYPE &rt, size_t regionIndex) const;
  protected:
      // Get the rate for this lepton -- real/fake for electron or muon
      // Specific rate depends on type of lepton supplied and the
      // RATE_TYPE parameter
      // for internal use
      float getRate(const MatrixLepton&,
                    RATE_TYPE,
                    size_t regionIndex,
                    float MetRel,
                    Systematic::Value syst = Systematic::SYS_NOM) const;
      float getRateSyst(const MatrixLepton&,
                        RATE_TYPE,
                        size_t regionIndex,
                        float MetRel,
                        Systematic::Value syst = Systematic::SYS_NOM) const;

      int getRateBin(const MatrixLepton& lep,
                     TH1* h_rate,
                     Parametrization::Value rate_param) const;

      // Get the fake/real contribution for this event -- for internal use
      float getTotalFake(const MatrixLepton& lep1,
                         const MatrixLepton& lep2,
                         size_t regionIndex,
                         float MetRel,
                         Systematic::Value syst = Systematic::SYS_NOM) const;
      float getRR(const MatrixLepton& lep1,
                  const MatrixLepton& lep2,
                  size_t regionIndex,
                  float MetRel,
                  Systematic::Value syst = Systematic::SYS_NOM) const;
      float getRF(const MatrixLepton& lep1,
                  const MatrixLepton& lep2,
                  size_t regionIndex,
                  float MetRel,
                  Systematic::Value syst = Systematic::SYS_NOM) const;
      float getFR(const MatrixLepton& lep1,
                  const MatrixLepton& lep2,
                  size_t regionIndex,
                  float MetRel,
                  Systematic::Value syst = Systematic::SYS_NOM) const;
      float getFF(const MatrixLepton& lep1,
                  const MatrixLepton& lep2,
                  size_t regionIndex,
                  float MetRel,
                  Systematic::Value syst = Systematic::SYS_NOM) const;

      // Additional methods needed for systematics
      float getStatError(const MatrixLepton&, RATE_TYPE, size_t regionIndex) const;
      float getRelStatError(const MatrixLepton&, RATE_TYPE, size_t regionIndex) const;

      // Determine if the event is TT/TL/LT/LL -- for internal use
      int getTT(const MatrixLepton& lep1, const MatrixLepton& lep2) const;
      int getTL(const MatrixLepton& lep1, const MatrixLepton& lep2) const;
      int getLT(const MatrixLepton& lep1, const MatrixLepton& lep2) const;
      int getLL(const MatrixLepton& lep1, const MatrixLepton& lep2) const;

      // Verbose output
      void printInfo(const MatrixLepton& lep1,
                     const MatrixLepton& lep2,
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
      bool getHistoAndParametrization(const MatrixLepton &lep,
                                      const size_t regionIndex,
                                      const RATE_TYPE &rt,
                                      TH1* &h,
                                      Parametrization::Value &rp) const;
      float getFracRelativeError(const MatrixLepton &lep, RATE_TYPE rt, size_t regionIndex,
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
