#ifndef SusyMatrixMethod_DiLeptonMatrixMethod_h
#define SusyMatrixMethod_DiLeptonMatrixMethod_h

/**
 * @author Brett Jackson <Brett.David.Jackson@cern.ch>
 * @data July 2012
 *
 * Class to perform the matrix method fake estimate
 */

#include <iostream>
#include <map>
#include <math.h>
#include <string>

#include "SusyMatrixMethod/MatrixLepton.h"
#include "SusyMatrixMethod/FakeRegions.h"

class TFile;
class TH1;
class TH1F;

namespace SusyMatrixMethod
{
  enum RATE_PARAM { PT_ETA
                  , PT
                  };

  enum SYSTEMATIC { SYS_NOM    = 0
                  , SYS_EL_RE_UP
                  , SYS_EL_RE_DOWN
                  , SYS_EL_FR_UP
                  , SYS_EL_FR_DOWN
                  , SYS_MU_RE_UP
                  , SYS_MU_RE_DOWN
                  , SYS_MU_FR_UP
                  , SYS_MU_FR_DOWN
                  , SYS_N_USER
		  , SYS_EL_DATAMC_UP    // Take into account data/mc differences in rates
		  , SYS_EL_DATAMC_DOWN  // Take into account data/mc differences in rates
		  , SYS_MU_DATAMC_UP    // Take into account data/mc differences in rates
		  , SYS_MU_DATAMC_DOWN  // Take into account data/mc differences in rates
		    //, SYS_EL_METREL       // Any metrel dep
		    //, SYS_MU_METREL       // Any metrel dep
		  , SYS_EL_FR_STAT_UP   // Stat error up
		  , SYS_EL_FR_STAT_DOWN // Stat error down
		  , SYS_MU_FR_STAT_UP   // Stat error up
		  , SYS_MU_FR_STAT_DOWN // Stat error down
		  , SYS_EL_ETA          // Any eta dep
		  , SYS_MU_ETA          // Any eta dep
                  , SYS_EL_HFLF_UP      // Needed for HF/LF diff for elec
                  , SYS_EL_HFLF_DOWN    // Needed for HF/LF diff for elec
		  , SYS_EL_REG_UP       // Error from percentage in weighted avg
		  , SYS_EL_REG_DOWN     // Error from percentage in weighted avg
		  , SYS_MU_REG_UP       // Error from percentage in weighted avg
		  , SYS_MU_REG_DOWN     // Error from percentage in weighted avg
                  , SYS_N
                  };

  static std::string systematic_names[] = { "NOM"
                                         , "EL_RE_UP"
                                         , "EL_RE_DOWN"
                                         , "EL_FR_UP"
                                         , "EL_FR_DOWN"
                                         , "MU_RE_UP"
                                         , "MU_RE_DOWN"
                                         , "MU_FR_UP"
                                         , "MU_FR_DOWN"
					 , "SYS_N_USER"
				         , "EL_DATAMC_UP"
				         , "EL_DATAMC_DOWN"
				         , "MU_DATAMC_UP"
				         , "MU_DATAMC_DOWN"
					    //, "EL_METREL"
					    //, "MU_METREL"
                                         , "EL_FR_STAT_UP"
                                         , "EL_FR_STAT_DOWN"
                                         , "MU_FR_STAT_UP"
                                         , "MU_FR_STAT_DOWN"
                                         , "EL_ETA"
                                         , "MU_ETA"
                                         , "EL_HFLF_UP"
					 , "EL_HFLF_DOWN"
					 , "EL_REG_UP"
				         , "EL_REG_DOWN"
				         , "MU_REG_UP"
				         , "MU_REG_DOWN"
                                         , "SYS_N"
                                         };

  class DiLeptonMatrixMethod
  {
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    public:
      /**
       * constructor/destructor
       */
      DiLeptonMatrixMethod();
      ~DiLeptonMatrixMethod();

      /**
       * congigure the tool with a histogram file
       */
      bool configure( std::string file_name
                    , RATE_PARAM rate_param_real_el
                    , RATE_PARAM rate_param_fake_el
                    , RATE_PARAM rate_param_real_mu
                    , RATE_PARAM rate_param_fake_mu
                    );

      /**
       * set the units for the rates file
       */
      // void setMeV();
      // void setGeV();

      /**
       * Get the total fake contribution for this event -- called by the user
       */
      float getTotalFake( bool isTight1, bool isElectron1, float pt1, float eta1
                        , bool isTight2, bool isElectron2, float pt2, float eta2
                        , susy::fake::Region region
                        , float MetRel
                        , SYSTEMATIC syst = SYS_NOM
                        ) const;

      /**
       * Get the RR contribution for this event -- called by the user
       */
      float getRR( bool isTight1, bool isElectron1, float pt1, float eta1
                 , bool isTight2, bool isElectron2, float pt2, float eta2
                 , susy::fake::Region region
                 , float MetRel
                 , SYSTEMATIC syst = SYS_NOM
                 ) const;

      /**
       * Get the RF contribution for this event -- called by the user
       */
      float getRF( bool isTight1, bool isElectron1, float pt1, float eta1
                 , bool isTight2, bool isElectron2, float pt2, float eta2
                 , susy::fake::Region region
                 , float MetRel
                 , SYSTEMATIC syst = SYS_NOM
                 ) const;

      /**
       * Get the FR contribution for this event -- called by the user
       */
      float getFR( bool isTight1, bool isElectron1, float pt1, float eta1
                 , bool isTight2, bool isElectron2, float pt2, float eta2
                 , susy::fake::Region region
                 , float MetRel
                 , SYSTEMATIC syst = SYS_NOM
                 ) const;

      /**
       * Get the FF contribution for this event -- called by the user
       */
      float getFF( bool isTight1, bool isElectron1, float pt1, float eta1
                 , bool isTight2, bool isElectron2, float pt2, float eta2
                 , susy::fake::Region region
                 , float MetRel
                 , SYSTEMATIC syst = SYS_NOM
                 ) const;
      //! given a region, determine the internal index used to store its histograms; abort if invalid
      static int getIndexRegion(susy::fake::Region region);
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      enum RATE_TYPE { REAL
                     , FAKE
                     };

      // Get the rate for this lepton -- real/fake for electron or muon
      // Specific rate depends on type of lepton supplied and the
      // RATE_TYPE parameter
      float getRate( bool isTight
                   , bool isElectron
                   , float pt
                   , float eta
                   , RATE_TYPE rate_type
                   , susy::fake::Region region
                   , float MetRel
                   , SYSTEMATIC syst = SYS_NOM
                   ) const;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  protected:
      // Get the rate for this lepton -- real/fake for electron or muon
      // Specific rate depends on type of lepton supplied and the
      // RATE_TYPE parameter
      // for internal use
      float getRate( const MatrixLepton&
                   , RATE_TYPE
                   , susy::fake::Region region
                   , float MetRel
                   , SYSTEMATIC syst = SYS_NOM
                   ) const;
      float getRateSyst( const MatrixLepton&
                       , RATE_TYPE
           , susy::fake::Region region
                       , float MetRel
                       , SYSTEMATIC syst = SYS_NOM
                       ) const;

      int getRateBin( const MatrixLepton& lep,
			TH1* h_rate,
		      RATE_PARAM rate_param
			) const;

      // Get the fake/real contribution for this event -- for internal use
      float getTotalFake( const MatrixLepton& lep1
                        , const MatrixLepton& lep2
                        , susy::fake::Region region
                        , float MetRel
                        , SYSTEMATIC syst = SYS_NOM
                        ) const;
      float getRR( const MatrixLepton& lep1
                 , const MatrixLepton& lep2
                 , susy::fake::Region region
                 , float MetRel
                 , SYSTEMATIC syst = SYS_NOM
                 ) const;
      float getRF( const MatrixLepton& lep1
                 , const MatrixLepton& lep2
                 , susy::fake::Region region
                 , float MetRel
                 , SYSTEMATIC syst = SYS_NOM
                 ) const;
      float getFR( const MatrixLepton& lep1
                 , const MatrixLepton& lep2
                 , susy::fake::Region region
                 , float MetRel
                 , SYSTEMATIC syst = SYS_NOM
                 ) const;
      float getFF( const MatrixLepton& lep1
                 , const MatrixLepton& lep2
                 , susy::fake::Region region
                 , float MetRel
                 , SYSTEMATIC syst = SYS_NOM
                 ) const;

      // Additional methods needed for systematics
      float getStatError(const MatrixLepton&, RATE_TYPE, susy::fake::Region) const;
      float getRelStatError(const MatrixLepton&, RATE_TYPE, susy::fake::Region) const;

      // Determine if the event is TT/TL/LT/LL -- for internal use
      int getTT(const MatrixLepton& lep1, const MatrixLepton& lep2) const;
      int getTL(const MatrixLepton& lep1, const MatrixLepton& lep2) const;
      int getLT(const MatrixLepton& lep1, const MatrixLepton& lep2) const;
      int getLL(const MatrixLepton& lep1, const MatrixLepton& lep2) const;

      // Verbose output
      void printInfo( const MatrixLepton& lep1
                    , const MatrixLepton& lep2
                    , susy::fake::Region region
                    , float MetRel
                    , SYSTEMATIC syst = SYS_NOM
                    ) const;

      // Get systematic value from file -- for internal use
      void loadSysFromFile();

      bool getHistoAndParametrization(const MatrixLepton &lep,
                                      const susy::fake::Region reg,
                                      const RATE_TYPE &rt,
                                      TH1* &h,
                                      RATE_PARAM &rp) const;
      // Histograms which hold the real efficiency and fake rates for leptons
      TFile* m_hist_file;
      TH1* m_el_real_eff [susy::fake::NumberOfSignalRegions];
      TH1* m_el_fake_rate[susy::fake::NumberOfSignalRegions];
      TH1* m_mu_real_eff [susy::fake::NumberOfSignalRegions];
      TH1* m_mu_fake_rate[susy::fake::NumberOfSignalRegions];

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
      TH1F* m_el_metrel;
      TH1F* m_mu_metrel;
      TH1F* m_el_eta;
      TH1F* m_mu_eta;
      double m_el_HFLFerr;

      // rate parameterization that describes how the rates are binned
      //RATE_PARAM m_rate_param;
      RATE_PARAM m_rate_param_real_el;
      RATE_PARAM m_rate_param_fake_el;
      RATE_PARAM m_rate_param_real_mu;
      RATE_PARAM m_rate_param_fake_mu;
  };
}

#endif /* end of include guard: SusyMatrixMethod_DiLeptonMatrixMethod_h */
