// -*- c++ -*-
#ifndef SUSY_FAKE_SYSTEMATIC_H
#define SUSY_FAKE_SYSTEMATIC_H

#include <string>

namespace susy{
namespace fake{

/// Possible fake matrix systematic variations
/**
   Note that some systematic variations might need additional input
   histograms. If they are not available, the DileptonMatrixMethod
   class will warn you when configuring.

   The ones indicated with 'histo' below require additional input
   histograms.

   @author  davide.gerbaudo@gmail.com
   @data July 2014
*/
struct Systematic {
    enum Value {
      SYS_NOM    = 0       ///< nominal
      ,SYS_EL_RE_UP        ///< electron real up
      ,SYS_EL_RE_DOWN      ///< electron real down
      ,SYS_EL_FR_UP        ///< electron fake up
      ,SYS_EL_FR_DOWN      ///< electron fake down
      ,SYS_MU_RE_UP        ///< muon real up
      ,SYS_MU_RE_DOWN      ///< muon real down
      ,SYS_MU_FR_UP        ///< muon fake up
      ,SYS_MU_FR_DOWN      ///< muon fake down
      ,SYS_EL_DATAMC_UP    ///< electron data/mc rate differences up
      ,SYS_EL_DATAMC_DOWN  ///< electron data/mc rate differences down
      ,SYS_MU_DATAMC_UP    ///< muon data/mc rate differences up
      ,SYS_MU_DATAMC_DOWN  ///< muon data/mc rate differences down
      //, SYS_EL_METREL       // Any metrel dep
      //, SYS_MU_METREL       // Any metrel dep
      ,SYS_EL_FR_STAT_UP   ///< electron stat error up
      ,SYS_EL_FR_STAT_DOWN ///< electron stat error down
      ,SYS_MU_FR_STAT_UP   ///< muon stat error up
      ,SYS_MU_FR_STAT_DOWN ///< muon stat error down
      ,SYS_EL_ETA          ///< electron eta dep (obsolete when using eta-dependent scale factors)
      ,SYS_MU_ETA          ///< muon eta dep
      ,SYS_EL_HFLF_UP      ///< electron diff between HF/LF up (obsolete when using HF/LF splitting)
      ,SYS_EL_HFLF_DOWN    ///< electron diff between HF/LF down
      ,SYS_EL_REG_UP       ///< Error from percentage in weighted avg (from Matt; obsolete?)
      ,SYS_EL_REG_DOWN     ///< Error from percentage in weighted avg
      ,SYS_MU_REG_UP       ///< Error from percentage in weighted avg
      ,SYS_MU_REG_DOWN     ///< Error from percentage in weighted avg
      ,SYS_EL_FRAC_UP      ///< electron uncertainty on the fraction determination (up)
      ,SYS_EL_FRAC_DO      ///< electron uncertainty on the fraction determination (down)
      ,SYS_MU_FRAC_UP      ///< muon uncertainty on the fraction determination (up)
      ,SYS_MU_FRAC_DO      ///< muon uncertainty on the fraction determination (down)
    };
  /// first valid value (useful for loop and validation)
  static const Value first() { return SYS_NOM; }
  /// last valid value (needs to be updated when enum is extended)
  static const Value last() { return SYS_MU_FRAC_DO; }
  /// whether p is a valid enum value
  static bool isValid(const Systematic::Value &p);
  /// string representation
  static std::string str(const Systematic::Value &p);
  /// whether this variation requires additional input histograms
  static bool requiresHistogram(const Systematic::Value &p);

};

} // fake
} // susy
#endif
