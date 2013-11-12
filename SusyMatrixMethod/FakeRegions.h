// emacs -*- C++ -*-
#ifndef SUSY_FAKE_FAKEREGIONS_H
#define SUSY_FAKE_FAKEREGIONS_H

#include <string>


// binnings used to derive the fake estimate

namespace susy {
namespace fake {

enum Region
{
  // DG These are needed to compute the SF and the rates (pseudo t&p)
  CR_Real = 0,       // Real Z window
  CR_SideLow,        // Real Side band Lower
  CR_SideHigh,       // Real Side band higher
  CR_HF,             // HF b-bbar tag and probe
  CR_HF_high,        // HF b-bbar tag and probe with met < 80
  CR_Conv,           // Zmumu tag and probe
  CR_MCConv,         // MC cr conversions
  CR_MCQCD,          // MC cr QCD = LF+HF
  CR_MCReal,         // MC cr Real
  // DG These are the SR, which we need b/c we want to compute the compositions.
  CR_SSInc,
  CR_SRWHSS,
  CR_CR8lpt,
  CR_CR8ee,
  CR_CR8mm,
  CR_CR8mmMtww,
  CR_CR8mmHt,
  CR_SsEwk,
  CR_SsEwkLoose,
};

const Region SignalRegions[] = {
    CR_SSInc, CR_SRWHSS,
    CR_CR8lpt, CR_CR8ee, CR_CR8mm, CR_CR8mmMtww, CR_CR8mmHt,
    CR_SsEwk, CR_SsEwkLoose
};
const int NumberOfSignalRegions = sizeof(SignalRegions) / sizeof(SignalRegions[0]);
const string RegionNames[] =
{
  "realCR",
  "realSideLow",
  "realSideHigh",
  "fakeHF",
  "fakeHF_high",
  "fakeConv",
  "convMC",
  "qcdMC",
  "realMC",
  "CR_SSInc",
  "CR_WHSS",
  "CR_CR8lpt",
  "CR_CR8ee",
  "CR_CR8mm",
  "CR_CR8mmMtww",
  "CR_CR8mmHt",
  "CR_SsEwk",
  "CR_SsEwkLoose",
};

inline std::string region2str(const Region &r) {return RegionNames[r];}

} // end namespace fake
} // end namespace susy

#endif
