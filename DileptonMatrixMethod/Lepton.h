// -*- c++ -*-
#ifndef SUSY_FAKE_LEPTON_H
#define SUSY_FAKE_LEPTON_H

#include <string>

namespace susy{
namespace fake{
/**
  Helper class to hold lepton properties for doing the matrix method
  This helper class is mostly borrowed from Steve Farrel's
  (<steven.farrell@cern.ch>) class MatrixLepton
  https://svnweb.cern.ch/trac/atlasinst/browser/Institutes/UCIrvine/sfarrell/MultiLepMatrixMethod/tags/MultiLepMatrixMethod-00-00-02/MultiLepMatrixMethod/MatrixLepton.h

  Cleanup, rewrite, davide.gerbaudo@gmail.com, July 2014

  @author Brett Jackson <Brett.David.Jackson@cern.ch>
  @data July 2012
*/
class Lepton
{
public:
    /// pt is in GeV
    Lepton(bool isTight, bool isElectron, float pt, float eta);
    bool isTight() const;
    bool isElectron() const;
    bool isMuon() const;
    float pt() const;
    float eta() const;
    bool operator< (const Lepton& rhs) const; ///< pt-based ordering
    bool operator> (const Lepton& rhs) const; ///< pt-based ordering
    Lepton& isTight(bool v); ///< isTight setter
    Lepton& isEl(bool v); ///< isEl setter
    Lepton& isMu(bool v); ///< isMu setter
    Lepton& pt(float v); ///< pt setter (in GeV)
    Lepton& eta(float v); ///< eta setter
    std::string str() const; ///< string representation
private:
    bool m_isTight; ///< whether it satisfies the tight requirements
    bool m_isElectron; ///< either it's an electron or a muon
    float m_pt; ///< pt in GeV
    float m_eta; ///< eta (signed)
};
} // fake
} // susy

#endif
