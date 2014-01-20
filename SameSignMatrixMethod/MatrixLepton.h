#ifndef SameSignMatrixMethod_MatrixLepton_h
#define SameSignMatrixMethod_MatrixLepton_h

/**
 * @author Brett Jackson <Brett.David.Jackson@cern.ch>
 * @data July 2012
 *
 * Helper class to hold lepton properties for doing the matrix method
 * This helper class is mostly borrowed from Steve Farrel's class of
 * the same name
 * Steve Farrel <steven.farrell@cern.ch>
 * https://svnweb.cern.ch/trac/atlasinst/browser/Institutes/UCIrvine/sfarrell/MultiLepMatrixMethod/tags/MultiLepMatrixMethod-00-00-02/MultiLepMatrixMethod/MatrixLepton.h
 * https://svnweb.cern.ch/trac/atlasinst/browser/Institutes/UCIrvine/sfarrell/MultiLepMatrixMethod/tags/MultiLepMatrixMethod-00-00-02/Root/MatrixLepton.cxx
 */

namespace SameSignMatrixMethod
{
  class MatrixLepton
  {
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    public:
      MatrixLepton(bool isTight, bool isElectron, float pt, float eta);
      ~MatrixLepton();

      bool isTight() const;

      // lepton flavor
      bool isElectron() const;
      bool isMuon() const;

      // kinematic variables
      float pt() const;
      float eta() const;

      bool operator< (const SameSignMatrixMethod::MatrixLepton& rhs) const;
      bool operator> (const SameSignMatrixMethod::MatrixLepton& rhs) const;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    private:
      // Is this lepton tight?
      bool m_isTight;

      // Used in parameterizing the fake rates
      bool m_isElectron;
      float m_pt;
      float m_eta;
  };
}

#endif
