#include "SameSignMatrixMethod/MatrixLepton.h"

// -----------------------------------------------------------------------------
SameSignMatrixMethod::MatrixLepton::MatrixLepton( bool isTight
                                            , bool isElectron
                                            , float pt
                                            , float eta
                                            )
                                            : m_isTight(isTight)
                                            , m_isElectron(isElectron)
                                            , m_pt(pt)
                                            , m_eta(eta)
{
  // do nothing
}

// -----------------------------------------------------------------------------
SameSignMatrixMethod::MatrixLepton::~MatrixLepton()
{
  // do nothing
}

// -----------------------------------------------------------------------------
bool SameSignMatrixMethod::MatrixLepton::isTight() const
{
  return m_isTight;
}

// -----------------------------------------------------------------------------
bool SameSignMatrixMethod::MatrixLepton::isElectron() const
{
  return m_isElectron;
}

// -----------------------------------------------------------------------------
bool SameSignMatrixMethod::MatrixLepton::isMuon() const
{
  return !m_isElectron;
}

// -----------------------------------------------------------------------------
float SameSignMatrixMethod::MatrixLepton::pt() const
{
  return m_pt;
}

// -----------------------------------------------------------------------------
float SameSignMatrixMethod::MatrixLepton::eta() const
{
  return m_eta;
}

// -----------------------------------------------------------------------------
bool SameSignMatrixMethod::MatrixLepton::operator<(
    const SameSignMatrixMethod::MatrixLepton& rhs) const
{
  return (pt() < rhs.pt());
}

// -----------------------------------------------------------------------------
bool SameSignMatrixMethod::MatrixLepton::operator>(
    const SameSignMatrixMethod::MatrixLepton& rhs) const
{
  return (pt() > rhs.pt());
}
