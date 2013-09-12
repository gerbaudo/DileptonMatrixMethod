#include "SusyMatrixMethod/MatrixLepton.h"

// -----------------------------------------------------------------------------
SusyMatrixMethod::MatrixLepton::MatrixLepton( bool isTight
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
SusyMatrixMethod::MatrixLepton::~MatrixLepton()
{
  // do nothing
}

// -----------------------------------------------------------------------------
bool SusyMatrixMethod::MatrixLepton::isTight() const
{
  return m_isTight;
}

// -----------------------------------------------------------------------------
bool SusyMatrixMethod::MatrixLepton::isElectron() const
{
  return m_isElectron;
}

// -----------------------------------------------------------------------------
bool SusyMatrixMethod::MatrixLepton::isMuon() const
{
  return !m_isElectron;
}

// -----------------------------------------------------------------------------
float SusyMatrixMethod::MatrixLepton::pt() const
{
  return m_pt;
}

// -----------------------------------------------------------------------------
float SusyMatrixMethod::MatrixLepton::eta() const
{
  return m_eta;
}

// -----------------------------------------------------------------------------
bool SusyMatrixMethod::MatrixLepton::operator<(
    const SusyMatrixMethod::MatrixLepton& rhs) const
{
  return (pt() < rhs.pt());
}

// -----------------------------------------------------------------------------
bool SusyMatrixMethod::MatrixLepton::operator>(
    const SusyMatrixMethod::MatrixLepton& rhs) const
{
  return (pt() > rhs.pt());
}
