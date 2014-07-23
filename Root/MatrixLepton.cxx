#include "DileptonMatrixMethod/MatrixLepton.h"

// -----------------------------------------------------------------------------
DileptonMatrixMethod::MatrixLepton::MatrixLepton( bool isTight
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
DileptonMatrixMethod::MatrixLepton::~MatrixLepton()
{
  // do nothing
}

// -----------------------------------------------------------------------------
bool DileptonMatrixMethod::MatrixLepton::isTight() const
{
  return m_isTight;
}

// -----------------------------------------------------------------------------
bool DileptonMatrixMethod::MatrixLepton::isElectron() const
{
  return m_isElectron;
}

// -----------------------------------------------------------------------------
bool DileptonMatrixMethod::MatrixLepton::isMuon() const
{
  return !m_isElectron;
}

// -----------------------------------------------------------------------------
float DileptonMatrixMethod::MatrixLepton::pt() const
{
  return m_pt;
}

// -----------------------------------------------------------------------------
float DileptonMatrixMethod::MatrixLepton::eta() const
{
  return m_eta;
}

// -----------------------------------------------------------------------------
bool DileptonMatrixMethod::MatrixLepton::operator<(
    const DileptonMatrixMethod::MatrixLepton& rhs) const
{
  return (pt() < rhs.pt());
}

// -----------------------------------------------------------------------------
bool DileptonMatrixMethod::MatrixLepton::operator>(
    const DileptonMatrixMethod::MatrixLepton& rhs) const
{
  return (pt() > rhs.pt());
}
