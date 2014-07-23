#include "DileptonMatrixMethod/MatrixLepton.h"

// -----------------------------------------------------------------------------
susy::fake::MatrixLepton::MatrixLepton( bool isTight
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
susy::fake::MatrixLepton::~MatrixLepton()
{
  // do nothing
}

// -----------------------------------------------------------------------------
bool susy::fake::MatrixLepton::isTight() const
{
  return m_isTight;
}

// -----------------------------------------------------------------------------
bool susy::fake::MatrixLepton::isElectron() const
{
  return m_isElectron;
}

// -----------------------------------------------------------------------------
bool susy::fake::MatrixLepton::isMuon() const
{
  return !m_isElectron;
}

// -----------------------------------------------------------------------------
float susy::fake::MatrixLepton::pt() const
{
  return m_pt;
}

// -----------------------------------------------------------------------------
float susy::fake::MatrixLepton::eta() const
{
  return m_eta;
}

// -----------------------------------------------------------------------------
bool susy::fake::MatrixLepton::operator<(
    const susy::fake::MatrixLepton& rhs) const
{
  return (pt() < rhs.pt());
}

// -----------------------------------------------------------------------------
bool susy::fake::MatrixLepton::operator>(
    const susy::fake::MatrixLepton& rhs) const
{
  return (pt() > rhs.pt());
}
