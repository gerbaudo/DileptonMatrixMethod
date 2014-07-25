#include "DileptonMatrixMethod/Lepton.h"

// -----------------------------------------------------------------------------
susy::fake::Lepton::Lepton( bool isTight
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
susy::fake::Lepton::~Lepton()
{
  // do nothing
}

// -----------------------------------------------------------------------------
bool susy::fake::Lepton::isTight() const
{
  return m_isTight;
}

// -----------------------------------------------------------------------------
bool susy::fake::Lepton::isElectron() const
{
  return m_isElectron;
}

// -----------------------------------------------------------------------------
bool susy::fake::Lepton::isMuon() const
{
  return !m_isElectron;
}

// -----------------------------------------------------------------------------
float susy::fake::Lepton::pt() const
{
  return m_pt;
}

// -----------------------------------------------------------------------------
float susy::fake::Lepton::eta() const
{
  return m_eta;
}

// -----------------------------------------------------------------------------
bool susy::fake::Lepton::operator<(
    const susy::fake::Lepton& rhs) const
{
  return (pt() < rhs.pt());
}

// -----------------------------------------------------------------------------
bool susy::fake::Lepton::operator>(
    const susy::fake::Lepton& rhs) const
{
  return (pt() > rhs.pt());
}
