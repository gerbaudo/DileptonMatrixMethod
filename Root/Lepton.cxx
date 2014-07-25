#include "DileptonMatrixMethod/Lepton.h"

using susy::fake::Lepton;
//----------------------------------------------------------
Lepton::Lepton(bool isTight, bool isElectron, float pt, float eta):
    m_isTight(isTight),
    m_isElectron(isElectron),
    m_pt(pt),
    m_eta(eta)
{
}
//----------------------------------------------------------
bool Lepton::isTight() const
{
    return m_isTight;
}
//----------------------------------------------------------
bool Lepton::isElectron() const
{
    return m_isElectron;
}
//----------------------------------------------------------
bool Lepton::isMuon() const
{
    return !m_isElectron;
}
//----------------------------------------------------------
float Lepton::pt() const
{
    return m_pt;
}
//----------------------------------------------------------
float Lepton::eta() const
{
    return m_eta;
}
//----------------------------------------------------------
bool Lepton::operator<(const Lepton& rhs) const
{
    return (pt() < rhs.pt());
}
//----------------------------------------------------------
bool Lepton::operator>(const Lepton& rhs) const
{
    return (pt() > rhs.pt());
}
//----------------------------------------------------------
Lepton& Lepton::isTight(bool v)
{
    m_isTight=v;
    return *this;
}
//----------------------------------------------------------
Lepton& Lepton::isEl(bool v)
{
    m_isElectron = v;
    return *this;
}
//----------------------------------------------------------
Lepton& Lepton::isMu(bool v)
{
    m_isElectron = !v;
    return *this;
}
//----------------------------------------------------------
Lepton& Lepton::pt(float v)
{
    m_pt = v;
    return *this;
}
//----------------------------------------------------------
Lepton& Lepton::eta(float v)
{
    m_eta = v;
    return *this;
}
//----------------------------------------------------------
