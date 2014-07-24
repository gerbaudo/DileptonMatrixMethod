#include "DileptonMatrixMethod/Systematic.h"

using susy::fake::Systematic;
//----------------------------------------------------------
bool Systematic::isValid(const Systematic::Value &p)
{
    return p>=Systematic::first() && p<=Systematic::last();
}
//----------------------------------------------------------
std::string Systematic::str(const Systematic::Value &p)
{
    std::string name="unknown";
    switch(p) {
    case Systematic::SYS_NOM             : name = "NOM"            ; break;
    case Systematic::SYS_EL_RE_UP        : name = "EL_RE_UP"       ; break;
    case Systematic::SYS_EL_RE_DOWN      : name = "EL_RE_DOWN"     ; break;
    case Systematic::SYS_EL_FR_UP        : name = "EL_FR_UP"       ; break;
    case Systematic::SYS_EL_FR_DOWN      : name = "EL_FR_DOWN"     ; break;
    case Systematic::SYS_MU_RE_UP        : name = "MU_RE_UP"       ; break;
    case Systematic::SYS_MU_RE_DOWN      : name = "MU_RE_DOWN"     ; break;
    case Systematic::SYS_MU_FR_UP        : name = "MU_FR_UP"       ; break;
    case Systematic::SYS_MU_FR_DOWN      : name = "MU_FR_DOWN"     ; break;
    case Systematic::SYS_EL_DATAMC_UP    : name = "EL_DATAMC_UP"   ; break;
    case Systematic::SYS_EL_DATAMC_DOWN  : name = "EL_DATAMC_DOWN" ; break;
    case Systematic::SYS_MU_DATAMC_UP    : name = "MU_DATAMC_UP"   ; break;
    case Systematic::SYS_MU_DATAMC_DOWN  : name = "MU_DATAMC_DOWN" ; break;
    case Systematic::SYS_EL_FR_STAT_UP   : name = "EL_FR_STAT_UP"  ; break;
    case Systematic::SYS_EL_FR_STAT_DOWN : name = "EL_FR_STAT_DOWN"; break;
    case Systematic::SYS_MU_FR_STAT_UP   : name = "MU_FR_STAT_UP"  ; break;
    case Systematic::SYS_MU_FR_STAT_DOWN : name = "MU_FR_STAT_DOWN"; break;
    case Systematic::SYS_EL_ETA          : name = "EL_ETA"         ; break;
    case Systematic::SYS_MU_ETA          : name = "MU_ETA"         ; break;
    case Systematic::SYS_EL_HFLF_UP      : name = "EL_HFLF_UP"     ; break;
    case Systematic::SYS_EL_HFLF_DOWN    : name = "EL_HFLF_DOWN"   ; break;
    case Systematic::SYS_EL_REG_UP       : name = "EL_REG_UP"      ; break;
    case Systematic::SYS_EL_REG_DOWN     : name = "EL_REG_DOWN"    ; break;
    case Systematic::SYS_MU_REG_UP       : name = "MU_REG_UP"      ; break;
    case Systematic::SYS_MU_REG_DOWN     : name = "MU_REG_DOWN"    ; break;
    case Systematic::SYS_EL_FRAC_UP      : name = "EL_FRAC_UP"     ; break;
    case Systematic::SYS_EL_FRAC_DO      : name = "EL_FRAC_DO"     ; break;
    case Systematic::SYS_MU_FRAC_UP      : name = "MU_FRAC_UP"     ; break;
    case Systematic::SYS_MU_FRAC_DO      : name = "MU_FRAC_DO"     ; break;
    // do not put default, so that the compiler will warn on forgotten cases
    }
    return name;
}
//----------------------------------------------------------
bool Systematic::requiresHistogram(const Systematic::Value &p)
{
    return (p==SYS_EL_FRAC_UP ||
            p==SYS_EL_FRAC_DO ||
            p==SYS_MU_FRAC_UP ||
            p==SYS_MU_FRAC_DO ||
            p==SYS_EL_ETA     ||
            p==SYS_MU_ETA);
}
//----------------------------------------------------------
