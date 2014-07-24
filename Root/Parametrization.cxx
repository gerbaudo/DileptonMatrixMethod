#include "DileptonMatrixMethod/Parametrization.h"

using susy::fake::Parametrization;
//----------------------------------------------------------
bool Parametrization::isValid(const Parametrization::Value &p)
{
    const Parametrization::Value first = Parametrization::PT_ETA;
    const Parametrization::Value last = Parametrization::PT;
    return p>=first && p<=last; // note this won't work if you assign values by hand
}
//----------------------------------------------------------
std::string Parametrization::str(const Parametrization::Value &p)
{
    std::string name="unknown";
    switch(p) {
    case Parametrization::PT_ETA : name = "PT_ETA"; break;
    case Parametrization::PT     : name = "PT";     break;
    // do not put default, so that the compiler will warn on forgotten cases
    }
    return name;
}
//----------------------------------------------------------
