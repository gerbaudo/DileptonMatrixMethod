// -*- c++ -*-
#ifndef SUSY_FAKE_PARAMETRIZATION_H
#define SUSY_FAKE_PARAMETRIZATION_H

#include <string>

namespace susy{
namespace fake{

/// Possible fake matrix parametrizations
/**
   Two values: vs. eta (1D) or vs. pt vs eta (2D).
   Based on the original code by Matt and Brett.

   Note to self: enclose enum in a struct to avoid collisions; see
   explanation on
   <A HREF="http://stackoverflow.com/questions/7090130/enum-in-a-namespace">stackoveflow</A>.
   @author  davide.gerbaudo@gmail.com
   @data July 2014
*/
struct Parametrization {
    enum Value {
        PT_ETA ///< 2D parametrization vs. eta vs. pt
        ,PT    ///< 1D parametrization vs. eta
    };
    /// whether p is a valid enum value
    static bool isValid(const Parametrization::Value &p);
    /// string representation
    static std::string str(const Parametrization::Value &p);
};

} // fake
} // susy
#endif
