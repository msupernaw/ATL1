/* 
 * File:   Primitive.hpp
 * Author: matthewsupernaw
 *
 * Created on December 8, 2014, 3:38 PM
 */

#ifndef PRIMITIVE_HPP
#define	PRIMITIVE_HPP

namespace atl {

    template<class T>
    bool IsPrimitive() {
        return false;
    }

#define DEFINE_IS_PRIMITIVE(TYPE) \
template<>                    \
bool IsPrimitive<TYPE>(){               \
    return true;                  \
}                              \

    DEFINE_IS_PRIMITIVE(char);
    DEFINE_IS_PRIMITIVE(unsigned char);
    DEFINE_IS_PRIMITIVE(short);
    DEFINE_IS_PRIMITIVE(unsigned short);
    DEFINE_IS_PRIMITIVE(int);
    DEFINE_IS_PRIMITIVE(unsigned int);
    DEFINE_IS_PRIMITIVE(long);
    DEFINE_IS_PRIMITIVE(unsigned long);
    DEFINE_IS_PRIMITIVE(float);
    DEFINE_IS_PRIMITIVE(double);
    DEFINE_IS_PRIMITIVE(long double);

}
#endif	/* PRIMITIVE_HPP */

