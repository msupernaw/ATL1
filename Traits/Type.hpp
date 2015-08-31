/* 
 * File:   Type.hpp
 * Author: matthewsupernaw
 *
 * Created on December 4, 2014, 3:54 PM
 */

#ifndef TYPE_HPP
#define	TYPE_HPP

namespace atl {

    template<typename T>
    struct VariableTrait {
        static bool is_variable;
        static size_t size;
        typedef  T BASE_TYPE;

        static bool is_recording() {
            return false;
        }
    };

    template<typename T>
    bool VariableTrait<T>::is_variable = false;

    template<typename T>
    size_t VariableTrait<T>::size = sizeof (T);

#define DECLARE_IS_VARIABLE_TRAIT(A, B, C,D,E) template<> struct VariableTrait<A> { \
    static bool is_variable;\
     static size_t size; \
     typedef  E BASE_TYPE; \
    static bool is_recording(){\
            return D;\
        }\
}; \
     bool VariableTrait<A>::is_variable = B; \
    size_t VariableTrait<A>::size = C; \
    //#define DECLARE_IS_VARIABLE_TRAIT(A, B) template<> struct IsVariable<A> { bool is_variable = B;};

    DECLARE_IS_VARIABLE_TRAIT(atl::Variable<long double>, true, sizeof (long double), atl::Variable<long double>::IsRecording(),long double);
    DECLARE_IS_VARIABLE_TRAIT(atl::Variable<double>, true, sizeof (double), atl::Variable<double>::IsRecording(),double);
    DECLARE_IS_VARIABLE_TRAIT(atl::Variable<float>, true, sizeof (float), atl::Variable<float>::IsRecording(),float);
    DECLARE_IS_VARIABLE_TRAIT(long double, false, sizeof (long double), false,long double);
    DECLARE_IS_VARIABLE_TRAIT(double, false, sizeof (double), false,double);
    DECLARE_IS_VARIABLE_TRAIT(float, false, sizeof (float), false,float);
}



#endif	/* TYPE_HPP */

