/* 
 * File:   Traits.hpp
 * Author: matthewsupernaw
 *
 * Created on November 14, 2014, 12:11 PM
 */

#ifndef TRAITS_HPP
#define	TRAITS_HPP
#include "../AutoDiff/AutoDiff.hpp"

//#define ATL_PROMOTE_BINARY_OPERATIONS

namespace atl {


    template<class REAL_T, int group>
    class Variable;

    template<typename T>
    struct VariableBase {
    };

    template<>
    struct VariableBase<atl::Variable<long double> > {
        typedef long double base_type;
    };

    template<typename T>
    struct VariableBase<atl::ExpressionBase<long double, T> > {
        typedef long double base_type;
    };

    template<>
    struct VariableBase<long double > {
        typedef long double base_type;
    };

    template<typename T>
    struct VariableBase<atl::ExpressionBase<double, T> > {
        typedef long double base_type;
    };

    template<>
    struct VariableBase<atl::Variable<double> > {
        typedef double base_type;
    };

    template<>
    struct VariableBase<double > {
        typedef double base_type;
    };

    template<typename T>
    struct VariableBase<atl::ExpressionBase<float, T> > {
        typedef long double base_type;
    };

    template<>
    struct VariableBase<atl::Variable<float> > {
        typedef float base_type;
    };

    template<>
    struct VariableBase<float > {
        typedef float base_type;
    };

    template<typename T1, typename T2>
    struct PromoteType {
        typedef T1 return_type;
    };



#define DECLARE_PROMOTION(A, B, C) template<> struct PromoteType<A, B> { typedef C return_type;};


    DECLARE_PROMOTION(atl::Variable<long double>, atl::Variable<long double>, atl::Variable<long double>)
    DECLARE_PROMOTION(atl::Variable<long double>, atl::Variable< double>, atl::Variable<long double>);
    DECLARE_PROMOTION(atl::Variable< double>, atl::Variable<long double>, atl::Variable<long double>);
    DECLARE_PROMOTION(atl::Variable<long double>, atl::Variable<float>, atl::Variable<long double>);
    DECLARE_PROMOTION(atl::Variable<float>, atl::Variable<long double>, atl::Variable<long double>);
    DECLARE_PROMOTION(atl::Variable<long double>, long double, atl::Variable<long double>);
    DECLARE_PROMOTION(long double, atl::Variable<long double>, atl::Variable<long double>);
    DECLARE_PROMOTION(atl::Variable<long double>, double, atl::Variable<long double>);
    DECLARE_PROMOTION(double, atl::Variable<long double>, atl::Variable<long double>);
    DECLARE_PROMOTION(atl::Variable<long double>, float, atl::Variable<long double>);
    DECLARE_PROMOTION(float, atl::Variable<long double>, atl::Variable<long double>);
    DECLARE_PROMOTION(atl::Variable<long double>, long, atl::Variable<long double>);
    DECLARE_PROMOTION(long, atl::Variable<long double>, atl::Variable<long double>);
    DECLARE_PROMOTION(atl::Variable<long double>, int, atl::Variable<long double>);
    DECLARE_PROMOTION(int, atl::Variable<long double>, atl::Variable<long double>);
    DECLARE_PROMOTION(atl::Variable<long double>, short, atl::Variable<long double>);
    DECLARE_PROMOTION(short, atl::Variable<long double>, atl::Variable<long double>);
    DECLARE_PROMOTION(atl::Variable< double>, atl::Variable< double>, atl::Variable< double>);
    DECLARE_PROMOTION(atl::Variable< double>, atl::Variable<float>, atl::Variable< double>);
    DECLARE_PROMOTION(atl::Variable<float>, atl::Variable< double>, atl::Variable< double>);
    DECLARE_PROMOTION(atl::Variable< double>, long double, atl::Variable< double>);
    DECLARE_PROMOTION(long double, atl::Variable< double>, atl::Variable< double>);
    DECLARE_PROMOTION(atl::Variable< double>, double, atl::Variable< double>);
    DECLARE_PROMOTION(double, atl::Variable< double>, atl::Variable< double>);
    DECLARE_PROMOTION(atl::Variable< double>, float, atl::Variable< double>);
    DECLARE_PROMOTION(float, atl::Variable< double>, atl::Variable< double>);
    DECLARE_PROMOTION(atl::Variable< double>, long, atl::Variable< double>);
    DECLARE_PROMOTION(long, atl::Variable< double>, atl::Variable< double>);
    DECLARE_PROMOTION(atl::Variable< double>, int, atl::Variable< double>);
    DECLARE_PROMOTION(int, atl::Variable< double>, atl::Variable< double>);
    DECLARE_PROMOTION(atl::Variable< double>, short, atl::Variable< double>);
    DECLARE_PROMOTION(short, atl::Variable< double>, atl::Variable< double>);
    DECLARE_PROMOTION(atl::Variable<float>, atl::Variable<float>, atl::Variable<float>);
    DECLARE_PROMOTION(atl::Variable<float>, long double, atl::Variable<float>);
    DECLARE_PROMOTION(long double, atl::Variable<float>, atl::Variable<float>);
    DECLARE_PROMOTION(atl::Variable<float>, double, atl::Variable<float>);
    DECLARE_PROMOTION(double, atl::Variable<float>, atl::Variable<float>);
    DECLARE_PROMOTION(atl::Variable<float>, float, atl::Variable<float>);
    DECLARE_PROMOTION(float, atl::Variable<float>, atl::Variable<float>);
    DECLARE_PROMOTION(atl::Variable<float>, long, atl::Variable<float>);
    DECLARE_PROMOTION(long, atl::Variable<float>, atl::Variable<float>);
    DECLARE_PROMOTION(atl::Variable<float>, int, atl::Variable<float>);
    DECLARE_PROMOTION(int, atl::Variable<float>, atl::Variable<float>);
    DECLARE_PROMOTION(atl::Variable<float>, short, atl::Variable<float>);
    DECLARE_PROMOTION(short, atl::Variable<float>, atl::Variable<float>);
    DECLARE_PROMOTION(long double, long double, long double);
    DECLARE_PROMOTION(long double, double, long double);
    DECLARE_PROMOTION(double, long double, long double);
    DECLARE_PROMOTION(long double, float, long double);
    DECLARE_PROMOTION(float, long double, long double);
    DECLARE_PROMOTION(long double, long, long double);
    DECLARE_PROMOTION(long, long double, long double);
    DECLARE_PROMOTION(long double, int, long double);
    DECLARE_PROMOTION(int, long double, long double);
    DECLARE_PROMOTION(long double, short, long double);
    DECLARE_PROMOTION(short, long double, long double);
    DECLARE_PROMOTION(double, double, double);
    DECLARE_PROMOTION(double, float, double);
    DECLARE_PROMOTION(float, double, double);
    DECLARE_PROMOTION(double, long, double);
    DECLARE_PROMOTION(long, double, double);
    DECLARE_PROMOTION(double, int, double);
    DECLARE_PROMOTION(int, double, double);
    DECLARE_PROMOTION(double, short, double);
    DECLARE_PROMOTION(short, double, double);
    DECLARE_PROMOTION(float, float, float);
    DECLARE_PROMOTION(float, long, float);
    DECLARE_PROMOTION(long, float, float);
    DECLARE_PROMOTION(float, int, float);
    DECLARE_PROMOTION(int, float, float);
    DECLARE_PROMOTION(float, short, float);
    DECLARE_PROMOTION(short, float, float);
    DECLARE_PROMOTION(long, long, long);
    DECLARE_PROMOTION(long, int, long);
    DECLARE_PROMOTION(int, long, long);
    DECLARE_PROMOTION(long, short, long);
    DECLARE_PROMOTION(short, long, long);
    DECLARE_PROMOTION(int, int, int);
    DECLARE_PROMOTION(int, short, int);
    DECLARE_PROMOTION(short, int, int);
    DECLARE_PROMOTION(short, short, short);

    /**
     * Used for promoting \ref variable
     * expression templates rather than \ref variable
     * types which will make things slow.
     */
    enum PromotableOperator {
        ANY = 0,
        ADD,
        SUBTRACT,
        MULTIPLY,
        DIVIDE
    };

    /**
     * Trait for operator types for containers. The idea
     * is to create traits so that \b ATL \ref variable
     * type operations are returned rather than \ref variable
     * types themselves. Array, Vector, and Matrix operations 
     * will result in many unnecessary temporary \ref variable types 
     * otherwise.
     */
    template<typename T1, typename T2, PromotableOperator PO>
    struct PromoteBinaryOpReturnType {
        typedef typename atl::PromoteType<T1, T2>::return_type return_type;
    };



#ifdef ATL_PROMOTE_BINARY_OPERATIONS

#define DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(A,B,PO)    \
template<> struct PromoteBinaryOpReturnType< A,B,PO >{typedef typename PromoteType< A, B>::return_type return_type;}; \
    
#define DECLARE_VARIABLE_VARIABLE_OPERATOR(PO, R)    \
template<class  T1, class T2> struct PromoteBinaryOpReturnType< T1,T2,PO >{typedef R<typename PromoteType<typename T1::BASE_TYPE,typename T2::BASE_TYPE>::return_type ,T1,T2 > return_type;}; \
  
#define DECLARE_SCALAR_VARIABLE_OPERATOR(A,PO, R)    \
template<typename T1 > struct PromoteBinaryOpReturnType<T1, A, PO>{typedef R<typename PromoteType<typename T1::BASE_TYPE, A>::return_type,T1 > return_type;}; \

#define DECLARE_VARIABLE_SCALAR_OPERATOR(A,PO, R)    \
template<typename T1> struct PromoteBinaryOpReturnType<A,T1, PO>{typedef  R<typename PromoteType<typename T1::BASE_TYPE, A>::return_type,T1 > return_type; }; \

    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(long double, long double, ADD);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(long double, double, ADD);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(double, long double, ADD);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(long double, float, ADD);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(float, long double, ADD);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(long double, long, ADD);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(long, long double, ADD);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(long double, int, ADD);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(int, long double, ADD);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(long double, short, ADD);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(short, long double, ADD);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(double, double, ADD);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(double, float, ADD);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(float, double, ADD);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(double, long, ADD);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(long, double, ADD);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(double, int, ADD);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(int, double, ADD);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(double, short, ADD);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(short, double, ADD);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(float, float, ADD);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(float, long, ADD);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(long, float, ADD);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(float, int, ADD);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(int, float, ADD);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(float, short, ADD);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(short, float, ADD);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(long, long, ADD);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(long, int, ADD);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(int, long, ADD);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(long, short, ADD);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(short, long, ADD);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(int, int, ADD);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(int, short, ADD);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(short, int, ADD);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(short, short, ADD);
    DECLARE_VARIABLE_VARIABLE_OPERATOR(ADD, atl::Add);
    DECLARE_SCALAR_VARIABLE_OPERATOR(long double, ADD, atl::AddConstant);
    DECLARE_SCALAR_VARIABLE_OPERATOR(double, ADD, atl::AddConstant);
    DECLARE_SCALAR_VARIABLE_OPERATOR(float, ADD, atl::AddConstant);
    DECLARE_VARIABLE_SCALAR_OPERATOR(long double, ADD, atl::ConstantAdd);
    DECLARE_VARIABLE_SCALAR_OPERATOR(double, ADD, atl::ConstantAdd);
    DECLARE_VARIABLE_SCALAR_OPERATOR(float, ADD, atl::ConstantAdd);

    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(long double, long double, SUBTRACT);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(long double, double, SUBTRACT);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(double, long double, SUBTRACT);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(long double, float, SUBTRACT);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(float, long double, SUBTRACT);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(long double, long, SUBTRACT);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(long, long double, SUBTRACT);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(long double, int, SUBTRACT);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(int, long double, SUBTRACT);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(long double, short, SUBTRACT);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(short, long double, SUBTRACT);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(double, double, SUBTRACT);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(double, float, SUBTRACT);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(float, double, SUBTRACT);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(double, long, SUBTRACT);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(long, double, SUBTRACT);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(double, int, SUBTRACT);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(int, double, SUBTRACT);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(double, short, SUBTRACT);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(short, double, SUBTRACT);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(float, float, SUBTRACT);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(float, long, SUBTRACT);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(long, float, SUBTRACT);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(float, int, SUBTRACT);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(int, float, SUBTRACT);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(float, short, SUBTRACT);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(short, float, SUBTRACT);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(long, long, SUBTRACT);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(long, int, SUBTRACT);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(int, long, SUBTRACT);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(long, short, SUBTRACT);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(short, long, SUBTRACT);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(int, int, SUBTRACT);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(int, short, SUBTRACT);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(short, int, SUBTRACT);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(short, short, SUBTRACT);
    DECLARE_VARIABLE_VARIABLE_OPERATOR(SUBTRACT, atl::Subtract);
    DECLARE_SCALAR_VARIABLE_OPERATOR(long double, SUBTRACT, atl::SubtractConstant);
    DECLARE_SCALAR_VARIABLE_OPERATOR(double, SUBTRACT, atl::SubtractConstant);
    DECLARE_SCALAR_VARIABLE_OPERATOR(float, SUBTRACT, atl::SubtractConstant);
    DECLARE_VARIABLE_SCALAR_OPERATOR(long double, SUBTRACT, atl::ConstantSubtract);
    DECLARE_VARIABLE_SCALAR_OPERATOR(double, SUBTRACT, atl::ConstantSubtract);
    DECLARE_VARIABLE_SCALAR_OPERATOR(float, SUBTRACT, atl::ConstantSubtract);

    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(long double, long double, MULTIPLY);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(long double, double, MULTIPLY);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(double, long double, MULTIPLY);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(long double, float, MULTIPLY);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(float, long double, MULTIPLY);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(long double, long, MULTIPLY);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(long, long double, MULTIPLY);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(long double, int, MULTIPLY);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(int, long double, MULTIPLY);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(long double, short, MULTIPLY);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(short, long double, MULTIPLY);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(double, double, MULTIPLY);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(double, float, MULTIPLY);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(float, double, MULTIPLY);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(double, long, MULTIPLY);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(long, double, MULTIPLY);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(double, int, MULTIPLY);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(int, double, MULTIPLY);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(double, short, MULTIPLY);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(short, double, MULTIPLY);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(float, float, MULTIPLY);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(float, long, MULTIPLY);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(long, float, MULTIPLY);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(float, int, MULTIPLY);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(int, float, MULTIPLY);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(float, short, MULTIPLY);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(short, float, MULTIPLY);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(long, long, MULTIPLY);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(long, int, MULTIPLY);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(int, long, MULTIPLY);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(long, short, MULTIPLY);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(short, long, MULTIPLY);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(int, int, MULTIPLY);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(int, short, MULTIPLY);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(short, int, MULTIPLY);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(short, short, MULTIPLY);
    DECLARE_VARIABLE_VARIABLE_OPERATOR(MULTIPLY, atl::Multiply);
    DECLARE_SCALAR_VARIABLE_OPERATOR(long double, MULTIPLY, atl::MultiplyConstant);
    DECLARE_SCALAR_VARIABLE_OPERATOR(double, MULTIPLY, atl::MultiplyConstant);
    DECLARE_SCALAR_VARIABLE_OPERATOR(float, MULTIPLY, atl::MultiplyConstant);
    DECLARE_VARIABLE_SCALAR_OPERATOR(long double, MULTIPLY, atl::ConstantMultiply);
    DECLARE_VARIABLE_SCALAR_OPERATOR(double, MULTIPLY, atl::ConstantMultiply);
    DECLARE_VARIABLE_SCALAR_OPERATOR(float, MULTIPLY, atl::ConstantMultiply);

    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(long double, long double, DIVIDE);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(long double, double, DIVIDE);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(double, long double, DIVIDE);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(long double, float, DIVIDE);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(float, long double, DIVIDE);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(long double, long, DIVIDE);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(long, long double, DIVIDE);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(long double, int, DIVIDE);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(int, long double, DIVIDE);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(long double, short, DIVIDE);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(short, long double, DIVIDE);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(double, double, DIVIDE);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(double, float, DIVIDE);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(float, double, DIVIDE);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(double, long, DIVIDE);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(long, double, DIVIDE);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(double, int, DIVIDE);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(int, double, DIVIDE);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(double, short, DIVIDE);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(short, double, DIVIDE);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(float, float, DIVIDE);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(float, long, DIVIDE);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(long, float, DIVIDE);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(float, int, DIVIDE);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(int, float, DIVIDE);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(float, short, DIVIDE);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(short, float, DIVIDE);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(long, long, DIVIDE);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(long, int, DIVIDE);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(int, long, DIVIDE);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(long, short, DIVIDE);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(short, long, DIVIDE);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(int, int, DIVIDE);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(int, short, DIVIDE);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(short, int, DIVIDE);
    DECLARE_PRIMITIVE_PRIMITIVE_OPERATOR(short, short, DIVIDE);
    DECLARE_VARIABLE_VARIABLE_OPERATOR(DIVIDE, atl::Divide);
    DECLARE_SCALAR_VARIABLE_OPERATOR(long double, DIVIDE, atl::DivideConstant);
    DECLARE_SCALAR_VARIABLE_OPERATOR(double, DIVIDE, atl::DivideConstant);
    DECLARE_SCALAR_VARIABLE_OPERATOR(float, DIVIDE, atl::DivideConstant);
    DECLARE_VARIABLE_SCALAR_OPERATOR(long double, DIVIDE, atl::ConstantDivide);
    DECLARE_VARIABLE_SCALAR_OPERATOR(double, DIVIDE, atl::ConstantDivide);
    DECLARE_VARIABLE_SCALAR_OPERATOR(float, DIVIDE, atl::ConstantDivide);

//    DECLARE_VARIABLE_VARIABLE_OPERATOR(POW, atl::Pow);
//    DECLARE_SCALAR_VARIABLE_OPERATOR(long double, POW, atl::PowConstant);
//    DECLARE_SCALAR_VARIABLE_OPERATOR(double, POW, atl::PowConstant);
//    DECLARE_SCALAR_VARIABLE_OPERATOR(float, POW, atl::PowConstant);
//    DECLARE_VARIABLE_SCALAR_OPERATOR(long double, POW, atl::ConstantPow);
//    DECLARE_VARIABLE_SCALAR_OPERATOR(double, POW, atl::ConstantPow);
//    DECLARE_VARIABLE_SCALAR_OPERATOR(float, POW, atl::ConstantPow);
#endif

    template<typename T>
    struct PromoteUnaryOp {
    };


    //    DECLARE_PROMOTION(short, short, short);
    //    DECLARE_PROMOTION(int, int, int);
    //    DECLARE_PROMOTION(short, int, int);
    //    DECLARE_PROMOTION(int, short, int);
    //    DECLARE_PROMOTION(long, long, long);
    //    DECLARE_PROMOTION(long, short, long);
    //    DECLARE_PROMOTION(short, long, long);
    //    DECLARE_PROMOTION(long, int, long);
    //    DECLARE_PROMOTION(int, long, long);
    //
    //
    //
    //
    //    DECLARE_PROMOTION(float, float, float);
    //    DECLARE_PROMOTION(float, short, float);
    //    DECLARE_PROMOTION(short, float, float);
    //    DECLARE_PROMOTION(float, int, float);
    //    DECLARE_PROMOTION(int, float, float);
    //    DECLARE_PROMOTION(float, long, float);
    //    DECLARE_PROMOTION(long, float, float);
    //
    //    DECLARE_PROMOTION(double, double, double);
    //    DECLARE_PROMOTION(double, short, double);
    //    DECLARE_PROMOTION(short, double, double);
    //    DECLARE_PROMOTION(double, int, double);
    //    DECLARE_PROMOTION(int, double, double);
    //    DECLARE_PROMOTION(double, long, double);
    //    DECLARE_PROMOTION(long, double, double);
    //    DECLARE_PROMOTION(double, float, double);
    //    DECLARE_PROMOTION(float, double, double);
    //
    //
    //    DECLARE_PROMOTION(long double, long double, long double);
    //    DECLARE_PROMOTION(long double, short, long double);
    //    DECLARE_PROMOTION(short, long double, long double);
    //    DECLARE_PROMOTION(long double, int, long double);
    //    DECLARE_PROMOTION(int, long double, long double);
    //    DECLARE_PROMOTION(long double, long, long double);
    //    DECLARE_PROMOTION(long, long double, long double);
    //    DECLARE_PROMOTION(long double, float, long double);
    //    DECLARE_PROMOTION(float, long double, long double);
    //    DECLARE_PROMOTION(long double, double, long double);
    //    DECLARE_PROMOTION(double, long double, long double);
    //
    //    DECLARE_PROMOTION(et4ad::Variable<long double>, long double, et4ad::Variable<long double>);
    //    DECLARE_PROMOTION(long double, et4ad::Variable<long double>, et4ad::Variable<long double>);
    //    DECLARE_PROMOTION(et4ad::Variable<long double>, double, et4ad::Variable<long double>);
    //    DECLARE_PROMOTION(double, et4ad::Variable<long double>, et4ad::Variable<long double>);
    //    DECLARE_PROMOTION(et4ad::Variable<long double>, float, et4ad::Variable<long double>);
    //    DECLARE_PROMOTION(float, et4ad::Variable<long double>, et4ad::Variable<long double>);
    //    DECLARE_PROMOTION(et4ad::Variable<long double>, long, et4ad::Variable<long double>);
    //    DECLARE_PROMOTION(long , et4ad::Variable<long double>, et4ad::Variable<long double>);
    //    DECLARE_PROMOTION(et4ad::Variable<long double>, int, et4ad::Variable<long double>);
    //    DECLARE_PROMOTION(int, et4ad::Variable<long double>, et4ad::Variable<long double>);
    //    DECLARE_PROMOTION(et4ad::Variable<long double>, short, et4ad::Variable<long double>);
    //    DECLARE_PROMOTION(long double, et4ad::Variable<long double>, et4ad::Variable<long double>);

    /**
     * 
     * Traits for intrinsic sse types.
     * false and size =0 are default.
     *  
     */
    template<typename T1>
    struct intrinsic_trait {
        typedef T1 TYPE;
        static bool is_intrinsic;
        static size_t size;
    };

    template<typename T1>
    bool intrinsic_trait<T1>::is_intrinsic = false;

    template<typename T1>
    size_t intrinsic_trait<T1>::size = 0;

    //
#define DECLARE_INTRINSIC(A, B) template<> struct intrinsic_trait<A> { \
typedef A TYPE;   \
    static bool is_intrinsic; \
     static size_t size; \
}; \
     bool intrinsic_trait<A>::is_intrinsic = true; \
     size_t intrinsic_trait<A>::size = sizeof (A)*8;\

    DECLARE_INTRINSIC(long double, true);
    DECLARE_INTRINSIC(double, true);
    DECLARE_INTRINSIC(float, true);

    
    template<typename T1>
    struct IntrinsicBaseType {
        typedef T1 TYPE;
    };
    
#define DECLARE_INTRINSIC_BASE(A, B) template<> struct IntrinsicBaseType<A> { \
typedef B TYPE;   \
}; \

    DECLARE_INTRINSIC_BASE(atl::Variable<float>, float)
    DECLARE_INTRINSIC_BASE(atl::Variable<double>, double)
    DECLARE_INTRINSIC_BASE(atl::Variable<long double>, long double)
    
}



#endif	/* TRAITS_HPP */

