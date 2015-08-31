/* 
 * File:   Norm.hpp
 * Author: matthewsupernaw
 *
 * Created on December 17, 2014, 9:26 AM
 */

#ifndef NORM_HPP
#define	NORM_HPP

#include "Vector.hpp"
#include "../AutoDiff/Expression.hpp"
#include "VectorExpressionBase.hpp"
#include "../Traits/Type.hpp"

namespace atl {

    template<class REAL_T, class T, class V>//RHS is vector, array, or matrix
    struct VectorVariableNorm2 : public atl::ExpressionBase< REAL_T, VectorVariableNorm2<REAL_T, T, V> > {
        const V& v_m;
        mutable REAL_T value_m;
        mutable atl::Variable<REAL_T>  var;
        REAL_T coefficient;
        mutable bool computed;

        //for adjoint-mode accumulation 
        atl::Adjoint<REAL_T> entry_m;
        std::shared_ptr<uint32_t> adjoint_method_id;


        //for tag accumulation
        std::vector<REAL_T> gradient;
        atl::IDSet ids_m;

        VectorVariableNorm2(const atl::VectorExpression<T, V>& v) : v_m(v.Cast()), computed(false), adjoint_method_id(VariableIdGenerator::instance()->next()) {
            //
            //            uint32_t size = v.Size(0);
            //            //do everything here;
            //            if (atl::VariableTrait<T>::is_variable) {
            //
            //
            //                if (T::ACCUMULATOR == atl::REVERSE) {
            //                    adjoint_method_id = VariableIdGenerator::instance()->next();
            //                } else {
            //
            //                }
            //                if (atl::VariableTrait<T>::is_recording()) {
            //                    std::cout << "i'm recording\n" << std::flush;
            //                    exit(0);
            //                    value_m = 0;
            //                    coefficient = 1.0;
            //                    for (int i = 0; i < v_m.Size(0); i++) {
            //                        HandleExpression(v_m(i) * v_m(i));
            //                    }
            //
            //
            //
            //
            //                }
            //            } else {
            //
            //                //                for(int i=0; i< v_m.Size(0); i++){
            //                //                    value_m+=v_m(i)*v_m(i);
            //                //                    std::cout<<value_m<<"\n";
            //                //                }
            //
            //            }

        }

        ~VectorVariableNorm2() {

            if (atl::Variable<REAL_T>::ACCUMULATOR == atl::REVERSE) {
                atl::Variable<REAL_T>::gradient_stack.Recycle(adjoint_method_id);
            }

        }

        //        operator REAL_T()const {
        //            return this->GetValue();
        //        }
        //
        //        operator REAL_T() {
        //            return this->GetValue();
        //        }

        template<class TT>
        void HandleExpression(atl::Adjoint<REAL_T>& entry_m, REAL_T& coefficient, const ExpressionBase<REAL_T, TT>& exp)const {
            this->value_m += exp.GetValue() * exp.GetValue();
            (var+exp).PushStackEntry(entry_m,coefficient);
//            entry_m.coefficients.push_back(std::pair<uint32_t, REAL_T>(*adjoint_method_id, coefficient));
//            std::cout << "->" << *adjoint_method_id << " , " << coefficient << "\n";
//            exp.PushStackEntry(entry_m, coefficient);
            var.SetValue(this->value_m);
        }

        inline const REAL_T GetValue() const {

            if (!computed) {
                for (int i = 0; i < v_m.Size(0); i++) {
                    value_m += (v_m(i) * v_m(i)).GetValue();
                }
            }
            return this->value_m;
        }

        inline const uint32_t GetId() const {
            return 0;
        }

        /**
         * Compute or get the stored derivative with respect to the unique 
         * identifier id.
         * 
         * @param id
         * @param found
         * @return 
         */
        const REAL_T Derivative(const uint32_t &id, bool &found) const {
            return 1;
        }

        /**
         * Compute or get the stored derivative with respect to the unique 
         * identifier id.
         * 
         * @param id
         * @param found
         * @return 
         */
        inline const REAL_T Derivative(const uint32_t & id) const {
            return 1;
        }

#ifdef ET4AD_USE_SSE

        /**
         * Using sse operations, compute or get the stored derivative with respect to the unique 
         * identifier id.
         * 
         * @param id
         * @return __m128d type
         */
        inline const __m128d SSEDoubleDerivative(uint32_t* id) const {
            return Cast().SSEDoubleDerivative(id);
        }

        /**
         * Using sse operations, compute or get the stored derivative with respect to the unique 
         * identifier id.
         * 
         * @param id
         * @return __m128 type
         */
        inline const __m128 SSEFloatDerivative(uint32_t* id) const {
            return Cast().SSEFloatDerivative(id);
        }

#endif

        /**
         * Complex Step value. Used for complex step derivative computations.
         * Complex step differentiation is a technique that employs complex 
         * arithmetic to obtain the numerical value of the first derivative of 
         * a real valued analytic function of a real variable, avoiding the 
         * loss of precision inherent in traditional finite differences.
         * 
         * For Example:
         * 
         * F(x+ih) = cos(x+ih)
         * dF/dx+ih = F(x+ih).imag()/h
         * 
         * @param id
         * @param h
         * @return 
         */
        inline const std::complex<REAL_T> ComplexStepValue(const uint32_t & id, REAL_T h = REAL_T(0.00000000000001)) const {
            return 1;
        }

        /**
         * Push the ids for independent variables to a storage object.
         * 
         * @param storage
         */
        inline void PushIds(atl::IDSet & ids) const {

        }

        /**
         * Push expression statement into a vector in postfix order. 
         * @param storage
         */
        inline void PushStatements(std::vector<Statement<REAL_T> > &storage) const {

        }

        /**
         * Push partials into a gradient stack entry. Forward derivatives are
         * accumulated and used during the reverse mode differentiation method
         * in the class GradientStack.
         * @param entry
         * @param coefficient
         */
        inline void PushStackEntry(Adjoint<REAL_T>& entry, REAL_T coefficient = 1.0) const {
//            var.SetValue(value_m);
//            atl::Variable<REAL_T> v;
            for (int i = 0; i < v_m.Size(0); i++) {
//                (var +v_m(i) * v_m(i)).PushStackEntry(entry,coefficient);
                HandleExpression(entry, coefficient, v_m(i) * v_m(i));
            }
//            v.PushStackEntry(entry);
//            value_m = var.GetValue();
//            var.PushStackEntry(entry,coefficient);
            this->computed = true;
        }

        /**
         * Push entries for other dependent variables. This is used for reverse 
         * mode derivative calculations. Each dependent variable maintains a 
         * list of StackEntry's and is only used if a gradient calculation is 
         * required. The element in the vector is a pair, that contains a vector
         * of entries and a index value representing the end of ownership, 
         * just in case the entered vector has been appended after this function
         * has been called.
         * 
         * @param entries
         */
        inline void PushStackEntry(std::vector<std::shared_ptr<std::vector<Adjoint<REAL_T> > > >& entries) const {

        }

        inline void PushStackEntry(std::shared_ptr<atl::GradientStack<REAL_T> >& gs) const {

        }

        /**
         * Used for merging dependent variable info in the encapsulated version
         * of the adjoint algorithm. 
         * @param exp
         * @return 
         */
        inline void MergeGradientStructures(atl::GradientStack<REAL_T>& gs) const {

        }





    };

    template<class T, class V>
    const inline VectorVariableNorm2<typename T::BASE_TYPE, T, V > Norm2(const atl::VectorExpression<T, V> & v) {
        return VectorVariableNorm2<typename T::BASE_TYPE, T, V > (v.Cast());
    }

    //#define MAKE_PRIMITIVE_NORM2(T)  \
//template<class V>                 \
//    const inline VectorVariableNorm2< T, T, V > Norm2(const atl::VectorExpression<T, V> & v) { \
//        return VectorVariableNorm2< T, T, V > (v.Cast()); \
//    }\
//
    //    MAKE_PRIMITIVE_NORM2(long double);
    //    MAKE_PRIMITIVE_NORM2(double);
    //    MAKE_PRIMITIVE_NORM2(float);
}



#endif	/* NORM_HPP */

