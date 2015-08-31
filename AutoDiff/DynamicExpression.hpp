/* 
 * File:   DynamicExpression.hpp
 * Author: matthewsupernaw
 *
 * Created on June 12, 2015, 11:00 AM
 */

#ifndef DYNAMICEXPRESSION_HPP
#define	DYNAMICEXPRESSION_HPP

#include <cmath>
#include <string>
#include <sstream>
#include <iostream>
#include <memory>
#include "Statement.hpp"
#include <atomic>
#include "../third_party/clfmalloc.h"
#include "../third_party/khash.hpp"
#include "GradientStructure.hpp"
#include "VariableInfo.hpp"

namespace atl {
    //
    //    template<class T, size_t initSize>
    //    class MemPool {
    //        std::vector<T> poolBuffer; // The memory pool
    //        std::vector<unsigned short> freeList; // Ring-buffer (indices of free items)
    //        std::atomic<unsigned short> nHead;
    //        std::atomic<unsigned short> nTail; // Ring-buffer management
    //        std::atomic<int> nCount; // Number of elements in ring-buffer
    //        std::mutex mutex_m;
    //        
    //        MemPool<T, initSize> *nextPool; // For expanding memory pool
    //        public:
    //        MemPool(){
    //            nCount = initSize-1;
    //            poolBuffer.resize(initSize);
    //            nHead =0;
    //            nTail =nCount;
    //            
    //        }
    //        
    //        T* alloc() {
    //            // Check the queue size.  If it's zero (or less) we need to pass on
    //            // to the next pool and/or allocate a new one.
    //            if (nCount <= 0) {
    //                return alloc_extend();
    //            }
    //
    //            int count = (nCount--);
    //            if (count <= 0) {
    //                T *mem = alloc_extend();
    //                (nCount++); // undo
    //                return mem;
    //            }
    //
    //            // We are guaranteed that there is at least 1 element in the list for us.
    //            // This will overflow naturally to achieve modulo by 65536.  You can only
    //            // deal with queue sizes that are a power of 2.  If you want 32768 values,
    //            // for example, you must do this: head &= 0x7fff;
    //            unsigned short head = (nHead++);
    //
    //            // Spin until the element is valid (use a reference)
    //            unsigned short & idx = freeList[head];
    ////            while (idx == 0xffff);
    //
    //            // Grab the pool item, and blitz the index from the queue
    //            T * mem = &poolBuffer[idx];
    //            idx = 0xffff;
    //
    //            return mem;
    //        }
    //
    //        T * alloc_extend() {
    //            if (nextPool == NULL) {
    //                mutex_m.lock();
    //                if (nextPool == NULL) nextPool = new MemPool<T,initSize>();
    //                mutex_m.unlock();
    //                if (nextPool == NULL) return NULL;
    //            }
    //            return nextPool->alloc();
    //        }
    //
    //        void free(T* mem) {
    //            // Find the right pool to free from.
    //            if (mem < &poolBuffer.front() || mem > &poolBuffer.back()) {
    //                if (nextPool) nextPool->free(mem);
    //                return;
    //            }
    //
    //            // You might want to maintain a bitset that indicates whether the memory has
    //            // actually been allocated so you don't corrupt your pool here, but I won't
    //            // do that in this example...
    //
    //            // Work out the index.  Hope my pointer arithmetic is correct here.
    //            unsigned short idx = (unsigned short) (mem - &poolBuffer.front());
    //
    //            // Push index back onto the queue.  As with alloc(), you might want to
    //            // use a mask on the tail to achieve modulo.
    //            int tail = (nTail++);
    //            freeList[tail] = idx;
    //
    //            // Don't need to check the list size.  We assume everything is sane. =)
    //            (nCount++);
    //        }
    //
    //    };

    template<class T>
    class MemoryPool {
        std::vector<T> pool; //actual heap of objects
        std::vector<T* > free_list; //available objects
        long long index;
        size_t size_m, next_size_m;
        MemoryPool<T >* next_m;

        inline void* extend() {
            if (next_m == NULL) {
                this->next_m = new MemoryPool<T> (this->next_size_m);
            }

            return next_m->malloc();
        }

    public:

        MemoryPool(uint32_t size) : next_size_m(size), size_m(size), pool(size), free_list(size), index(size - 1), next_m(NULL) {

            for (int i = 0; i < size; i++) {
                free_list[i] = (&pool[i]);
            }
        }

        ~MemoryPool() {

            if (next_m != NULL) {
                delete next_m;
            }
        }

        void SetResize(uint32_t resize) {
            this->next_size_m = resize;
        }

        inline void* malloc() {
            if (index <= 0) {
                return this->extend();
            }
            return free_list[index--];
        }

        inline void free(void* ptr) {

            if (ptr < &pool[0] || ptr > &pool[size_m - 1]) {
                if (next_m != NULL) {
                    next_m->free(ptr);
                } else {
                    std::cout << "pointer not freed...\n";
                }

            } else {
                if (ptr)
                    free_list[++index] = (T*) ptr;
            }

        }

    };


    //    enum Operation {
    //        MINUS = 0,
    //        PLUS,
    //        MULTIPLY,
    //        DIVIDE,
    //        SIN,
    //        COS,
    //        TAN,
    //        ASIN,
    //        ACOS,
    //        ATAN,
    //        ATAN2, //atan(adnumber,adnumber)
    //        ATAN3, //atan(T,adnumber)
    //        ATAN4, //atan(adnumber,T)
    //        SQRT,
    //        POW, //pow(adnumber,adnumber)
    //        POW1, //pow(T,adnumber)
    //        POW2, //pow(adnumber,T)
    //        LOG,
    //        LOG10,
    //        EXP,
    //        SINH,
    //        COSH,
    //        TANH,
    //        ABS,
    //        FABS,
    //        FLOOR,
    //        CONSTANT,
    //        VARIABLE,
    //        NONE
    //    };

    template<class T>
    class DynamicExpression;

    template<class T>
    static T EvalConstant(DynamicExpression<T>* exp);

    template<class T>
    static T EvalVariable(DynamicExpression<T>* exp);

    template<class T>
    static T EvalPlus(DynamicExpression<T>* exp);

    template<class T>
    static T EvalMinus(DynamicExpression<T>* exp);

    template<class T>
    static T EvalMultiply(DynamicExpression<T>* exp);

    template<class T>
    static T EvalDivide(DynamicExpression<T>* exp);

    //    template<class T>
    //    T (*EvalFunctions[4]) = {EvalConstant<T>, EvalVariable<T>,EvalPlus<T>,EvalMinus<T>};
    //    
#warning need to make DynamicExpression::GetValue and DynamicExpression::PushHessianEntry function pointers.

    /**
     * class DynamicExpression
     * 
     * Used when higher order derivative information is required. Since 
     * Expression templates can easily become large and cause stack overflow
     * when differentiated, dynamic memory is used via the DynamicExpression 
     * class.
     * 
     */
    template<class T>
    class DynamicExpression {
    public:
        static MemoryPool<DynamicExpression<T> > pool_m;

        typedef std::vector<std::pair<VariableInfo<T>*, T> > HessianInfo;
        T value_m;
        DynamicExpression<T>* left_m;
        DynamicExpression<T>* right_m;
        Operation op_m;
        atl::VariableInfo<T>* info;

        DynamicExpression()
        : right_m(NULL),
        left_m(NULL),
        op_m(CONSTANT),
        value_m(T(0.0)),
        info(NULL) {


        }

        DynamicExpression(const T &value,
                const Operation &op)
        : op_m(op),
        info(NULL),
        right_m(NULL),
        left_m(NULL),
        value_m(value) {

        }

        DynamicExpression(const T &value,
                atl::VariableInfo<T>* info,
                const Operation &op)
        : op_m(op),
        info(info),
        right_m(NULL),
        left_m(NULL),
        value_m(value) {

        }

        DynamicExpression(const T &value,
                atl::VariableInfo<T>* info,
                const Operation &op,
                DynamicExpression<T>* left,
                DynamicExpression<T>* right)
        : op_m(op),
        right_m(right),
        left_m(left),
        info(info),
        value_m(value) {

        }

        DynamicExpression(const Operation &op,
                DynamicExpression<T>* left,
                DynamicExpression<T>* right)
        : op_m(op),
        right_m(right),
        left_m(left),
        info(NULL),
        value_m(0) {

        }

        DynamicExpression(const Operation &op,
                DynamicExpression<T>* left)
        : op_m(op),
        left_m(left),
        right_m(NULL),
        info(NULL),
        value_m(0) {

        }

        ~DynamicExpression() {
            if (left_m != NULL) {
                delete left_m;
            }

            if (right_m != NULL) {
                delete right_m;
            }
        }

        void* operator new(size_t size) {
            return malloc(size); //DynamicExpression<T>::pool_m.malloc();
        }

        void operator delete(void* ptr) {
            free(ptr); //DynamicExpression<T>::pool_m.free((DynamicExpression<T>*)ptr);
        }

        T GetValue() const {

            switch (op_m) {
                case CONSTANT:
                    return value_m;
                case VARIABLE:
                    return value_m;
                case PLUS:
                    return left_m->GetValue() + right_m->GetValue();
                case MINUS:
                    return left_m->GetValue() - right_m->GetValue();
                case MULTIPLY:
                    return left_m->GetValue() * right_m->GetValue();
                case DIVIDE:
                    return left_m->GetValue() / right_m->GetValue();
                default:
                    return 0.0;
            }

        }

        void PushHessianEntry(HessianInfo& hessian_row, T coefficient = 1.0) {
            switch (op_m) {

                case CONSTANT:
                    break;
                case VARIABLE:
//                    std::cout << "pushing " << info->id << ", " << coefficient << "\n";
                    //                    hessian_row[info] =  coefficient;
//                    if (coefficient != 0.0)
                        hessian_row.push_back(std::make_pair(info, coefficient));
                    break;
                case PLUS:
                    this->left_m->PushHessianEntry(hessian_row, coefficient);
                    this->right_m->PushHessianEntry(hessian_row, coefficient);
                    break;
                case MINUS:
                    this->left_m->PushHessianEntry(hessian_row, coefficient);
                    this->right_m->PushHessianEntry(hessian_row, -1.0 * coefficient);
                    break;
                case MULTIPLY:
                    this->left_m->PushHessianEntry(hessian_row, coefficient * this->GetRight()->GetValue());
                    this->right_m->PushHessianEntry(hessian_row, coefficient * this->GetLeft()->GetValue());
                    break;
                case DIVIDE:
                    this->left_m->PushHessianEntry(hessian_row, coefficient * (1.0 / GetRight()->GetValue()));
                    this->right_m->PushHessianEntry(hessian_row, -1.0*coefficient * GetValue()* (1.0 / GetRight()->GetValue()));
                    break;
                default:
                    std::cout << "Operation not yet implemented!\n";
                    break;

            }
        }

        DynamicExpression<T>* GetLeft() const {
            return left_m;
        }

        void SetLeft(DynamicExpression<T>* left) {
            this->left_m = left;
        }

        DynamicExpression<T>* GetRight() const {
            return right_m;
        }

        void SetRight(DynamicExpression<T>* right) {
            this->right_m = right;
        }

        Operation GetOp() const {
            return op_m;
        }

        void SetOp(Operation op) {
            this->op_m = op;
        }

        /*!
         * Represent this expression as a string. ADNumbers are represented
         * in wkt format by default unless latex flag is set to true.
         * Constants are represented by value.
         *
         */
        std::string ToString(bool latex = false) {

            //            if (Size() > 1000) {
            //                return std::string("Error: ToString...expression too long.");
            //            }
            //
            std::stringstream ss;
            std::stringstream temp;

            std::string l, r;



            if (GetRight() != NULL) {
                r = GetRight()->ToString(latex);
            }

            if (GetLeft() != NULL) {
                l = GetLeft()->ToString(latex);
            }


            temp.str("");



            switch (op_m) {
                case CONSTANT:
                    if (latex) {
                        ss << GetValue();
                    } else {
                        ss << "CONST[" << GetValue() << "]";
                    }
                    break;
                case VARIABLE:

                    ss << "VAR[" << GetValue() << ",ID[" << info->id << "]" << "]";

                    break;
                case MINUS:
                    ss << "(" << l << " - " << r << ")";
                    break;
                case PLUS:
                    ss << "(" << l << " + " << r << ")";
                    break;
                case DIVIDE:
                    if (latex) {
                        ss << "\\frac{" << l << "}{" << r << "}";
                    } else {
                        ss << l << " / " << r;
                    }
                    break;
                case MULTIPLY:

                    if (latex && GetRight()->GetOp() == POW
                            && GetRight()->GetRight()->GetOp() == CONSTANT
                            && GetRight()->GetRight()->GetValue() == T(-1)) {
                        ss << "\\frac{" << l << "}{" << GetRight()->GetLeft()->ToString(latex) << "}";
                    } else {

                        ss << "(" << l << " * " << r << ")";
                    }
                    break;
                case SIN:
                    ss << "sin(" << l << ")";
                    break;
                case COS:
                    ss << "cos(" << l << ")";
                    break;
                case TAN:
                    ss << "tan(" << l << ")";
                    break;
                case ASIN:
                    ss << "asin(" << l << ")";
                    break;
                case ACOS:
                    ss << "acos(" << l << ")";
                    break;
                case ATAN:
                    ss << "atan(" << l << ")";
                    break;
                case ATAN2:
                    ss << "atan(" << l << "," << r << ")";
                    break;
                case ATAN3:
                    ss << "atan(" << l << "," << r << ")";
                    break;
                case ATAN4:
                    ss << "atan(" << l << "," << r << ")";
                    break;
                case SQRT:
                    ss << "sqrt(" << l << ")";
                    break;
                case POW:


                case POW1:

                case POW2:
                    if (latex) {
                        ss << l << "^{" << r << "}";
                    } else {
                        ss << "pow(" << l << "," << r << ")";
                    }
                    break;
                case LOG:
                    if (latex) {
                        ss << "ln(" << l << ")";
                    } else {
                        ss << "log(" << l << ")";
                    }
                    break;
                case LOG10:
                    ss << "log10(" << l << ")";
                    break;
                case EXP:
                    ss << "exp(" << l << ")";
                    break;
                case SINH:
                    ss << "sinh(" << l << ")";
                    break;
                case COSH:
                    ss << "cosh(" << l << ")";
                    break;
                case TANH:
                    ss << "tanh(" << l << ")";
                    break;
                case FABS:
                    ss << "fabs(" << l << ")";
                case FLOOR:
                    ss << "floor(" << l << ")";
                case NONE:
                    break;
                default:
                    break;
            }

            //ss << ")";

            return ss.str();

        }


    };
    template<class T>
    MemoryPool<DynamicExpression<T> > DynamicExpression<T>::pool_m(10000);


    //    enum Operation {
    //        MINUS = 0,
    //        PLUS,
    //        MULTIPLY,
    //        DIVIDE,
    //        SIN,
    //        COS,
    //        TAN,
    //        ASIN,
    //        ACOS,
    //        ATAN,
    //        ATAN2, //atan(adnumber,adnumber)
    //        ATAN3, //atan(T,adnumber)
    //        ATAN4, //atan(adnumber,T)
    //        SQRT,
    //        POW, //pow(adnumber,adnumber)
    //        POW1, //pow(T,adnumber)
    //        POW2, //pow(adnumber,T)
    //        LOG,
    //        LOG10,
    //        EXP,
    //        SINH,
    //        COSH,
    //        TANH,
    //        ABS,
    //        FABS,
    //        FLOOR,
    //        CONSTANT,
    //        VARIABLE,
    //        NONE
    //    };

    //        template<class T>
    //        class DynamicExpression {
    //        public:
    //            typedef khmap_t<VariableInfo<T>*, T> HessianInfo;
    //            T value_m;
    //            DynamicExpression<T>* left_m;
    //            DynamicExpression<T>* right_m;
    //            unsigned long id_m;
    //            Operation op_m;
    //            std::string name_m;
    //    
    //    
    //            static std::atomic<uint32_t> counter;
    //            static std::vector<DynamicExpression<T> > pool;
    //    
    //            DynamicExpression()
    //            : right_m(NULL),
    //            left_m(NULL),
    //            op_m(CONSTANT),
    //            id_m(0),
    //            value_m(T(0.0)) {
    //    
    //    
    //            }
    //    
    //            DynamicExpression(const T &value, const unsigned long &id, const Operation &op, DynamicExpression<T>* left, DynamicExpression<T>* right)
    //            : op_m(op),
    //            id_m(id),
    //            right_m(right),
    //            left_m(left),
    //            value_m(value) {
    //    
    //            }
    //    
    //            ~DynamicExpression() {
    //    
    //                //            if (this->left_m != NULL) {
    //                //                delete this->left_m;
    //                //            }
    //                //
    //                //            if (this->right_m != NULL) {
    //                //                delete this->right_m;
    //                //            }
    //            }
    //    
    //            static inline DynamicExpression<T>* NextStatic() {
    //    
    //                uint32_t c = counter++;
    //                if (c >= DynamicExpression<T>::pool.size()) {
    //                    std::cout << "Dynamic Expression pool exhausted!!\n";
    //                    exit(0);
    //                }
    //                return &DynamicExpression<T>::pool[counter++];
    //            }
    //    
    //            static void Reset() {
    //                counter = 0;
    //            }
    //    
    //            void* operator new(size_t size) {
    //                return malloc(size);
    //            }
    //    
    //            void operator delete(void* ptr) {
    //                free(ptr);
    //            }
    //            //        
    //    
    //            bool HasId(const uint32_t &id) {
    //                //     std::cout << this->id_ << " ?= " << id << "\n";
    //                if (this->id_m == id) {
    //                    return true;
    //                }
    //                if (this->left_m) {
    //                    if (this->left_m->HasId(id)) {
    //                        return true;
    //                    }
    //                }
    //    
    //                if (this->right_m) {
    //                    if (this->right_m->HasId(id)) {
    //                        return true;
    //                    }
    //                }
    //    
    //                return false;
    //            }
    //    
    //    
    //    
    //            //         void PushStackEntry(atl::Adjoint<T>* entry) {}
    //    
    //            //            switch (op_m) {
    //            //
    //            //                case CONSTANT:
    //            //
    //            //
    //            //                case VARIABLE:
    //            //
    //            //                case MINUS:
    //            //
    //            //
    //            //
    //            //                case PLUS:
    //            //
    //            //
    //            //
    //            //                case DIVIDE:
    //            //
    //            //                case MULTIPLY:
    //            //
    //            //                case SIN:
    //            //
    //            //                case COS:
    //            //
    //            //                case TAN:
    //            //
    //            //                case ASIN:
    //            //
    //            //                case ACOS:
    //            //
    //            //
    //            //                case ATAN:
    //            //
    //            //                case ATAN2:
    //            //
    //            //                case ATAN3:
    //            //
    //            //                    //can be removed.
    //            //                    break;
    //            //
    //            //                case ATAN4:
    //            //                    break;
    //            //                case SQRT:
    //            //
    //            //                case POW:
    //            //                case LOG:
    //            //
    //            //                case LOG10:
    //            //
    //            //                case EXP:
    //            //
    //            //                case SINH:
    //            //
    //            //                case COSH:
    //            //
    //            //                case TANH:
    //            //                    //f(x) = tanh(x)
    //            //                    //f'(x) =1- tanh(x)*tanh(x)
    //            //
    //            //
    //            //
    //            //                case FABS:
    //            //
    //            //                case FLOOR:
    //            //
    //            //                case NONE://shouldn't happen.
    //            //
    //            //                default:
    //            //                    break;
    //            //            }
    //            //
    //    
    //    
    //    
    //            //        }
    //    
    //            /*!
    //             * Builds a expression tree representing the derivative with respect to 
    //             * some ADNumber via its id.(reverse mode) 
    //             * 
    //             * @return Expression<T>
    //             */
    //            DynamicExpression<T>* Differentiate(const uint32_t &id) {
    //                //#warning need to check partial derivatives....
    //    
    //                DynamicExpression<T>* ret = new DynamicExpression<T > ();
    //    
    //                //            std::cout<<this->ToString()<<std::endl;
    //                switch (op_m) {
    //    
    //                    case CONSTANT:
    //                        //f(x) = C
    //                        //f'(x) = 0
    //    
    //                        ret->op_m = CONSTANT;
    //                        ret->value_m = T(0); //this->value_;
    //    
    //    
    //                        return ret;
    //    
    //                    case VARIABLE:
    //                        if (this->id_m == id) {
    //                            //f(x) = x
    //                            //f'(x) = 1
    //    
    //                            ret->op_m = CONSTANT;
    //                            ret->value_m = T(1.0);
    //    
    //    
    //                            return ret;
    //                        } else {//constant
    //                            //f(x) = C
    //                            //f'(x) = 0
    //                            ret->op_m = CONSTANT;
    //                            ret->value_m = T(0.0);
    //                            return ret;
    //                        }
    //                    case MINUS:
    //    
    //                        //f(x) = g(x) - h(x)
    //                        //f'(x) = g'(x) - h'(x)
    //    
    //                        ret->op_m = MINUS;
    //                        if (this->left_m) {
    //                            ret->left_m = this->left_m->Differentiate(id);
    //    
    //                        }
    //    
    //                        if (this->right_m) {
    //                            ret->right_m = this->right_m->Differentiate(id);
    //                        }
    //    
    //                        return ret;
    //    
    //                    case PLUS:
    //    
    //                        //f(x) = g(x) + h(x)
    //                        //f'(x) = g'(x) + h'(x)
    //    
    //                        ret->op_m = PLUS;
    //                        if (this->left_m) {
    //                            ret->left_m = this->left_m->Differentiate(id);
    //                        }
    //    
    //                        if (this->right_m) {
    //                            ret->right_m = this->right_m->Differentiate(id);
    //                        }
    //    
    //    
    //                        return ret;
    //    
    //                    case DIVIDE:
    //    
    //                        //f(x) = g(x)/h(x);
    //                        //f'(x) = (g'(x)h(x) - g(x)h'(x))/h(x)^2
    //    
    //    
    //                        ret->op_m = DIVIDE;
    //    
    //                        ret->left_m = new DynamicExpression<T > (); //g'(x)h(x) - g(x)h'(x)
    //                        ret->left_m->op_m = MINUS;
    //    
    //    
    //                        ret->left_m->left_m = new DynamicExpression<T > (); //g'(x)h(x)
    //                        ret->left_m->left_m->op_m = MULTIPLY;
    //                        if (this->left_m) {
    //                            ret->left_m->left_m->left_m = this->left_m->Differentiate(id);
    //                        }
    //                        ret->left_m->left_m->right_m = this->right_m; //->Clone();
    //    
    //                        ret->left_m->right_m = new DynamicExpression<T > (); //g(x)h'(x)
    //                        ret->left_m->right_m->op_m = MULTIPLY;
    //                        ret->left_m->right_m->left_m = this->left_m; //->Clone();
    //                        if (this->right_m) {
    //                            ret->left_m->right_m->right_m = this->right_m->Differentiate(id);
    //                        }
    //    
    //    
    //                        ret->right_m = new DynamicExpression<T > ();
    //                        ret->right_m->op_m = MULTIPLY;
    //                        ret->right_m->left_m = this->right_m; //->Clone();
    //                        ret->right_m->right_m = this->right_m; //->Clone();
    //    
    //    
    //                        return ret;
    //    
    //                    case MULTIPLY:
    //                        //f(x) = g(x)h(x);
    //                        //f'(x) = g'(x)h(x) + g(x)h'(x)
    //    
    //                        if (this->left_m->op_m == CONSTANT
    //                                && this->right_m->op_m != CONSTANT) {
    //                            ret->op_m = MULTIPLY;
    //                            if (this->left_m) {
    //                                ret->left_m = this->left_m; //->Clone();
    //                            }
    //                            if (this->right_m) {
    //                                ret->right_m = this->right_m->Differentiate(id);
    //                            }
    //    
    //    
    //                        } else if (this->right_m->op_m == CONSTANT
    //                                && this->left_m->op_m != CONSTANT) {
    //                            ret->op_m = MULTIPLY;
    //                            if (this->left_m) {
    //                                ret->left_m = this->left_m->Differentiate(id);
    //                            }
    //                            if (this->right_m) {
    //                                ret->right_m = this->right_m; //->Clone();
    //                            }
    //                        } else {
    //    
    //    
    //    
    //                            ret->op_m = PLUS;
    //    
    //                            ret->left_m = new DynamicExpression<T > ();
    //                            ret->left_m->op_m = MULTIPLY;
    //    
    //                            ret->left_m->right_m = this->right_m; //->Clone();
    //    
    //                            if (this->right_m != NULL) {
    //                                ret->left_m->left_m = this->left_m->Differentiate(id);
    //                            }
    //                            ret->right_m = new DynamicExpression<T > ();
    //                            ret->right_m->op_m = MULTIPLY;
    //    
    //                            ret->right_m->left_m = this->left_m; //->Clone();
    //                            if (this->left_m != NULL) {
    //                                ret->right_m->right_m = this->right_m->Differentiate(id);
    //                            }
    //    
    //    
    //    
    //                        }
    //                        return ret;
    //    
    //                    case SIN:
    //    
    //                        if (this->left_m->HasId(id)) {
    //                            //f'(x) = cos(x)
    //    
    //                            ret->op_m = MULTIPLY;
    //                            ret->left_m = this->left_m->Differentiate(id);
    //                            ret->right_m = new DynamicExpression<T > ();
    //                            ret->right_m->op_m = COS;
    //                            ret->right_m->left_m = this->left_m; //->Clone();
    //    
    //                            return ret;
    //                        } else {
    //                            ret->op_m = CONSTANT;
    //                            ret->value_m = T(0.0);
    //    
    //    
    //                            return ret;
    //                        }
    //    
    //                    case COS:
    //                        if (this->left_m->HasId(id)) {
    //                            //f'(x) = -sin(x)
    //    
    //                            ret->op_m = MULTIPLY;
    //    
    //    
    //                            ret->left_m = this->left_m->Differentiate(id);
    //                            ret->right_m = new DynamicExpression<T > ();
    //    
    //                            ret->right_m->op_m = MULTIPLY;
    //                            ret->right_m->left_m = new DynamicExpression<T > ();
    //                            ret->right_m->left_m->op_m = CONSTANT;
    //                            ret->right_m->left_m->value_m = T(-1.0);
    //    
    //                            ret->right_m->right_m = new DynamicExpression<T > ();
    //                            ret->right_m->right_m->op_m = SIN;
    //                            ret->right_m->right_m->left_m = this->left_m; //->Clone();
    //    
    //    
    //                            return ret;
    //    
    //                        } else {
    //                            ret->op_m = CONSTANT;
    //                            ret->value_m = T(0.0);
    //    
    //    
    //                            return ret;
    //                        }
    //                    case TAN:
    //                        if (this->left_m->HasId(id)) {
    //                            //f'(x) = 1/cos(x)
    //    
    //                            ret->op_m = MULTIPLY;
    //                            ret->left_m = this->left_m->Differentiate(id);
    //    
    //    
    //                            ret->right_m = new DynamicExpression<T > ();
    //                            ret->right_m->op_m = MULTIPLY;
    //    
    //                            ret->right_m->left_m = new DynamicExpression<T > ();
    //                            ret->right_m->left_m->op_m = DIVIDE;
    //    
    //    
    //                            ret->right_m->left_m->left_m = new DynamicExpression<T > ();
    //                            ret->right_m->left_m->left_m->op_m = CONSTANT;
    //                            ret->right_m->left_m->left_m->value_m = T(1.0);
    //    
    //    
    //                            ret->right_m->left_m->right_m = new DynamicExpression<T > ();
    //                            ret->right_m->left_m->right_m->op_m = COS;
    //                            ret->right_m->left_m->right_m->left_m = this->left_m; //->Clone();
    //    
    //    
    //                            ret->right_m->right_m = new DynamicExpression<T > ();
    //                            ret->right_m->right_m->op_m = DIVIDE;
    //    
    //    
    //                            ret->right_m->right_m->left_m = new DynamicExpression<T > ();
    //                            ret->right_m->right_m->left_m->op_m = CONSTANT;
    //                            ret->right_m->right_m->left_m->value_m = T(1.0);
    //    
    //    
    //                            ret->right_m->right_m->right_m = new DynamicExpression<T > ();
    //                            ret->right_m->right_m->right_m->op_m = COS;
    //                            ret->right_m->right_m->right_m->left_m = this->left_m; //->Clone();
    //    
    //    
    //                            return ret;
    //                        } else {
    //                            ret->op_m = CONSTANT;
    //                            ret->value_m = T(0.0);
    //    
    //    
    //                            return ret;
    //                        }
    //                    case ASIN:
    //    
    //                        if (this->left_m->HasId(id)) {
    //                            //f(x) = asin(x)
    //                            //f'(x) = 1/(2 sqrt(1-x^2)= 1/(pow((1-pow(x,2)),0.5)
    //    
    //                            ret->op_m = MULTIPLY;
    //                            ret->left_m = this->left_m->Differentiate(id);
    //    
    //    
    //                            ret->right_m = new DynamicExpression<T > ();
    //                            ret->right_m->op_m = DIVIDE;
    //    
    //                            ret->right_m->left_m = new DynamicExpression<T > ();
    //                            ret->right_m->left_m->op_m = CONSTANT;
    //                            ret->right_m->left_m ->value_m = T(1.0);
    //    
    //                            ret->right_m->right_m = new DynamicExpression<T > ();
    //                            ret->right_m->right_m->op_m = POW;
    //    
    //                            ret->right_m->right_m->left_m = new DynamicExpression<T > ();
    //                            ret->right_m->right_m->left_m->op_m = MINUS;
    //    
    //                            ret->right_m->right_m->left_m->left_m = new DynamicExpression<T > ();
    //                            ret->right_m->right_m->left_m->left_m->op_m = CONSTANT;
    //                            ret->right_m->right_m->left_m->left_m->value_m = T(1.0);
    //    
    //                            ret->right_m->right_m->left_m->right_m = new DynamicExpression<T > ();
    //                            ret->right_m->right_m->left_m->right_m->op_m = POW;
    //                            ret->right_m->right_m->left_m->right_m->left_m = this->left_m; //->Clone();
    //    
    //                            ret->right_m->right_m->left_m->right_m->right_m = new DynamicExpression<T > ();
    //                            ret->right_m->right_m->left_m->right_m->right_m->op_m = CONSTANT;
    //                            ret->right_m->right_m->left_m->right_m->right_m->value_m = T(2.0);
    //    
    //                            ret->right_m->right_m->right_m = new DynamicExpression<T > ();
    //                            ret->right_m->right_m->right_m->op_m = CONSTANT;
    //                            ret->right_m->right_m->right_m->value_m = T(0.5);
    //    
    //                            return ret;
    //                        } else {
    //                            ret->op_m = CONSTANT;
    //                            ret->value_m = T(0.0);
    //    
    //    
    //                            return ret;
    //                        }
    //                    case ACOS:
    //    
    //                        if (this->left_m->HasId(id)) {
    //                            //f(x) = acos(x)
    //                            //f'(x) = -1/(sqrt(1-x^2) = -1/(pow((1-pow(x,2)),0.5)
    //                            //-1/sqrt(1-x^2)
    //                            ret->op_m = MULTIPLY;
    //                            ret->left_m = new DynamicExpression<T > ();
    //                            ret->left_m->op_m = MULTIPLY;
    //                            ret->left_m->left_m = new DynamicExpression<T > ();
    //    
    //                            ret->left_m->left_m->op_m = CONSTANT;
    //                            ret->left_m->left_m->value_m = T(-1.0);
    //    
    //    
    //                            ret->left_m->right_m = this->left_m->Differentiate(id);
    //    
    //                            ret->right_m = new DynamicExpression<T > ();
    //                            ret->right_m->op_m = DIVIDE;
    //    
    //                            ret->right_m->left_m = new DynamicExpression<T > ();
    //                            ret->right_m->left_m->op_m = CONSTANT;
    //                            ret->right_m->left_m ->value_m = T(1.0);
    //    
    //                            ret->right_m->right_m = new DynamicExpression<T > ();
    //                            ret->right_m->right_m->op_m = POW;
    //    
    //                            ret->right_m->right_m->left_m = new DynamicExpression<T > ();
    //                            ret->right_m->right_m->left_m->op_m = MINUS;
    //    
    //                            ret->right_m->right_m->left_m->left_m = new DynamicExpression<T > ();
    //                            ret->right_m->right_m->left_m->left_m->op_m = CONSTANT;
    //                            ret->right_m->right_m->left_m->left_m->value_m = T(1.0);
    //    
    //                            ret->right_m->right_m->left_m->right_m = new DynamicExpression<T > ();
    //                            ret->right_m->right_m->left_m->right_m->op_m = POW;
    //                            ret->right_m->right_m->left_m->right_m->left_m = this->left_m; //->Clone();
    //    
    //                            ret->right_m->right_m->left_m->right_m->right_m = new DynamicExpression<T > ();
    //                            ret->right_m->right_m->left_m->right_m->right_m->op_m = CONSTANT;
    //                            ret->right_m->right_m->left_m->right_m->right_m->value_m = T(2.0);
    //    
    //                            ret->right_m->right_m->right_m = new DynamicExpression<T > ();
    //                            ret->right_m->right_m->right_m->op_m = CONSTANT;
    //                            ret->right_m->right_m->right_m->value_m = T(0.5);
    //    
    //                            return ret;
    //                        } else {
    //                            ret->op_m = CONSTANT;
    //                            ret->value_m = T(0.0);
    //    
    //    
    //                            return ret;
    //                        }
    //                    case ATAN:
    //                        if (this->left_m->HasId(id)) {
    //                            //f(x) = atan(x)
    //                            //f'(x) 1/(x^2+1)
    //    
    //                            ret->op_m = DIVIDE;
    //                            ret->left_m = new DynamicExpression<T > ();
    //                            ret->left_m->op_m = MULTIPLY;
    //                            ret->left_m->right_m = new DynamicExpression<T > ();
    //    
    //                            ret->left_m->right_m->op_m = CONSTANT;
    //                            ret->left_m->right_m->value_m = T(1.0);
    //    
    //    
    //                            ret->left_m->left_m = this->left_m->Differentiate(id);
    //    
    //                            ret->right_m = new DynamicExpression<T > ();
    //                            ret->right_m->op_m = PLUS;
    //    
    //                            ret->right_m->left_m = new DynamicExpression<T > ();
    //                            ret->right_m->left_m->op_m = MULTIPLY;
    //                            ret->right_m->left_m->left_m = this->left_m; //->Clone();
    //                            ret->right_m->left_m->right_m = this->left_m; //->Clone();
    //    
    //    
    //                            ret->right_m->right_m = new DynamicExpression<T > ();
    //                            ret->right_m->right_m->op_m = CONSTANT;
    //                            ret->right_m->right_m->value_m = T(1.0);
    //    
    //    
    //                            return ret;
    //                        } else {
    //                            ret->op_m = CONSTANT;
    //                            ret->value_m = T(0.0);
    //    
    //    
    //                            return ret;
    //    
    //                        }
    //                    case ATAN2:
    //                        //if w.r.t. check both expressions for id
    //                        if (this->left_m->HasId(id)) {
    //                            //f(x) = atan2(x,y)
    //                            //f'(x) y/(x^2+y^2)
    //    
    //                            ret->op_m = DIVIDE;
    //                            ret->left_m = new DynamicExpression<T > ();
    //                            ret->left_m->op_m = MULTIPLY;
    //                            ret->left_m->left_m = this->right_m; //->Clone(); //y
    //                            ret->left_m->right_m = left_m->Differentiate(id);
    //    
    //    
    //                            ret->right_m = new DynamicExpression<T > ();
    //                            ret->right_m->op_m = PLUS;
    //    
    //                            ret->right_m->left_m = new DynamicExpression<T > ();
    //                            ret->right_m->left_m->op_m = MULTIPLY;
    //                            ret->right_m->left_m->left_m = this->left_m; //->Clone();
    //                            ret->right_m->left_m->right_m = this->left_m; //->Clone();
    //    
    //    
    //                            ret->right_m->right_m = new DynamicExpression<T > ();
    //                            ret->right_m->right_m->op_m = MULTIPLY;
    //                            ret->right_m->right_m->left_m = this->right_m; //->Clone();
    //                            ret->right_m->right_m->right_m = this->right_m; //->Clone();
    //    
    //    
    //                            return ret;
    //                        } else {
    //                            ret->op_m = CONSTANT;
    //                            ret->value_m = T(0.0);
    //    
    //    
    //                            return ret;
    //                        }
    //                    case ATAN3:
    //    
    //                        //can be removed.
    //                        break;
    //    
    //                    case ATAN4:
    //                        break;
    //                    case SQRT:
    //                        if (this->left_m->HasId(id)) {
    //                            //f(x) = sqrt(x)
    //                            //f'(x) = .5/sqrt(x)
    //    
    //                            ret->op_m = DIVIDE;
    //                            ret->left_m = new DynamicExpression<T > ();
    //                            ret->left_m->op_m = MULTIPLY;
    //    
    //                            ret->left_m->right_m = new DynamicExpression<T > ();
    //                            ret->left_m->right_m->value_m = T(0.5);
    //    
    //                            ret->left_m->left_m = this->left_m->Differentiate(id);
    //    
    //    
    //                            ret->right_m = new DynamicExpression<T > ();
    //                            ret->right_m->op_m = SQRT;
    //                            ret->right_m->left_m = this->left_m; //->Clone();
    //    
    //                            //std::cout<<ret->ToString();
    //    
    //                            return ret;
    //                        } else {
    //                            ret->op_m = CONSTANT;
    //                            ret->value_m = T(0.0);
    //    
    //    
    //                            return ret;
    //                        }
    //                    case POW:
    //    
    //                        if (this->left_m->HasId(id)) {
    //                            //f(x) =  x^y
    //                            //f'(x) = yx^y-1
    //    
    //                            ret->op_m = MULTIPLY;
    //    
    //                            ret->left_m = new DynamicExpression<T > ();
    //                            ret->left_m->op_m = MULTIPLY;
    //                            ret->left_m->left_m = this->left_m->Differentiate(id);
    //                            ret->left_m->right_m = this->right_m; //->Clone();
    //    
    //    
    //    
    //                            ret->right_m = new DynamicExpression<T > ();
    //                            ret->right_m->op_m = POW;
    //    
    //    
    //                            ret->right_m->left_m = this->left_m; //->Clone();
    //    
    //    
    //                            ret->right_m->right_m = new DynamicExpression<T > ();
    //                            ret->right_m->right_m->op_m = MINUS;
    //                            ret->right_m->right_m->left_m = this->right_m; //->Clone();
    //    
    //                            ret->right_m->right_m->right_m = new DynamicExpression<T > ();
    //                            ret->right_m->right_m->right_m->op_m = CONSTANT;
    //                            ret->right_m->right_m->right_m->value_m = T(1.0);
    //    
    //    
    //    
    //                            return ret;
    //                        } else {
    //                            ret->op_m = CONSTANT;
    //                            ret->value_m = T(0.0);
    //    
    //    
    //                            return ret;
    //                        }
    //                        //                case POW1:
    //                        //
    //                        //                    break;
    //                        //                    //                return pow(this->left_value_, this->right_value_ - T(1.0));
    //                        //                case POW2:
    //                        //                    break;
    //                        //                    //                return pow(this->left_value_, this->right_value_ - T(1.0));
    //                    case LOG:
    //                        if (this->left_m->HasId(id)) {
    //                            //f(x) = log(x)
    //                            //f'(x) = 1/x
    //    
    //                            ret->op_m = DIVIDE;
    //                            ret->left_m = new DynamicExpression<T > ();
    //                            ret->left_m->op_m = MULTIPLY;
    //                            ret->left_m->left_m = new DynamicExpression<T > ();
    //                            ret->left_m->left_m->op_m = CONSTANT;
    //                            ret->left_m->left_m->value_m = T(1.0);
    //                            ret->left_m->right_m = this->left_m->Differentiate(id);
    //    
    //                            ret->right_m = this->left_m; //->Clone();
    //    
    //    
    //    
    //                            return ret;
    //                        } else {
    //                            ret->op_m = CONSTANT;
    //                            ret->value_m = T(0.0);
    //    
    //    
    //                            return ret;
    //                        }
    //                    case LOG10:
    //                        //f(x) = log10(x)
    //                        //f'(x) = 1/(xlog(10))
    //    
    //                        if (this->left_m->HasId(id)) {
    //    
    //    
    //    
    //                            ret->op_m = DIVIDE;
    //    
    //                            ret->left_m = new DynamicExpression<T > ();
    //                            ret->left_m->op_m = MULTIPLY;
    //    
    //                            ret->left_m->left_m = new DynamicExpression<T > ();
    //                            ret->left_m->left_m->op_m = CONSTANT;
    //                            ret->left_m->left_m->value_m = T(1.0);
    //    
    //                            ret->left_m->right_m = this->left_m->Differentiate(id);
    //    
    //                            ret->right_m = new DynamicExpression<T > ();
    //                            ret->right_m->op_m = MULTIPLY;
    //    
    //                            ret->right_m->left_m = this->left_m; //->Clone();
    //    
    //                            ret->right_m->right_m = new DynamicExpression<T > ();
    //                            ret->right_m->right_m->op_m = CONSTANT;
    //                            ret->right_m->right_m->value_m = log(T(10.0));
    //    
    //    
    //                            return ret;
    //                        } else {
    //                            ret->op_m = LOG;
    //                            ret->left_m = this; //->Clone();
    //    
    //    
    //                            return ret;
    //                        }
    //                    case EXP:
    //                        //f(x) = e^x
    //                        //f'(x) =e^x
    //    
    //                        if (this->left_m->HasId(id)) {
    //    
    //                            ret->op_m = MULTIPLY;
    //                            ret->left_m = this->left_m->Differentiate(id);
    //    
    //    
    //                            ret->right_m = new DynamicExpression<T > ();
    //                            ret->right_m->op_m = EXP;
    //                            ret->right_m->left_m = this->left_m; //->Clone();
    //    
    //    
    //    
    //                            return ret;
    //                        } else {
    //                            ret->op_m = CONSTANT;
    //                            ret->value_m = T(0.0);
    //    
    //    
    //                            return ret;
    //                        }
    //                    case SINH:
    //                        if (this->left_m->HasId(id)) {
    //                            //f(x) = sinh(x)
    //                            //f'(x) = cosh(x)
    //    
    //                            ret->op_m = MULTIPLY;
    //                            ret->left_m = this->left_m->Differentiate(id);
    //    
    //                            ret->right_m = new DynamicExpression<T > ();
    //                            ret->right_m->op_m = COSH;
    //                            ret->right_m->left_m = this->left_m; //->Clone();
    //    
    //    
    //                            return ret;
    //                        } else {
    //                            ret->op_m = CONSTANT;
    //                            ret->value_m = T(0.0);
    //    
    //    
    //                            return ret;
    //                        }
    //                    case COSH:
    //                        if (this->left_m->HasId(id)) {
    //    
    //                            ret->op_m = MULTIPLY;
    //                            ret->left_m = this->left_m->Differentiate(id);
    //    
    //                            ret->right_m = new DynamicExpression<T > ();
    //                            ret->right_m->op_m = SINH;
    //                            ret->right_m->left_m = this->left_m; //->Clone();
    //    
    //    
    //                            return ret;
    //                        } else {
    //                            ret->op_m = CONSTANT;
    //                            ret->value_m = T(0.0);
    //    
    //    
    //                            return ret;
    //                        }
    //                    case TANH:
    //                        //f(x) = tanh(x)
    //                        //f'(x) =1- tanh(x)*tanh(x)
    //    
    //    
    //                        if (this->left_m->HasId(id)) {
    //    
    //                            ret->op_m = MULTIPLY;
    //    
    //                            ret->left_m = this->left_m->Differentiate(id);
    //    
    //                            ret->right_m = new DynamicExpression<T > ();
    //                            ret->right_m->op_m = MULTIPLY;
    //                            ret->right_m->left_m = new DynamicExpression<T > ();
    //    
    //    
    //                            ret->right_m->left_m->op_m = DIVIDE;
    //                            ret->right_m->left_m->left_m = new DynamicExpression<T > ();
    //                            ret->right_m->left_m->left_m->op_m = CONSTANT;
    //                            ret->right_m->left_m->left_m->value_m = T(1.0);
    //    
    //    
    //                            ret->right_m->left_m->right_m = new DynamicExpression<T > ();
    //                            ret->right_m->left_m->right_m->op_m = COSH;
    //                            ret->right_m->left_m->right_m->left_m = this->left_m; //->Clone();
    //    
    //    
    //                            ret->right_m->right_m = ret->right_m->left_m; //->Clone();
    //    
    //    
    //                            return ret;
    //                        } else {
    //                            ret->op_m = CONSTANT;
    //                            ret->value_m = T(0.0);
    //    
    //    
    //                            return ret;
    //                        }
    //    
    //                    case FABS:
    //    
    //                        if (this->left_m->HasId(id)) {
    //    
    //                            ret->op_m = DIVIDE;
    //                            ret->left_m = new DynamicExpression<T > ();
    //                            ret->left_m->op_m = MULTIPLY;
    //    
    //                            ret->left_m->left_m = this->left_m->Differentiate(id);
    //                            ret->left_m->right_m = this->left_m; //->Clone();;
    //    
    //    
    //                            ret->right_m = new DynamicExpression<T > ();
    //                            ret->right_m->op_m = FABS;
    //                            ret->right_m->left_m = this->left_m; //->Clone();
    //    
    //    
    //                            return ret;
    //                        } else {
    //                            ret->op_m = CONSTANT;
    //                            ret->value_m = T(0.0);
    //    
    //    
    //                            return ret;
    //                        }
    //                    case FLOOR:
    //                        if (this->left_m->id_m == id) {
    //    
    //    
    //    
    //                            ret->op_m = MULTIPLY;
    //    
    //                            ret->left_m = this->left_m->Differentiate(id);
    //    
    //                            ret->right_m = new DynamicExpression<T > ();
    //                            ret->right_m->op_m = FLOOR;
    //                            ret->right_m->left_m = this->left_m; //->Clone();
    //    
    //    
    //                            return ret;
    //                        } else {
    //                            ret->op_m = CONSTANT;
    //                            ret->value_m = T(0.0);
    //    
    //    
    //                            return ret;
    //                        }
    //                    case NONE://shouldn't happen.
    //                        return this; //->Clone();
    //    
    //                    default:
    //                        return NULL;
    //                }
    //                return NULL;
    //            }
    //    
    //            /*!
    //             * Evaluate this expression. 
    //             * @return 
    //             */
    //            const T Evaluate() const {
    //    
    //                T l = T(0);
    //                T r = T(0);
    //    
    //                if (this->right_m != NULL) {
    //                    r = this->right_m->Evaluate();
    //                }
    //    
    //                if (this->left_m != NULL) {
    //                    l = this->left_m->Evaluate();
    //                }
    //    
    //    
    //    
    //    
    //                switch (op_m) {
    //                    case CONSTANT:
    //    
    //                        return this->value_m;
    //                    case VARIABLE:
    //    
    //                        return this->value_m;
    //                    case MINUS:
    //    
    //                        return (l - r);
    //                    case PLUS:
    //    
    //                        return (l + r);
    //                    case DIVIDE:
    //    
    //                        return (l / r);
    //                    case MULTIPLY:
    //    
    //                        return (l * r);
    //                    case SIN:
    //    
    //                        return sin(l);
    //                    case COS:
    //    
    //                        return cos(l);
    //                    case TAN:
    //    
    //                        return tan(l);
    //                    case ASIN:
    //    
    //                        return asin(l);
    //                    case ACOS:
    //    
    //                        return acos(l);
    //                    case ATAN:
    //    
    //                        return atan(l);
    //                    case ATAN2:
    //    
    //                        return atan2(l, r);
    //                        //                case ATAN3:
    //                        //                    break;
    //                        //                case ATAN4:
    //                        //                    break;
    //                    case SQRT:
    //    
    //                        return sqrt(l);
    //                    case POW:
    //    
    //                        return pow(l, r);
    //                        //                case POW1:
    //                        //                    break;
    //                        //                case POW2:
    //                        //                    break;
    //                    case LOG:
    //    
    //                        return log(l);
    //                    case LOG10:
    //    
    //                        return log10(l);
    //                    case EXP:
    //    
    //                        return exp(l);
    //                    case SINH:
    //    
    //                        return sinh(l);
    //                    case COSH:
    //    
    //                        return cosh(l);
    //                    case TANH:
    //    
    //                        return tanh(l);
    //                    case FABS:
    //    
    //                        return fabs(l);
    //                    case ABS:
    //    
    //                        return abs(l);
    //                    case FLOOR:
    //    
    //                        return floor(l);
    //                    case NONE:
    //    
    //                        return this->value_m;
    //                    default:
    //                        return T(0);
    //                }
    //                return T(0);
    //            }
    //    
    //            /**
    //             * Returns the evaluated derivative of this expression tree. While the
    //             * derivative is computed, no expression tree manipulations are made.
    //             * @param id
    //             * @return 
    //             */
    //            T EvaluateDerivative(const uint32_t &id, bool &has_id) {
    //                //#warning need to check partial derivatives....
    //    
    //                T ret, g, h = T(-999.0);
    //    
    //                T left_derivative = T(0);
    //                T right_derivative = T(0);
    //    
    //                if (this->GetLeft()) {
    //                    left_derivative = this->GetLeft()->EvaluateDerivative(id, has_id);
    //                }
    //    
    //                if (this->GetRight()) {
    //                    right_derivative = this->GetRight()->EvaluateDerivative(id, has_id);
    //                }
    //    
    //                //            has_id = this->HasId(id);
    //    
    //    
    //                switch (op_m) {
    //    
    //                    case CONSTANT:
    //                        //f(x) = C
    //                        //f'(x) = 0
    //                        return T(0);
    //    
    //                    case VARIABLE:
    //                        if (this->id_m == id) {
    //                            //f(x) = x
    //                            //f'(x) = 1
    //                            has_id = true;
    //    
    //                            return T(1.0);
    //                        } else {//constant
    //                            //f(x) = C
    //                            //f'(x) = 0
    //    
    //                            return T(0.0);
    //                        }
    //                    case MINUS:
    //    
    //                        //f(x) = g(x) - h(x)
    //                        //f'(x) = g'(x) - h'(x)
    //    
    //    
    //                        return left_derivative - right_derivative;
    //    
    //                    case PLUS:
    //    
    //                        //f(x) = g(x) + h(x)
    //                        //f'(x) = g'(x) + h'(x)
    //    
    //    
    //                        return left_derivative + right_derivative;
    //    
    //    
    //                    case DIVIDE:
    //    
    //                        //f(x) = g(x)/h(x);
    //                        //f'(x) = (g'(x)h(x) - g(x)h'(x))/h(x)^2
    //    
    //    
    //                        ret = (left_derivative * this->right_m->Evaluate() -
    //                                this->left_m->Evaluate() * right_derivative) /
    //                                (this->right_m->Evaluate() * this->right_m->Evaluate());
    //    
    //    
    //                        return ret;
    //    
    //                    case MULTIPLY:
    //                        //f(x) = g(x)h(x);
    //                        //f'(x) = g'(x)h(x) + g(x)h'(x)
    //    
    //                        if (this->left_m->op_m == CONSTANT
    //                                && this->right_m->op_m != CONSTANT) {
    //    
    //                            ret = this->left_m->Evaluate() * right_derivative;
    //    
    //                        } else if (this->right_m->op_m == CONSTANT
    //                                && this->left_m->op_m != CONSTANT) {
    //    
    //                            ret = left_derivative * this->right_m->Evaluate();
    //                            //                        std::cout<<"ret = "<<ret<<"\n";
    //                        } else {
    //    
    //                            //g'(x)h(x) + g(x)h'(x)
    //    
    //                            ret = left_derivative * this->right_m->Evaluate() +
    //                                    this->left_m->Evaluate() * right_derivative;
    //    
    //    
    //                        }
    //                        return ret;
    //    
    //                    case SIN:
    //    
    //                        if (has_id) {
    //                            //f'(x) = cos(x)
    //                            ret = left_derivative *
    //                                    std::cos(this->left_m->Evaluate());
    //    
    //                            return ret;
    //                        } else {
    //                            return T(0.0);
    //                        }
    //    
    //                    case COS:
    //    
    //                        if (has_id) {
    //                            //f'(x) = -sin(x)
    //    
    //    
    //                            g = left_derivative;
    //    
    //                            ret = g * T(-1.0) * std::sin(this->left_m->Evaluate());
    //    
    //                            return ret;
    //    
    //                        } else {
    //    
    //                            return T(0.0);
    //                        }
    //                    case TAN:
    //                        if (has_id) {
    //                            //f'(x) = 1/cos(x)
    //    
    //    
    //                            g = left_derivative;
    //    
    //                            ret = g * ((T(1.0) / std::cos(this->left_m->Evaluate()))*(T(1.0) / std::cos(this->left_m->Evaluate())));
    //    
    //    
    //                            return ret;
    //                        } else {
    //    
    //                            return T(0.0);
    //                        }
    //                    case ASIN:
    //    
    //                        if (has_id) {
    //    
    //    
    //                            //f(x) = asin(x)
    //                            //f'(x) = 1/(2 sqrt(1-x^2)= 1/(pow((1-pow(x,2)),0.5)
    //    
    //    
    //                            g = left_derivative;
    //    
    //                            ret = (g * T(1.0) / std::pow((T(1.0) - std::pow(this->left_m->Evaluate(), T(2.0))), T(0.5)));
    //    
    //                            return ret;
    //                        } else {
    //                            return T(0.0);
    //                        }
    //                    case ACOS:
    //    
    //                        if (has_id) {
    //                            g = left_derivative;
    //    
    //                            ret = (g * T(-1.0) / std::pow((T(1.0) - std::pow(this->left_m->Evaluate(), T(2.0))), T(0.5)));
    //    
    //                            return ret;
    //                        } else {
    //    
    //                            return T(0.0);
    //                        }
    //                    case ATAN:
    //                        if (has_id) {
    //                            g = left_derivative;
    //                            ret = (g * T(1.0) / (this->left_m->Evaluate() * this->left_m->Evaluate() + T(1.0)));
    //    
    //                            return ret;
    //                        } else {
    //                            //                        ret->op_m = CONSTANT;
    //                            //                        ret->value_m = T(0.0);
    //                            return T(0.0);
    //    
    //                        }
    //                    case ATAN2:
    //                        //if w.r.t. check both expressions for id
    //                        if (has_id) {
    //                            //f(x) = atan2(x,y)
    //                            //f'(x) y/(x^2+y^2)
    //    
    //                            g = left_derivative;
    //                            ret = (this->right_m->Evaluate() * g / (this->left_m->Evaluate() * this->left_m->Evaluate()+(this->right_m->Evaluate() * this->right_m->Evaluate())));
    //    
    //                            return ret;
    //                        } else {
    //    
    //                            return T(0.0);
    //                        }
    //                    case ATAN3:
    //    
    //                        //can be removed.
    //                        break;
    //    
    //                    case ATAN4:
    //                        break;
    //                    case SQRT:
    //                        if (has_id) {
    //                            //f(x) = sqrt(x)
    //                            //f'(x) = .5/sqrt(x)
    //                            g = left_derivative;
    //                            ret = g * T(.5) / std::sqrt(this->left_m->Evaluate());
    //    
    //    
    //                            return ret;
    //                        } else {
    //                            return T(0.0);
    //                        }
    //                    case POW:
    //    
    //                        if (has_id) {
    //                            //f(x) =  x^y
    //                            //f'(x) = yx^y-1
    //                            ret = (left_derivative * this->right_m->Evaluate()) *
    //                                    std::pow(this->left_m->Evaluate(), (this->right_m->Evaluate() - T(1.0)));
    //    
    //                            return ret;
    //                        } else {
    //    
    //                            return T(0.0);
    //                        }
    //    
    //                    case LOG:
    //                        if (has_id) {
    //                            //f(x) = log(x)
    //                            //f'(x) = 1/x
    //                            ret = (left_derivative * T(1.0)) / this->left_m->Evaluate();
    //    
    //                            return ret;
    //                        } else {
    //    
    //                            return T(0.0);
    //                        }
    //                    case LOG10:
    //                        //f(x) = log10(x)
    //                        //f'(x) = 1/(xlog(10))
    //    
    //                        if (has_id) {
    //    
    //                            ret = (left_derivative * T(1.0)) / (this->left_m->Evaluate() * std::log(T(10.0)));
    //    
    //                            return ret;
    //                        } else {
    //                            return T(0.0);
    //                        }
    //                    case EXP:
    //                        //f(x) = e^x
    //                        //f'(x) =e^x
    //    
    //                        if (has_id) {
    //                            ret = left_derivative * std::exp(this->left_m->Evaluate());
    //    
    //                            return ret;
    //                        } else {
    //    
    //                            return T(0.0);
    //                        }
    //                    case SINH:
    //                        if (has_id) {
    //                            //f(x) = sinh(x)
    //                            //f'(x) = cosh(x)
    //                            return left_derivative * std::cosh(this->left_m->Evaluate());
    //    
    //                            return ret;
    //                        } else {
    //    
    //                            return T(0.0);
    //                        }
    //                    case COSH:
    //                        if (has_id) {
    //                            return left_derivative * std::sinh(this->left_m->Evaluate());
    //    
    //                            return ret;
    //                        } else {
    //    
    //                            return ret;
    //                        }
    //                    case TANH:
    //                        //f(x) = tanh(x)
    //                        //f'(x) =1- tanh(x)*tanh(x)
    //    
    //    
    //                        if (has_id) {
    //    
    //                            ret = left_derivative * (T(1.0) / std::cosh(this->left_m->Evaluate()))*(T(1.0) / std::cosh(this->left_m->Evaluate()));
    //    
    //    
    //                            return ret;
    //                        } else {
    //    
    //                            return T(0.0);
    //                        }
    //    
    //                    case FABS:
    //    
    //                        if (has_id) {
    //    
    //                            ret = (left_derivative * this->left_m->Evaluate()) /
    //                                    std::fabs(this->left_m->Evaluate());
    //    
    //                            return ret;
    //                        } else {
    //    
    //                            return T(0.0);
    //                        }
    //                    case FLOOR:
    //                        if (has_id) {
    //    
    //                            ret = left_derivative * std::floor(this->left_m->Evaluate());
    //    
    //                            return ret;
    //                        } else {
    //    
    //                            return ret;
    //                        }
    //                    case NONE://shouldn't happen.
    //                        return ret;
    //    
    //                    default:
    //                        return ret;
    //                }
    //    
    //                return 0.0;
    //            }
    //    
    //            const unsigned long GetId() const {
    //                return id_m;
    //            }
    //    
    //            void SetId(const unsigned long &id) {
    //                id_m = id;
    //            }
    //    
    //            DynamicExpression<T>* GetLeft() const {
    //                return left_m;
    //            }
    //    
    //            void SetLeft(DynamicExpression<T>* left) {
    //                left_m = left;
    //            }
    //    
    //            const Operation GetOp() const {
    //                return op_m;
    //            }
    //    
    //            void SetOp(const Operation &op) {
    //                op_m = op;
    //            }
    //    
    //            DynamicExpression<T>* GetRight() const {
    //                return right_m;
    //            }
    //    
    //            void SetRight(DynamicExpression<T>* right) {
    //                right_m = right;
    //            }
    //    
    //            const T GetValue() const {
    //                return value_m;
    //            }
    //    
    //            void SetValue(T value) {
    //                value_m = value;
    //            }
    //    
    //            const std::string GetName() const {
    //    
    //                if (name_m == "") {
    //                    std::stringstream ss;
    //                    ss << "x" << GetId();
    //                    return ss.str();
    //                }
    //                return name_m;
    //            }
    //    
    //            void SetName(const std::string &name) {
    //                name_m = name;
    //            }
    //    
    //            /*!
    //             * Represent this expression as a string. ADNumbers are represented
    //             * in wkt format by default unless latex flag is set to true.
    //             * Constants are represented by value.
    //             *
    //             */
    //            std::string ToString(bool latex = false) {
    //    
    //                //            if (Size() > 1000) {
    //                //                return std::string("Error: ToString...expression too long.");
    //                //            }
    //                //
    //                std::stringstream ss;
    //                std::stringstream temp;
    //    
    //                std::string l, r;
    //    
    //    
    //    
    //                if (GetRight() != NULL) {
    //                    r = GetRight()->ToString(latex);
    //                }
    //    
    //                if (GetLeft() != NULL) {
    //                    l = GetLeft()->ToString(latex);
    //                }
    //    
    //    
    //                temp.str("");
    //    
    //    
    //    
    //                switch (GetOp()) {
    //                    case CONSTANT:
    //                        if (latex) {
    //                            ss << GetValue();
    //                        } else {
    //                            ss << "CONST[" << GetValue() << "]";
    //                        }
    //                        break;
    //                    case VARIABLE:
    //                        if (latex) {
    //                            if (GetName() == "") {
    //                                ss << "x" << GetId();
    //                            } else {
    //                                ss << GetName();
    //                            }
    //                        } else {
    //                            ss << "VAR[" << GetValue() << ",ID[" << GetId() << "]" << "]";
    //                        }
    //                        break;
    //                    case MINUS:
    //                        ss << "(" << l << " - " << r << ")";
    //                        break;
    //                    case PLUS:
    //                        ss << "(" << l << " + " << r << ")";
    //                        break;
    //                    case DIVIDE:
    //                        if (latex) {
    //                            ss << "\\frac{" << l << "}{" << r << "}";
    //                        } else {
    //                            ss << l << " / " << r;
    //                        }
    //                        break;
    //                    case MULTIPLY:
    //    
    //                        if (latex && GetRight()->GetOp() == POW
    //                                && GetRight()->GetRight()->GetOp() == CONSTANT
    //                                && GetRight()->GetRight()->GetValue() == T(-1)) {
    //                            ss << "\\frac{" << l << "}{" << GetRight()->GetLeft()->ToString(latex) << "}";
    //                        } else {
    //    
    //                            ss << "(" << l << " * " << r << ")";
    //                        }
    //                        break;
    //                    case SIN:
    //                        ss << "sin(" << l << ")";
    //                        break;
    //                    case COS:
    //                        ss << "cos(" << l << ")";
    //                        break;
    //                    case TAN:
    //                        ss << "tan(" << l << ")";
    //                        break;
    //                    case ASIN:
    //                        ss << "asin(" << l << ")";
    //                        break;
    //                    case ACOS:
    //                        ss << "acos(" << l << ")";
    //                        break;
    //                    case ATAN:
    //                        ss << "atan(" << l << ")";
    //                        break;
    //                    case ATAN2:
    //                        ss << "atan(" << l << "," << r << ")";
    //                        break;
    //                    case ATAN3:
    //                        ss << "atan(" << l << "," << r << ")";
    //                        break;
    //                    case ATAN4:
    //                        ss << "atan(" << l << "," << r << ")";
    //                        break;
    //                    case SQRT:
    //                        ss << "sqrt(" << l << ")";
    //                        break;
    //                    case POW:
    //    
    //    
    //                    case POW1:
    //    
    //                    case POW2:
    //                        if (latex) {
    //                            ss << l << "^{" << r << "}";
    //                        } else {
    //                            ss << "pow(" << l << "," << r << ")";
    //                        }
    //                        break;
    //                    case LOG:
    //                        if (latex) {
    //                            ss << "ln(" << l << ")";
    //                        } else {
    //                            ss << "log(" << l << ")";
    //                        }
    //                        break;
    //                    case LOG10:
    //                        ss << "log10(" << l << ")";
    //                        break;
    //                    case EXP:
    //                        ss << "exp(" << l << ")";
    //                        break;
    //                    case SINH:
    //                        ss << "sinh(" << l << ")";
    //                        break;
    //                    case COSH:
    //                        ss << "cosh(" << l << ")";
    //                        break;
    //                    case TANH:
    //                        ss << "tanh(" << l << ")";
    //                        break;
    //                    case FABS:
    //                        ss << "fabs(" << l << ")";
    //                    case FLOOR:
    //                        ss << "floor(" << l << ")";
    //                    case NONE:
    //                        break;
    //                    default:
    //                        break;
    //                }
    //    
    //                //ss << ")";
    //    
    //                return ss.str();
    //    
    //            }
    //    
    //    
    //        };
    //    
    //        template<class T>
    //        std::vector<DynamicExpression<T> > DynamicExpression<T>::pool(1);
    //    
    //    
    //        template<class T>
    //        std::atomic<uint32_t> DynamicExpression<T>::counter;

}





#endif	/* DYNAMICEXPRESSION_HPP */

