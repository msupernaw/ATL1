#ifndef MATRIXMATH_HPP
#define MATRIXMATH_HPP

#ifdef ENABLE_BOUNDS_CHECKING
#include <assert.h>
#endif




namespace atl{

 template<class C>
    struct MatrixACos : atl::MatrixExpression<typename C::RET_TYPE, MatrixACos<C> > {
        const C& c_m;
        typedef typename C::RET_TYPE RET_TYPE;
        typedef typename C::BASE_TYPE BASE_TYPE;

        inline explicit MatrixACos(const atl::MatrixExpression<typename C::RET_TYPE,C>& c) : c_m(c.Cast()) {
        inline explicit MatrixACos(const atl::MatrixExpression<typename C::RET_TYPE,C>& c) : c_m(c.Cast()) {

        }

        inline const size_t Size(const int32_t & dimension) const {
            switch (dimension) {
                case 0:
                    return c_m.Size(1);
                case 1:
                    return c_m.Size(0);
                default:
                    return 0;
            }
        }

        inline const RET_TYPE operator()(size_t i, size_t j) const {
#ifdef ENABLE_BOUNDS_CHECKING
            assert(i < c_m.Size(0));
            assert(j < c_m.Size(1));
#endif   

            return ::acos(c_m(i,j));

        }
        
          
        /*
         *
         *Returns the first valid index for this matrix given dimension d. 
         */
        inline const size_t IndexMin(uint32_t d) const {
            return c_m.IndexMin(d);
        }

        /*
         *
         *Returns the last valid index for this matrix given dimension d. 
         */
        inline const size_t IndexMax(uint32_t d) const {
            return c_m.IndexMax(d);
        }

        /**
         * Get a value based on the raw indices for the underlying 
         * data. 
         * @param i
         * @return 
         */
        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j) const {
            return ::acos(c_m.AtRaw(i,j));
        }

        /**
         * Used to determine if the resulting matrix is used in the expression.
         * 
         * @param aliased
         * @param ptr
         */
        inline void IsAliased(bool& aliased, void* ptr) const {
            c_m.IsAliased(aliased, ptr);
        }

        inline void ExpressionLength(uint32_t& length)const {
            length++;
            c_m.ExpressionLength(length);
        }


    };

    template<class C>
    inline const MatrixACos<C> acos(const atl::MatrixExpression<typename C::RET_TYPE,C>& c) {
        return MatrixACos<C>(c);
    }
 template<class C>
    struct MatrixASin : atl::MatrixExpression<typename C::RET_TYPE, MatrixASin<C> > {
        const C& c_m;
        typedef typename C::RET_TYPE RET_TYPE;
        typedef typename C::BASE_TYPE BASE_TYPE;

        inline explicit MatrixASin(const atl::MatrixExpression<typename C::RET_TYPE,C>& c) : c_m(c.Cast()) {
        inline explicit MatrixASin(const atl::MatrixExpression<typename C::RET_TYPE,C>& c) : c_m(c.Cast()) {

        }

        inline const size_t Size(const int32_t & dimension) const {
            switch (dimension) {
                case 0:
                    return c_m.Size(1);
                case 1:
                    return c_m.Size(0);
                default:
                    return 0;
            }
        }

        inline const RET_TYPE operator()(size_t i, size_t j) const {
#ifdef ENABLE_BOUNDS_CHECKING
            assert(i < c_m.Size(0));
            assert(j < c_m.Size(1));
#endif   

            return ::asin(c_m(i,j));

        }
        
          
        /*
         *
         *Returns the first valid index for this matrix given dimension d. 
         */
        inline const size_t IndexMin(uint32_t d) const {
            return c_m.IndexMin(d);
        }

        /*
         *
         *Returns the last valid index for this matrix given dimension d. 
         */
        inline const size_t IndexMax(uint32_t d) const {
            return c_m.IndexMax(d);
        }

        /**
         * Get a value based on the raw indices for the underlying 
         * data. 
         * @param i
         * @return 
         */
        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j) const {
            return ::asin(c_m.AtRaw(i,j));
        }

        /**
         * Used to determine if the resulting matrix is used in the expression.
         * 
         * @param aliased
         * @param ptr
         */
        inline void IsAliased(bool& aliased, void* ptr) const {
            c_m.IsAliased(aliased, ptr);
        }

        inline void ExpressionLength(uint32_t& length)const {
            length++;
            c_m.ExpressionLength(length);
        }


    };

    template<class C>
    inline const MatrixASin<C> asin(const atl::MatrixExpression<typename C::RET_TYPE,C>& c) {
        return MatrixASin<C>(c);
    }
 template<class C>
    struct MatrixATan : atl::MatrixExpression<typename C::RET_TYPE, MatrixATan<C> > {
        const C& c_m;
        typedef typename C::RET_TYPE RET_TYPE;
        typedef typename C::BASE_TYPE BASE_TYPE;

        inline explicit MatrixATan(const atl::MatrixExpression<typename C::RET_TYPE,C>& c) : c_m(c.Cast()) {
        inline explicit MatrixATan(const atl::MatrixExpression<typename C::RET_TYPE,C>& c) : c_m(c.Cast()) {

        }

        inline const size_t Size(const int32_t & dimension) const {
            switch (dimension) {
                case 0:
                    return c_m.Size(1);
                case 1:
                    return c_m.Size(0);
                default:
                    return 0;
            }
        }

        inline const RET_TYPE operator()(size_t i, size_t j) const {
#ifdef ENABLE_BOUNDS_CHECKING
            assert(i < c_m.Size(0));
            assert(j < c_m.Size(1));
#endif   

            return ::atan(c_m(i,j));

        }
        
          
        /*
         *
         *Returns the first valid index for this matrix given dimension d. 
         */
        inline const size_t IndexMin(uint32_t d) const {
            return c_m.IndexMin(d);
        }

        /*
         *
         *Returns the last valid index for this matrix given dimension d. 
         */
        inline const size_t IndexMax(uint32_t d) const {
            return c_m.IndexMax(d);
        }

        /**
         * Get a value based on the raw indices for the underlying 
         * data. 
         * @param i
         * @return 
         */
        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j) const {
            return ::atan(c_m.AtRaw(i,j));
        }

        /**
         * Used to determine if the resulting matrix is used in the expression.
         * 
         * @param aliased
         * @param ptr
         */
        inline void IsAliased(bool& aliased, void* ptr) const {
            c_m.IsAliased(aliased, ptr);
        }

        inline void ExpressionLength(uint32_t& length)const {
            length++;
            c_m.ExpressionLength(length);
        }


    };

    template<class C>
    inline const MatrixATan<C> atan(const atl::MatrixExpression<typename C::RET_TYPE,C>& c) {
        return MatrixATan<C>(c);
    }
 template<class C>
    struct MatrixACeil : atl::MatrixExpression<typename C::RET_TYPE, MatrixACeil<C> > {
        const C& c_m;
        typedef typename C::RET_TYPE RET_TYPE;
        typedef typename C::BASE_TYPE BASE_TYPE;

        inline explicit MatrixACeil(const atl::MatrixExpression<typename C::RET_TYPE,C>& c) : c_m(c.Cast()) {
        inline explicit MatrixACeil(const atl::MatrixExpression<typename C::RET_TYPE,C>& c) : c_m(c.Cast()) {

        }

        inline const size_t Size(const int32_t & dimension) const {
            switch (dimension) {
                case 0:
                    return c_m.Size(1);
                case 1:
                    return c_m.Size(0);
                default:
                    return 0;
            }
        }

        inline const RET_TYPE operator()(size_t i, size_t j) const {
#ifdef ENABLE_BOUNDS_CHECKING
            assert(i < c_m.Size(0));
            assert(j < c_m.Size(1));
#endif   

            return ::ceil(c_m(i,j));

        }
        
          
        /*
         *
         *Returns the first valid index for this matrix given dimension d. 
         */
        inline const size_t IndexMin(uint32_t d) const {
            return c_m.IndexMin(d);
        }

        /*
         *
         *Returns the last valid index for this matrix given dimension d. 
         */
        inline const size_t IndexMax(uint32_t d) const {
            return c_m.IndexMax(d);
        }

        /**
         * Get a value based on the raw indices for the underlying 
         * data. 
         * @param i
         * @return 
         */
        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j) const {
            return ::ceil(c_m.AtRaw(i,j));
        }

        /**
         * Used to determine if the resulting matrix is used in the expression.
         * 
         * @param aliased
         * @param ptr
         */
        inline void IsAliased(bool& aliased, void* ptr) const {
            c_m.IsAliased(aliased, ptr);
        }

        inline void ExpressionLength(uint32_t& length)const {
            length++;
            c_m.ExpressionLength(length);
        }


    };

    template<class C>
    inline const MatrixACeil<C> ceil(const atl::MatrixExpression<typename C::RET_TYPE,C>& c) {
        return MatrixACeil<C>(c);
    }
 template<class C>
    struct MatrixCos : atl::MatrixExpression<typename C::RET_TYPE, MatrixCos<C> > {
        const C& c_m;
        typedef typename C::RET_TYPE RET_TYPE;
        typedef typename C::BASE_TYPE BASE_TYPE;

        inline explicit MatrixCos(const atl::MatrixExpression<typename C::RET_TYPE,C>& c) : c_m(c.Cast()) {
        inline explicit MatrixCos(const atl::MatrixExpression<typename C::RET_TYPE,C>& c) : c_m(c.Cast()) {

        }

        inline const size_t Size(const int32_t & dimension) const {
            switch (dimension) {
                case 0:
                    return c_m.Size(1);
                case 1:
                    return c_m.Size(0);
                default:
                    return 0;
            }
        }

        inline const RET_TYPE operator()(size_t i, size_t j) const {
#ifdef ENABLE_BOUNDS_CHECKING
            assert(i < c_m.Size(0));
            assert(j < c_m.Size(1));
#endif   

            return ::cos(c_m(i,j));

        }
        
          
        /*
         *
         *Returns the first valid index for this matrix given dimension d. 
         */
        inline const size_t IndexMin(uint32_t d) const {
            return c_m.IndexMin(d);
        }

        /*
         *
         *Returns the last valid index for this matrix given dimension d. 
         */
        inline const size_t IndexMax(uint32_t d) const {
            return c_m.IndexMax(d);
        }

        /**
         * Get a value based on the raw indices for the underlying 
         * data. 
         * @param i
         * @return 
         */
        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j) const {
            return ::cos(c_m.AtRaw(i,j));
        }

        /**
         * Used to determine if the resulting matrix is used in the expression.
         * 
         * @param aliased
         * @param ptr
         */
        inline void IsAliased(bool& aliased, void* ptr) const {
            c_m.IsAliased(aliased, ptr);
        }

        inline void ExpressionLength(uint32_t& length)const {
            length++;
            c_m.ExpressionLength(length);
        }


    };

    template<class C>
    inline const MatrixCos<C> cos(const atl::MatrixExpression<typename C::RET_TYPE,C>& c) {
        return MatrixCos<C>(c);
    }
 template<class C>
    struct MatrixCosh : atl::MatrixExpression<typename C::RET_TYPE, MatrixCosh<C> > {
        const C& c_m;
        typedef typename C::RET_TYPE RET_TYPE;
        typedef typename C::BASE_TYPE BASE_TYPE;

        inline explicit MatrixCosh(const atl::MatrixExpression<typename C::RET_TYPE,C>& c) : c_m(c.Cast()) {
        inline explicit MatrixCosh(const atl::MatrixExpression<typename C::RET_TYPE,C>& c) : c_m(c.Cast()) {

        }

        inline const size_t Size(const int32_t & dimension) const {
            switch (dimension) {
                case 0:
                    return c_m.Size(1);
                case 1:
                    return c_m.Size(0);
                default:
                    return 0;
            }
        }

        inline const RET_TYPE operator()(size_t i, size_t j) const {
#ifdef ENABLE_BOUNDS_CHECKING
            assert(i < c_m.Size(0));
            assert(j < c_m.Size(1));
#endif   

            return ::cosh(c_m(i,j));

        }
        
          
        /*
         *
         *Returns the first valid index for this matrix given dimension d. 
         */
        inline const size_t IndexMin(uint32_t d) const {
            return c_m.IndexMin(d);
        }

        /*
         *
         *Returns the last valid index for this matrix given dimension d. 
         */
        inline const size_t IndexMax(uint32_t d) const {
            return c_m.IndexMax(d);
        }

        /**
         * Get a value based on the raw indices for the underlying 
         * data. 
         * @param i
         * @return 
         */
        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j) const {
            return ::cosh(c_m.AtRaw(i,j));
        }

        /**
         * Used to determine if the resulting matrix is used in the expression.
         * 
         * @param aliased
         * @param ptr
         */
        inline void IsAliased(bool& aliased, void* ptr) const {
            c_m.IsAliased(aliased, ptr);
        }

        inline void ExpressionLength(uint32_t& length)const {
            length++;
            c_m.ExpressionLength(length);
        }


    };

    template<class C>
    inline const MatrixCosh<C> cosh(const atl::MatrixExpression<typename C::RET_TYPE,C>& c) {
        return MatrixCosh<C>(c);
    }
 template<class C>
    struct MatrixExp : atl::MatrixExpression<typename C::RET_TYPE, MatrixExp<C> > {
        const C& c_m;
        typedef typename C::RET_TYPE RET_TYPE;
        typedef typename C::BASE_TYPE BASE_TYPE;

        inline explicit MatrixExp(const atl::MatrixExpression<typename C::RET_TYPE,C>& c) : c_m(c.Cast()) {
        inline explicit MatrixExp(const atl::MatrixExpression<typename C::RET_TYPE,C>& c) : c_m(c.Cast()) {

        }

        inline const size_t Size(const int32_t & dimension) const {
            switch (dimension) {
                case 0:
                    return c_m.Size(1);
                case 1:
                    return c_m.Size(0);
                default:
                    return 0;
            }
        }

        inline const RET_TYPE operator()(size_t i, size_t j) const {
#ifdef ENABLE_BOUNDS_CHECKING
            assert(i < c_m.Size(0));
            assert(j < c_m.Size(1));
#endif   

            return ::exp(c_m(i,j));

        }
        
          
        /*
         *
         *Returns the first valid index for this matrix given dimension d. 
         */
        inline const size_t IndexMin(uint32_t d) const {
            return c_m.IndexMin(d);
        }

        /*
         *
         *Returns the last valid index for this matrix given dimension d. 
         */
        inline const size_t IndexMax(uint32_t d) const {
            return c_m.IndexMax(d);
        }

        /**
         * Get a value based on the raw indices for the underlying 
         * data. 
         * @param i
         * @return 
         */
        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j) const {
            return ::exp(c_m.AtRaw(i,j));
        }

        /**
         * Used to determine if the resulting matrix is used in the expression.
         * 
         * @param aliased
         * @param ptr
         */
        inline void IsAliased(bool& aliased, void* ptr) const {
            c_m.IsAliased(aliased, ptr);
        }

        inline void ExpressionLength(uint32_t& length)const {
            length++;
            c_m.ExpressionLength(length);
        }


    };

    template<class C>
    inline const MatrixExp<C> exp(const atl::MatrixExpression<typename C::RET_TYPE,C>& c) {
        return MatrixExp<C>(c);
    }
 template<class C>
    struct MatrixFabs : atl::MatrixExpression<typename C::RET_TYPE, MatrixFabs<C> > {
        const C& c_m;
        typedef typename C::RET_TYPE RET_TYPE;
        typedef typename C::BASE_TYPE BASE_TYPE;

        inline explicit MatrixFabs(const atl::MatrixExpression<typename C::RET_TYPE,C>& c) : c_m(c.Cast()) {
        inline explicit MatrixFabs(const atl::MatrixExpression<typename C::RET_TYPE,C>& c) : c_m(c.Cast()) {

        }

        inline const size_t Size(const int32_t & dimension) const {
            switch (dimension) {
                case 0:
                    return c_m.Size(1);
                case 1:
                    return c_m.Size(0);
                default:
                    return 0;
            }
        }

        inline const RET_TYPE operator()(size_t i, size_t j) const {
#ifdef ENABLE_BOUNDS_CHECKING
            assert(i < c_m.Size(0));
            assert(j < c_m.Size(1));
#endif   

            return ::fabs(c_m(i,j));

        }
        
          
        /*
         *
         *Returns the first valid index for this matrix given dimension d. 
         */
        inline const size_t IndexMin(uint32_t d) const {
            return c_m.IndexMin(d);
        }

        /*
         *
         *Returns the last valid index for this matrix given dimension d. 
         */
        inline const size_t IndexMax(uint32_t d) const {
            return c_m.IndexMax(d);
        }

        /**
         * Get a value based on the raw indices for the underlying 
         * data. 
         * @param i
         * @return 
         */
        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j) const {
            return ::fabs(c_m.AtRaw(i,j));
        }

        /**
         * Used to determine if the resulting matrix is used in the expression.
         * 
         * @param aliased
         * @param ptr
         */
        inline void IsAliased(bool& aliased, void* ptr) const {
            c_m.IsAliased(aliased, ptr);
        }

        inline void ExpressionLength(uint32_t& length)const {
            length++;
            c_m.ExpressionLength(length);
        }


    };

    template<class C>
    inline const MatrixFabs<C> fabs(const atl::MatrixExpression<typename C::RET_TYPE,C>& c) {
        return MatrixFabs<C>(c);
    }
 template<class C>
    struct MatrixFloor : atl::MatrixExpression<typename C::RET_TYPE, MatrixFloor<C> > {
        const C& c_m;
        typedef typename C::RET_TYPE RET_TYPE;
        typedef typename C::BASE_TYPE BASE_TYPE;

        inline explicit MatrixFloor(const atl::MatrixExpression<typename C::RET_TYPE,C>& c) : c_m(c.Cast()) {
        inline explicit MatrixFloor(const atl::MatrixExpression<typename C::RET_TYPE,C>& c) : c_m(c.Cast()) {

        }

        inline const size_t Size(const int32_t & dimension) const {
            switch (dimension) {
                case 0:
                    return c_m.Size(1);
                case 1:
                    return c_m.Size(0);
                default:
                    return 0;
            }
        }

        inline const RET_TYPE operator()(size_t i, size_t j) const {
#ifdef ENABLE_BOUNDS_CHECKING
            assert(i < c_m.Size(0));
            assert(j < c_m.Size(1));
#endif   

            return ::floor(c_m(i,j));

        }
        
          
        /*
         *
         *Returns the first valid index for this matrix given dimension d. 
         */
        inline const size_t IndexMin(uint32_t d) const {
            return c_m.IndexMin(d);
        }

        /*
         *
         *Returns the last valid index for this matrix given dimension d. 
         */
        inline const size_t IndexMax(uint32_t d) const {
            return c_m.IndexMax(d);
        }

        /**
         * Get a value based on the raw indices for the underlying 
         * data. 
         * @param i
         * @return 
         */
        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j) const {
            return ::floor(c_m.AtRaw(i,j));
        }

        /**
         * Used to determine if the resulting matrix is used in the expression.
         * 
         * @param aliased
         * @param ptr
         */
        inline void IsAliased(bool& aliased, void* ptr) const {
            c_m.IsAliased(aliased, ptr);
        }

        inline void ExpressionLength(uint32_t& length)const {
            length++;
            c_m.ExpressionLength(length);
        }


    };

    template<class C>
    inline const MatrixFloor<C> floor(const atl::MatrixExpression<typename C::RET_TYPE,C>& c) {
        return MatrixFloor<C>(c);
    }
 template<class C>
    struct MatrixLog : atl::MatrixExpression<typename C::RET_TYPE, MatrixLog<C> > {
        const C& c_m;
        typedef typename C::RET_TYPE RET_TYPE;
        typedef typename C::BASE_TYPE BASE_TYPE;

        inline explicit MatrixLog(const atl::MatrixExpression<typename C::RET_TYPE,C>& c) : c_m(c.Cast()) {
        inline explicit MatrixLog(const atl::MatrixExpression<typename C::RET_TYPE,C>& c) : c_m(c.Cast()) {

        }

        inline const size_t Size(const int32_t & dimension) const {
            switch (dimension) {
                case 0:
                    return c_m.Size(1);
                case 1:
                    return c_m.Size(0);
                default:
                    return 0;
            }
        }

        inline const RET_TYPE operator()(size_t i, size_t j) const {
#ifdef ENABLE_BOUNDS_CHECKING
            assert(i < c_m.Size(0));
            assert(j < c_m.Size(1));
#endif   

            return ::log(c_m(i,j));

        }
        
          
        /*
         *
         *Returns the first valid index for this matrix given dimension d. 
         */
        inline const size_t IndexMin(uint32_t d) const {
            return c_m.IndexMin(d);
        }

        /*
         *
         *Returns the last valid index for this matrix given dimension d. 
         */
        inline const size_t IndexMax(uint32_t d) const {
            return c_m.IndexMax(d);
        }

        /**
         * Get a value based on the raw indices for the underlying 
         * data. 
         * @param i
         * @return 
         */
        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j) const {
            return ::log(c_m.AtRaw(i,j));
        }

        /**
         * Used to determine if the resulting matrix is used in the expression.
         * 
         * @param aliased
         * @param ptr
         */
        inline void IsAliased(bool& aliased, void* ptr) const {
            c_m.IsAliased(aliased, ptr);
        }

        inline void ExpressionLength(uint32_t& length)const {
            length++;
            c_m.ExpressionLength(length);
        }


    };

    template<class C>
    inline const MatrixLog<C> log(const atl::MatrixExpression<typename C::RET_TYPE,C>& c) {
        return MatrixLog<C>(c);
    }
 template<class C>
    struct MatrixLog10 : atl::MatrixExpression<typename C::RET_TYPE, MatrixLog10<C> > {
        const C& c_m;
        typedef typename C::RET_TYPE RET_TYPE;
        typedef typename C::BASE_TYPE BASE_TYPE;

        inline explicit MatrixLog10(const atl::MatrixExpression<typename C::RET_TYPE,C>& c) : c_m(c.Cast()) {
        inline explicit MatrixLog10(const atl::MatrixExpression<typename C::RET_TYPE,C>& c) : c_m(c.Cast()) {

        }

        inline const size_t Size(const int32_t & dimension) const {
            switch (dimension) {
                case 0:
                    return c_m.Size(1);
                case 1:
                    return c_m.Size(0);
                default:
                    return 0;
            }
        }

        inline const RET_TYPE operator()(size_t i, size_t j) const {
#ifdef ENABLE_BOUNDS_CHECKING
            assert(i < c_m.Size(0));
            assert(j < c_m.Size(1));
#endif   

            return ::log10(c_m(i,j));

        }
        
          
        /*
         *
         *Returns the first valid index for this matrix given dimension d. 
         */
        inline const size_t IndexMin(uint32_t d) const {
            return c_m.IndexMin(d);
        }

        /*
         *
         *Returns the last valid index for this matrix given dimension d. 
         */
        inline const size_t IndexMax(uint32_t d) const {
            return c_m.IndexMax(d);
        }

        /**
         * Get a value based on the raw indices for the underlying 
         * data. 
         * @param i
         * @return 
         */
        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j) const {
            return ::log10(c_m.AtRaw(i,j));
        }

        /**
         * Used to determine if the resulting matrix is used in the expression.
         * 
         * @param aliased
         * @param ptr
         */
        inline void IsAliased(bool& aliased, void* ptr) const {
            c_m.IsAliased(aliased, ptr);
        }

        inline void ExpressionLength(uint32_t& length)const {
            length++;
            c_m.ExpressionLength(length);
        }


    };

    template<class C>
    inline const MatrixLog10<C> log10(const atl::MatrixExpression<typename C::RET_TYPE,C>& c) {
        return MatrixLog10<C>(c);
    }
 template<class C>
    struct MatrixSin : atl::MatrixExpression<typename C::RET_TYPE, MatrixSin<C> > {
        const C& c_m;
        typedef typename C::RET_TYPE RET_TYPE;
        typedef typename C::BASE_TYPE BASE_TYPE;

        inline explicit MatrixSin(const atl::MatrixExpression<typename C::RET_TYPE,C>& c) : c_m(c.Cast()) {
        inline explicit MatrixSin(const atl::MatrixExpression<typename C::RET_TYPE,C>& c) : c_m(c.Cast()) {

        }

        inline const size_t Size(const int32_t & dimension) const {
            switch (dimension) {
                case 0:
                    return c_m.Size(1);
                case 1:
                    return c_m.Size(0);
                default:
                    return 0;
            }
        }

        inline const RET_TYPE operator()(size_t i, size_t j) const {
#ifdef ENABLE_BOUNDS_CHECKING
            assert(i < c_m.Size(0));
            assert(j < c_m.Size(1));
#endif   

            return ::sin(c_m(i,j));

        }
        
          
        /*
         *
         *Returns the first valid index for this matrix given dimension d. 
         */
        inline const size_t IndexMin(uint32_t d) const {
            return c_m.IndexMin(d);
        }

        /*
         *
         *Returns the last valid index for this matrix given dimension d. 
         */
        inline const size_t IndexMax(uint32_t d) const {
            return c_m.IndexMax(d);
        }

        /**
         * Get a value based on the raw indices for the underlying 
         * data. 
         * @param i
         * @return 
         */
        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j) const {
            return ::sin(c_m.AtRaw(i,j));
        }

        /**
         * Used to determine if the resulting matrix is used in the expression.
         * 
         * @param aliased
         * @param ptr
         */
        inline void IsAliased(bool& aliased, void* ptr) const {
            c_m.IsAliased(aliased, ptr);
        }

        inline void ExpressionLength(uint32_t& length)const {
            length++;
            c_m.ExpressionLength(length);
        }


    };

    template<class C>
    inline const MatrixSin<C> sin(const atl::MatrixExpression<typename C::RET_TYPE,C>& c) {
        return MatrixSin<C>(c);
    }
 template<class C>
    struct MatrixSinh : atl::MatrixExpression<typename C::RET_TYPE, MatrixSinh<C> > {
        const C& c_m;
        typedef typename C::RET_TYPE RET_TYPE;
        typedef typename C::BASE_TYPE BASE_TYPE;

        inline explicit MatrixSinh(const atl::MatrixExpression<typename C::RET_TYPE,C>& c) : c_m(c.Cast()) {
        inline explicit MatrixSinh(const atl::MatrixExpression<typename C::RET_TYPE,C>& c) : c_m(c.Cast()) {

        }

        inline const size_t Size(const int32_t & dimension) const {
            switch (dimension) {
                case 0:
                    return c_m.Size(1);
                case 1:
                    return c_m.Size(0);
                default:
                    return 0;
            }
        }

        inline const RET_TYPE operator()(size_t i, size_t j) const {
#ifdef ENABLE_BOUNDS_CHECKING
            assert(i < c_m.Size(0));
            assert(j < c_m.Size(1));
#endif   

            return ::sinh(c_m(i,j));

        }
        
          
        /*
         *
         *Returns the first valid index for this matrix given dimension d. 
         */
        inline const size_t IndexMin(uint32_t d) const {
            return c_m.IndexMin(d);
        }

        /*
         *
         *Returns the last valid index for this matrix given dimension d. 
         */
        inline const size_t IndexMax(uint32_t d) const {
            return c_m.IndexMax(d);
        }

        /**
         * Get a value based on the raw indices for the underlying 
         * data. 
         * @param i
         * @return 
         */
        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j) const {
            return ::sinh(c_m.AtRaw(i,j));
        }

        /**
         * Used to determine if the resulting matrix is used in the expression.
         * 
         * @param aliased
         * @param ptr
         */
        inline void IsAliased(bool& aliased, void* ptr) const {
            c_m.IsAliased(aliased, ptr);
        }

        inline void ExpressionLength(uint32_t& length)const {
            length++;
            c_m.ExpressionLength(length);
        }


    };

    template<class C>
    inline const MatrixSinh<C> sinh(const atl::MatrixExpression<typename C::RET_TYPE,C>& c) {
        return MatrixSinh<C>(c);
    }
 template<class C>
    struct MatrixSqrt : atl::MatrixExpression<typename C::RET_TYPE, MatrixSqrt<C> > {
        const C& c_m;
        typedef typename C::RET_TYPE RET_TYPE;
        typedef typename C::BASE_TYPE BASE_TYPE;

        inline explicit MatrixSqrt(const atl::MatrixExpression<typename C::RET_TYPE,C>& c) : c_m(c.Cast()) {
        inline explicit MatrixSqrt(const atl::MatrixExpression<typename C::RET_TYPE,C>& c) : c_m(c.Cast()) {

        }

        inline const size_t Size(const int32_t & dimension) const {
            switch (dimension) {
                case 0:
                    return c_m.Size(1);
                case 1:
                    return c_m.Size(0);
                default:
                    return 0;
            }
        }

        inline const RET_TYPE operator()(size_t i, size_t j) const {
#ifdef ENABLE_BOUNDS_CHECKING
            assert(i < c_m.Size(0));
            assert(j < c_m.Size(1));
#endif   

            return ::sqrt(c_m(i,j));

        }
        
          
        /*
         *
         *Returns the first valid index for this matrix given dimension d. 
         */
        inline const size_t IndexMin(uint32_t d) const {
            return c_m.IndexMin(d);
        }

        /*
         *
         *Returns the last valid index for this matrix given dimension d. 
         */
        inline const size_t IndexMax(uint32_t d) const {
            return c_m.IndexMax(d);
        }

        /**
         * Get a value based on the raw indices for the underlying 
         * data. 
         * @param i
         * @return 
         */
        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j) const {
            return ::sqrt(c_m.AtRaw(i,j));
        }

        /**
         * Used to determine if the resulting matrix is used in the expression.
         * 
         * @param aliased
         * @param ptr
         */
        inline void IsAliased(bool& aliased, void* ptr) const {
            c_m.IsAliased(aliased, ptr);
        }

        inline void ExpressionLength(uint32_t& length)const {
            length++;
            c_m.ExpressionLength(length);
        }


    };

    template<class C>
    inline const MatrixSqrt<C> sqrt(const atl::MatrixExpression<typename C::RET_TYPE,C>& c) {
        return MatrixSqrt<C>(c);
    }
 template<class C>
    struct MatrixTan : atl::MatrixExpression<typename C::RET_TYPE, MatrixTan<C> > {
        const C& c_m;
        typedef typename C::RET_TYPE RET_TYPE;
        typedef typename C::BASE_TYPE BASE_TYPE;

        inline explicit MatrixTan(const atl::MatrixExpression<typename C::RET_TYPE,C>& c) : c_m(c.Cast()) {
        inline explicit MatrixTan(const atl::MatrixExpression<typename C::RET_TYPE,C>& c) : c_m(c.Cast()) {

        }

        inline const size_t Size(const int32_t & dimension) const {
            switch (dimension) {
                case 0:
                    return c_m.Size(1);
                case 1:
                    return c_m.Size(0);
                default:
                    return 0;
            }
        }

        inline const RET_TYPE operator()(size_t i, size_t j) const {
#ifdef ENABLE_BOUNDS_CHECKING
            assert(i < c_m.Size(0));
            assert(j < c_m.Size(1));
#endif   

            return ::tan(c_m(i,j));

        }
        
          
        /*
         *
         *Returns the first valid index for this matrix given dimension d. 
         */
        inline const size_t IndexMin(uint32_t d) const {
            return c_m.IndexMin(d);
        }

        /*
         *
         *Returns the last valid index for this matrix given dimension d. 
         */
        inline const size_t IndexMax(uint32_t d) const {
            return c_m.IndexMax(d);
        }

        /**
         * Get a value based on the raw indices for the underlying 
         * data. 
         * @param i
         * @return 
         */
        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j) const {
            return ::tan(c_m.AtRaw(i,j));
        }

        /**
         * Used to determine if the resulting matrix is used in the expression.
         * 
         * @param aliased
         * @param ptr
         */
        inline void IsAliased(bool& aliased, void* ptr) const {
            c_m.IsAliased(aliased, ptr);
        }

        inline void ExpressionLength(uint32_t& length)const {
            length++;
            c_m.ExpressionLength(length);
        }


    };

    template<class C>
    inline const MatrixTan<C> tan(const atl::MatrixExpression<typename C::RET_TYPE,C>& c) {
        return MatrixTan<C>(c);
    }
 template<class C>
    struct MatrixTanh : atl::MatrixExpression<typename C::RET_TYPE, MatrixTanh<C> > {
        const C& c_m;
        typedef typename C::RET_TYPE RET_TYPE;
        typedef typename C::BASE_TYPE BASE_TYPE;

        inline explicit MatrixTanh(const atl::MatrixExpression<typename C::RET_TYPE,C>& c) : c_m(c.Cast()) {
        inline explicit MatrixTanh(const atl::MatrixExpression<typename C::RET_TYPE,C>& c) : c_m(c.Cast()) {

        }

        inline const size_t Size(const int32_t & dimension) const {
            switch (dimension) {
                case 0:
                    return c_m.Size(1);
                case 1:
                    return c_m.Size(0);
                default:
                    return 0;
            }
        }

        inline const RET_TYPE operator()(size_t i, size_t j) const {
#ifdef ENABLE_BOUNDS_CHECKING
            assert(i < c_m.Size(0));
            assert(j < c_m.Size(1));
#endif   

            return ::tanh(c_m(i,j));

        }
        
          
        /*
         *
         *Returns the first valid index for this matrix given dimension d. 
         */
        inline const size_t IndexMin(uint32_t d) const {
            return c_m.IndexMin(d);
        }

        /*
         *
         *Returns the last valid index for this matrix given dimension d. 
         */
        inline const size_t IndexMax(uint32_t d) const {
            return c_m.IndexMax(d);
        }

        /**
         * Get a value based on the raw indices for the underlying 
         * data. 
         * @param i
         * @return 
         */
        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j) const {
            return ::tanh(c_m.AtRaw(i,j));
        }

        /**
         * Used to determine if the resulting matrix is used in the expression.
         * 
         * @param aliased
         * @param ptr
         */
        inline void IsAliased(bool& aliased, void* ptr) const {
            c_m.IsAliased(aliased, ptr);
        }

        inline void ExpressionLength(uint32_t& length)const {
            length++;
            c_m.ExpressionLength(length);
        }


    };

    template<class C>
    inline const MatrixTanh<C> tanh(const atl::MatrixExpression<typename C::RET_TYPE,C>& c) {
        return MatrixTanh<C>(c);
    }
 template< class LHS, class RHS>
    struct MatrixPow : MatrixExpression<typename atl::PromoteType<typename LHS::RET_TYPE, typename RHS::RET_TYPE >::return_type, MatrixPow<LHS, RHS> > {
        typedef typename LHS::RET_TYPE RET_TYPEL;
        typedef typename RHS::RET_TYPE RET_TYPER;

        typedef typename atl::PromoteType<RET_TYPEL, RET_TYPER>::return_type RET_TYPE;
        typedef typename atl::PromoteType<typename LHS::BASE_TYPE, typename RHS::BASE_TYPE>::return_type BASE_TYPE;

        const LHS& lhs_m;
        const RHS& rhs_m;

        inline explicit MatrixPow(const atl::MatrixExpression<typename LHS::RET_TYPE,LHS>& lhs, const atl::MatrixExpression<typename RHS::RET_TYPE,RHS>& rhs) : lhs_m(lhs.Cast()), rhs_m(rhs.Cast()) {



        }

        inline const size_t Size(const int32_t & dimension) const {
            return lhs_m.Size(dimension) < rhs_m.Size(dimension) ? lhs_m.Size(dimension) : rhs_m.Size(dimension);
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j) const {
            return ::pow(lhs_m(i, j), rhs_m(i, j));
        }

        /*
         *
         *Returns the first valid index for this matrix given dimension d. 
         */
        inline const size_t IndexMin(uint32_t d) const {
            return std::max(lhs_m.IndexMin(d), rhs_m.IndexMin(d));
        }

        /*
         *
         *Returns the last valid index for this matrix given dimension d. 
         */
        inline const size_t IndexMax(uint32_t d) const {
            return std::min(lhs_m.IndexMax(d), rhs_m.IndexMax(d));
        }

        /**
         * Get a value based on the raw indices for the underlying 
         * data. 
         * @param i
         * @return 
         */
        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j) const {
            return ::pow(lhs_m.AtRaw(i, j), rhs_m.AtRaw(i, j));
        }

        /**
         * Used to determine if the resulting matrix is used in the expression.
         * 
         * @param aliased
         * @param ptr
         */
        inline void IsAliased(bool& aliased, void* ptr) const {
            lhs_m.IsAliased(aliased, ptr);
            rhs_m.IsAliased(aliased, ptr);
        }

        inline void ExpressionLength(uint32_t& length)const {
            length++;
            lhs_m.ExpressionLength(length);
            rhs_m.ExpressionLength(length);
        }


    };
    //

    template< class LHS, class T >
    struct MatrixPowScalar : MatrixExpression<typename PromoteType<typename LHS::RET_TYPE, T>::return_type, MatrixPowScalar< LHS, T> > {
        const LHS& lhs_m;
        const T& rhs_m;

        typedef typename atl::PromoteType<typename LHS::RET_TYPE, T>::return_type RET_TYPE;
        typedef typename atl::PromoteType<typename LHS::BASE_TYPE, T>::return_type BASE_TYPE;

        inline explicit MatrixPowScalar(const MatrixExpression<typename LHS::RET_TYPE, LHS>& lhs, const T & rhs) : lhs_m(lhs.Cast()), rhs_m(rhs) {



        }

        inline const size_t Size(const int32_t & dimension) const {
            return lhs_m.Size(dimension);
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j) const {
            return ::pow(lhs_m(i, j),rhs_m);
        }

        /*
         *
         *Returns the first valid index for this matrix given dimension d. 
         */
        inline const size_t IndexMin(uint32_t d) const {
            return (lhs_m.IndexMin(d));
        }

        /*
         *
         *Returns the last valid index for this matrix given dimension d. 
         */
        inline const size_t IndexMax(uint32_t d) const {
            return (lhs_m.IndexMax(d));
        }

        /**
         * Get a value based on the raw indices for the underlying 
         * data. 
         * @param i
         * @return 
         */
        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j) const {
            return ::pow(lhs_m.AtRaw(i, j) , rhs_m);
        }

        /**
         * Used to determine if the resulting matrix is used in the expression.
         * 
         * @param aliased
         * @param ptr
         */
        inline void IsAliased(bool& aliased, void* ptr) const{
            lhs_m.IsAliased(aliased, ptr);
        }

        inline void ExpressionLength(uint32_t& length)const {
            length++;
            lhs_m.ExpressionLength(length);
        }


    };

    template<class T, class RHS>
    struct ScalarMatrixPow : MatrixExpression<typename PromoteType<typename RHS::RET_TYPE, T>::return_type, ScalarMatrixPow<T, RHS> > {
        typedef typename atl::PromoteType<T, typename RHS::RET_TYPE>::return_type RET_TYPE;
        typedef typename atl::PromoteType<T, typename RHS::BASE_TYPE>::return_type BASE_TYPE;

        const T& lhs_m;
        const RHS& rhs_m;

        inline explicit ScalarMatrixPow(const T& lhs, const MatrixExpression<typename RHS::RET_TYPE, RHS> & rhs) : lhs_m(lhs), rhs_m(rhs.Cast()) {



        }

        inline const size_t Size(const int32_t & dimension) const {
            return rhs_m.Size(dimension);
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j) const {
            return ::pow(lhs_m, rhs_m(i, j));
        }

        /*
         *
         *Returns the first valid index for this matrix given dimension d. 
         */
        inline const size_t IndexMin(uint32_t d) const {
            return (rhs_m.IndexMin(d));
        }

        /*
         *
         *Returns the last valid index for this matrix given dimension d. 
         */
        inline const size_t IndexMax(uint32_t d) const {
            return (rhs_m.IndexMax(d));
        }

        /**
         * Get a value based on the raw indices for the underlying 
         * data. 
         * @param i
         * @return 
         */
        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j) const {
            return ::pow(lhs_m, rhs_m.AtRaw(i, j));
        }

        /**
         * Used to determine if the resulting matrix is used in the expression.
         * 
         * @param aliased
         * @param ptr
         */
        inline void IsAliased(bool& aliased, void* ptr) const{
            rhs_m.IsAliased(aliased, ptr);
        }

        inline void ExpressionLength(uint32_t& length)const {
            length++;
            rhs_m.ExpressionLength(length);
        }



    };
    //
    //typename et4ad::promote_trait<typename LHS::RET_TYPE , typename RHS::RET_TYPE >::return_type

    template <class LHS, class RHS>
    inline const MatrixPow< LHS, RHS> operator+(const MatrixExpression<typename LHS::RET_TYPE, LHS>& a,
            const MatrixExpression<typename RHS::RET_TYPE, RHS>& b) {
        return MatrixPow< LHS, RHS > (a.Cast(), b.Cast());
    }



#define ATL_MATRIX_MatrixPow_SCALAR(TYPE) \
    template< class LHS>      \
        inline const MatrixPowScalar<LHS,TYPE> \
    pow(const MatrixExpression<typename LHS::RET_TYPE, LHS>& a, const TYPE& b) {\
            return MatrixPowScalar<LHS,TYPE > (a.Cast(), b);\
        } \
    
    ATL_MATRIX_MatrixPow_SCALAR(short)
    ATL_MATRIX_MatrixPow_SCALAR(unsigned short)
    ATL_MATRIX_MatrixPow_SCALAR(int)
    ATL_MATRIX_MatrixPow_SCALAR(unsigned int)
    ATL_MATRIX_MatrixPow_SCALAR(long)
    ATL_MATRIX_MatrixPow_SCALAR(unsigned long)
    ATL_MATRIX_MatrixPow_SCALAR(float)
    ATL_MATRIX_MatrixPow_SCALAR(double)
    ATL_MATRIX_MatrixPow_SCALAR(long double)
    ATL_MATRIX_MatrixPow_SCALAR(atl::Variable<float>)
    ATL_MATRIX_MatrixPow_SCALAR(atl::Variable<double>)
    ATL_MATRIX_MatrixPow_SCALAR(atl::Variable<long double>)


#define ATL_MatrixPow_SCALAR_MATRIX(TYPE) \
    template< class RHS>      \
        inline const ScalarMatrixPow<TYPE,RHS> \
    pow(const TYPE& a, const MatrixExpression<typename RHS::RET_TYPE, RHS>& b ) {\
            return ScalarMatrixPow<TYPE,RHS > (a, b.Cast());\
        } \
            
    ATL_MatrixPow_SCALAR_MATRIX(short)
    ATL_MatrixPow_SCALAR_MATRIX(unsigned short)
    ATL_MatrixPow_SCALAR_MATRIX(int)
    ATL_MatrixPow_SCALAR_MATRIX(unsigned int)
    ATL_MatrixPow_SCALAR_MATRIX(long)
    ATL_MatrixPow_SCALAR_MATRIX(unsigned long)
    ATL_MatrixPow_SCALAR_MATRIX(float)
    ATL_MatrixPow_SCALAR_MATRIX(double)
    ATL_MatrixPow_SCALAR_MATRIX(long double)
    ATL_MatrixPow_SCALAR_MATRIX(atl::Variable<float>)
    ATL_MatrixPow_SCALAR_MATRIX(atl::Variable<double>)
    ATL_MatrixPow_SCALAR_MATRIX(atl::Variable<long double>)



 template< class LHS, class RHS>
    struct MatrixATan2 : MatrixExpression<typename atl::PromoteType<typename LHS::RET_TYPE, typename RHS::RET_TYPE >::return_type, MatrixATan2<LHS, RHS> > {
        typedef typename LHS::RET_TYPE RET_TYPEL;
        typedef typename RHS::RET_TYPE RET_TYPER;

        typedef typename atl::PromoteType<RET_TYPEL, RET_TYPER>::return_type RET_TYPE;
        typedef typename atl::PromoteType<typename LHS::BASE_TYPE, typename RHS::BASE_TYPE>::return_type BASE_TYPE;

        const LHS& lhs_m;
        const RHS& rhs_m;

        inline explicit MatrixATan2(const atl::MatrixExpression<typename LHS::RET_TYPE,LHS>& lhs, const atl::MatrixExpression<typename RHS::RET_TYPE,RHS>& rhs) : lhs_m(lhs.Cast()), rhs_m(rhs.Cast()) {



        }

        inline const size_t Size(const int32_t & dimension) const {
            return lhs_m.Size(dimension) < rhs_m.Size(dimension) ? lhs_m.Size(dimension) : rhs_m.Size(dimension);
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j) const {
            return ::atan2(lhs_m(i, j), rhs_m(i, j));
        }

        /*
         *
         *Returns the first valid index for this matrix given dimension d. 
         */
        inline const size_t IndexMin(uint32_t d) const {
            return std::max(lhs_m.IndexMin(d), rhs_m.IndexMin(d));
        }

        /*
         *
         *Returns the last valid index for this matrix given dimension d. 
         */
        inline const size_t IndexMax(uint32_t d) const {
            return std::min(lhs_m.IndexMax(d), rhs_m.IndexMax(d));
        }

        /**
         * Get a value based on the raw indices for the underlying 
         * data. 
         * @param i
         * @return 
         */
        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j) const {
            return ::atan2(lhs_m.AtRaw(i, j), rhs_m.AtRaw(i, j));
        }

        /**
         * Used to determine if the resulting matrix is used in the expression.
         * 
         * @param aliased
         * @param ptr
         */
        inline void IsAliased(bool& aliased, void* ptr) const {
            lhs_m.IsAliased(aliased, ptr);
            rhs_m.IsAliased(aliased, ptr);
        }

        inline void ExpressionLength(uint32_t& length)const {
            length++;
            lhs_m.ExpressionLength(length);
            rhs_m.ExpressionLength(length);
        }


    };
    //

    template< class LHS, class T >
    struct MatrixATan2Scalar : MatrixExpression<typename PromoteType<typename LHS::RET_TYPE, T>::return_type, MatrixATan2Scalar< LHS, T> > {
        const LHS& lhs_m;
        const T& rhs_m;

        typedef typename atl::PromoteType<typename LHS::RET_TYPE, T>::return_type RET_TYPE;
        typedef typename atl::PromoteType<typename LHS::BASE_TYPE, T>::return_type BASE_TYPE;

        inline explicit MatrixATan2Scalar(const MatrixExpression<typename LHS::RET_TYPE, LHS>& lhs, const T & rhs) : lhs_m(lhs.Cast()), rhs_m(rhs) {



        }

        inline const size_t Size(const int32_t & dimension) const {
            return lhs_m.Size(dimension);
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j) const {
            return ::atan2(lhs_m(i, j),rhs_m);
        }

        /*
         *
         *Returns the first valid index for this matrix given dimension d. 
         */
        inline const size_t IndexMin(uint32_t d) const {
            return (lhs_m.IndexMin(d));
        }

        /*
         *
         *Returns the last valid index for this matrix given dimension d. 
         */
        inline const size_t IndexMax(uint32_t d) const {
            return (lhs_m.IndexMax(d));
        }

        /**
         * Get a value based on the raw indices for the underlying 
         * data. 
         * @param i
         * @return 
         */
        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j) const {
            return ::atan2(lhs_m.AtRaw(i, j) , rhs_m);
        }

        /**
         * Used to determine if the resulting matrix is used in the expression.
         * 
         * @param aliased
         * @param ptr
         */
        inline void IsAliased(bool& aliased, void* ptr) const{
            lhs_m.IsAliased(aliased, ptr);
        }

        inline void ExpressionLength(uint32_t& length)const {
            length++;
            lhs_m.ExpressionLength(length);
        }


    };

    template<class T, class RHS>
    struct ScalarMatrixATan2 : MatrixExpression<typename PromoteType<typename RHS::RET_TYPE, T>::return_type, ScalarMatrixATan2<T, RHS> > {
        typedef typename atl::PromoteType<T, typename RHS::RET_TYPE>::return_type RET_TYPE;
        typedef typename atl::PromoteType<T, typename RHS::BASE_TYPE>::return_type BASE_TYPE;

        const T& lhs_m;
        const RHS& rhs_m;

        inline explicit ScalarMatrixATan2(const T& lhs, const MatrixExpression<typename RHS::RET_TYPE, RHS> & rhs) : lhs_m(lhs), rhs_m(rhs.Cast()) {



        }

        inline const size_t Size(const int32_t & dimension) const {
            return rhs_m.Size(dimension);
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j) const {
            return ::atan2(lhs_m, rhs_m(i, j));
        }

        /*
         *
         *Returns the first valid index for this matrix given dimension d. 
         */
        inline const size_t IndexMin(uint32_t d) const {
            return (rhs_m.IndexMin(d));
        }

        /*
         *
         *Returns the last valid index for this matrix given dimension d. 
         */
        inline const size_t IndexMax(uint32_t d) const {
            return (rhs_m.IndexMax(d));
        }

        /**
         * Get a value based on the raw indices for the underlying 
         * data. 
         * @param i
         * @return 
         */
        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j) const {
            return ::atan2(lhs_m, rhs_m.AtRaw(i, j));
        }

        /**
         * Used to determine if the resulting matrix is used in the expression.
         * 
         * @param aliased
         * @param ptr
         */
        inline void IsAliased(bool& aliased, void* ptr) const{
            rhs_m.IsAliased(aliased, ptr);
        }

        inline void ExpressionLength(uint32_t& length)const {
            length++;
            rhs_m.ExpressionLength(length);
        }



    };
    //
    //typename et4ad::promote_trait<typename LHS::RET_TYPE , typename RHS::RET_TYPE >::return_type

    template <class LHS, class RHS>
    inline const MatrixATan2< LHS, RHS> operator+(const MatrixExpression<typename LHS::RET_TYPE, LHS>& a,
            const MatrixExpression<typename RHS::RET_TYPE, RHS>& b) {
        return MatrixATan2< LHS, RHS > (a.Cast(), b.Cast());
    }



#define ATL_MATRIX_MatrixATan2_SCALAR(TYPE) \
    template< class LHS>      \
        inline const MatrixATan2Scalar<LHS,TYPE> \
    atan2(const MatrixExpression<typename LHS::RET_TYPE, LHS>& a, const TYPE& b) {\
            return MatrixATan2Scalar<LHS,TYPE > (a.Cast(), b);\
        } \
    
    ATL_MATRIX_MatrixATan2_SCALAR(short)
    ATL_MATRIX_MatrixATan2_SCALAR(unsigned short)
    ATL_MATRIX_MatrixATan2_SCALAR(int)
    ATL_MATRIX_MatrixATan2_SCALAR(unsigned int)
    ATL_MATRIX_MatrixATan2_SCALAR(long)
    ATL_MATRIX_MatrixATan2_SCALAR(unsigned long)
    ATL_MATRIX_MatrixATan2_SCALAR(float)
    ATL_MATRIX_MatrixATan2_SCALAR(double)
    ATL_MATRIX_MatrixATan2_SCALAR(long double)
    ATL_MATRIX_MatrixATan2_SCALAR(atl::Variable<float>)
    ATL_MATRIX_MatrixATan2_SCALAR(atl::Variable<double>)
    ATL_MATRIX_MatrixATan2_SCALAR(atl::Variable<long double>)


#define ATL_MatrixATan2_SCALAR_MATRIX(TYPE) \
    template< class RHS>      \
        inline const ScalarMatrixATan2<TYPE,RHS> \
    atan2(const TYPE& a, const MatrixExpression<typename RHS::RET_TYPE, RHS>& b ) {\
            return ScalarMatrixATan2<TYPE,RHS > (a, b.Cast());\
        } \
            
    ATL_MatrixATan2_SCALAR_MATRIX(short)
    ATL_MatrixATan2_SCALAR_MATRIX(unsigned short)
    ATL_MatrixATan2_SCALAR_MATRIX(int)
    ATL_MatrixATan2_SCALAR_MATRIX(unsigned int)
    ATL_MatrixATan2_SCALAR_MATRIX(long)
    ATL_MatrixATan2_SCALAR_MATRIX(unsigned long)
    ATL_MatrixATan2_SCALAR_MATRIX(float)
    ATL_MatrixATan2_SCALAR_MATRIX(double)
    ATL_MatrixATan2_SCALAR_MATRIX(long double)
    ATL_MatrixATan2_SCALAR_MATRIX(atl::Variable<float>)
    ATL_MatrixATan2_SCALAR_MATRIX(atl::Variable<double>)
    ATL_MatrixATan2_SCALAR_MATRIX(atl::Variable<long double>)





}
#endif