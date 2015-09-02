#ifndef VECTORMATH_HPP
#define VECTORMATH_HPP

#ifdef ENABLE_BOUNDS_CHECKING
#include <assert.h>
#endif




namespace atl {

    template<class C>
    struct VectorACos : atl::VectorExpression<typename C::RET_TYPE, VectorACos<C> > {
        const C& c_m;
        typedef typename C::RET_TYPE RET_TYPE;
        typedef typename C::BASE_TYPE BASE_TYPE;

        inline explicit VectorACos(const atl::VectorExpression<typename C::RET_TYPE, C>& c) : c_m(c.Cast()) {

        }

        inline const size_t Size(const int32_t & dimension) const {
            switch (dimension) {
                case 0:
                    return c_m.Size(0);
                default:
                    return 0;
            }
        }

        inline const RET_TYPE operator()(size_t i) const {
#ifdef ENABLE_BOUNDS_CHECKING
            assert(i < c_m.Size(0));
#endif   

            return ::acos(c_m(i));

        }

        /*
         *
         *Returns the first valid index for this container given dimension d. 
         */
        inline const size_t IndexMin(uint32_t d) const {
            return c_m.IndexMin(d);
        }

        /*
         *
         *Returns the last valid index for this container given dimension d. 
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
        inline const RET_TYPE AtRaw(const uint32_t& i) const {
            return ::acos(c_m.AtRaw(i));
        }

        /**
         * Used to determine if the resulting container is used in the expression.
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
    inline const VectorACos<C> acos(const atl::VectorExpression<typename C::RET_TYPE, C>& c) {
        return VectorACos<C>(c);
    }

    template<class C>
    struct VectorASin : atl::VectorExpression<typename C::RET_TYPE, VectorASin<C> > {
        const C& c_m;
        typedef typename C::RET_TYPE RET_TYPE;
        typedef typename C::BASE_TYPE BASE_TYPE;

        inline explicit VectorASin(const atl::VectorExpression<typename C::RET_TYPE, C>& c) : c_m(c.Cast()) {

        }

        inline const size_t Size(const int32_t & dimension) const {
            switch (dimension) {
                case 0:
                    return c_m.Size(0);
                default:
                    return 0;
            }
        }

        inline const RET_TYPE operator()(size_t i) const {
#ifdef ENABLE_BOUNDS_CHECKING
            assert(i < c_m.Size(0));
#endif   

            return ::asin(c_m(i));

        }

        /*
         *
         *Returns the first valid index for this container given dimension d. 
         */
        inline const size_t IndexMin(uint32_t d) const {
            return c_m.IndexMin(d);
        }

        /*
         *
         *Returns the last valid index for this container given dimension d. 
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
        inline const RET_TYPE AtRaw(const uint32_t& i) const {
            return ::asin(c_m.AtRaw(i));
        }

        /**
         * Used to determine if the resulting container is used in the expression.
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
    inline const VectorASin<C> asin(const atl::VectorExpression<typename C::RET_TYPE, C>& c) {
        return VectorASin<C>(c);
    }

    template<class C>
    struct VectorATan : atl::VectorExpression<typename C::RET_TYPE, VectorATan<C> > {
        const C& c_m;
        typedef typename C::RET_TYPE RET_TYPE;
        typedef typename C::BASE_TYPE BASE_TYPE;

        inline explicit VectorATan(const atl::VectorExpression<typename C::RET_TYPE, C>& c) : c_m(c.Cast()) {

        }

        inline const size_t Size(const int32_t & dimension) const {
            switch (dimension) {
                case 0:
                    return c_m.Size(0);
                default:
                    return 0;
            }
        }

        inline const RET_TYPE operator()(size_t i) const {
#ifdef ENABLE_BOUNDS_CHECKING
            assert(i < c_m.Size(0));
#endif   

            return ::atan(c_m(i));

        }

        /*
         *
         *Returns the first valid index for this container given dimension d. 
         */
        inline const size_t IndexMin(uint32_t d) const {
            return c_m.IndexMin(d);
        }

        /*
         *
         *Returns the last valid index for this container given dimension d. 
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
        inline const RET_TYPE AtRaw(const uint32_t& i) const {
            return ::atan(c_m.AtRaw(i));
        }

        /**
         * Used to determine if the resulting container is used in the expression.
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
    inline const VectorATan<C> atan(const atl::VectorExpression<typename C::RET_TYPE, C>& c) {
        return VectorATan<C>(c);
    }

    template<class C>
    struct VectorACeil : atl::VectorExpression<typename C::RET_TYPE, VectorACeil<C> > {
        const C& c_m;
        typedef typename C::RET_TYPE RET_TYPE;
        typedef typename C::BASE_TYPE BASE_TYPE;

        inline explicit VectorACeil(const atl::VectorExpression<typename C::RET_TYPE, C>& c) : c_m(c.Cast()) {

        }

        inline const size_t Size(const int32_t & dimension) const {
            switch (dimension) {
                case 0:
                    return c_m.Size(0);
                default:
                    return 0;
            }
        }

        inline const RET_TYPE operator()(size_t i) const {
#ifdef ENABLE_BOUNDS_CHECKING
            assert(i < c_m.Size(0));
#endif   

            return ::ceil(c_m(i));

        }

        /*
         *
         *Returns the first valid index for this container given dimension d. 
         */
        inline const size_t IndexMin(uint32_t d) const {
            return c_m.IndexMin(d);
        }

        /*
         *
         *Returns the last valid index for this container given dimension d. 
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
        inline const RET_TYPE AtRaw(const uint32_t& i) const {
            return ::ceil(c_m.AtRaw(i));
        }

        /**
         * Used to determine if the resulting container is used in the expression.
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
    inline const VectorACeil<C> ceil(const atl::VectorExpression<typename C::RET_TYPE, C>& c) {
        return VectorACeil<C>(c);
    }

    template<class C>
    struct VectorCos : atl::VectorExpression<typename C::RET_TYPE, VectorCos<C> > {
        const C& c_m;
        typedef typename C::RET_TYPE RET_TYPE;
        typedef typename C::BASE_TYPE BASE_TYPE;

        inline explicit VectorCos(const atl::VectorExpression<typename C::RET_TYPE, C>& c) : c_m(c.Cast()) {

        }

        inline const size_t Size(const int32_t & dimension) const {
            switch (dimension) {
                case 0:
                    return c_m.Size(0);
                default:
                    return 0;
            }
        }

        inline const RET_TYPE operator()(size_t i) const {
#ifdef ENABLE_BOUNDS_CHECKING
            assert(i < c_m.Size(0));
#endif   

            return ::cos(c_m(i));

        }

        /*
         *
         *Returns the first valid index for this container given dimension d. 
         */
        inline const size_t IndexMin(uint32_t d) const {
            return c_m.IndexMin(d);
        }

        /*
         *
         *Returns the last valid index for this container given dimension d. 
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
        inline const RET_TYPE AtRaw(const uint32_t& i) const {
            return ::cos(c_m.AtRaw(i));
        }

        /**
         * Used to determine if the resulting container is used in the expression.
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
    inline const VectorCos<C> cos(const atl::VectorExpression<typename C::RET_TYPE, C>& c) {
        return VectorCos<C>(c);
    }

    template<class C>
    struct VectorCosh : atl::VectorExpression<typename C::RET_TYPE, VectorCosh<C> > {
        const C& c_m;
        typedef typename C::RET_TYPE RET_TYPE;
        typedef typename C::BASE_TYPE BASE_TYPE;

        inline explicit VectorCosh(const atl::VectorExpression<typename C::RET_TYPE, C>& c) : c_m(c.Cast()) {

        }

        inline const size_t Size(const int32_t & dimension) const {
            switch (dimension) {
                case 0:
                    return c_m.Size(0);
                default:
                    return 0;
            }
        }

        inline const RET_TYPE operator()(size_t i) const {
#ifdef ENABLE_BOUNDS_CHECKING
            assert(i < c_m.Size(0));
#endif   

            return ::cosh(c_m(i));

        }

        /*
         *
         *Returns the first valid index for this container given dimension d. 
         */
        inline const size_t IndexMin(uint32_t d) const {
            return c_m.IndexMin(d);
        }

        /*
         *
         *Returns the last valid index for this container given dimension d. 
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
        inline const RET_TYPE AtRaw(const uint32_t& i) const {
            return ::cosh(c_m.AtRaw(i));
        }

        /**
         * Used to determine if the resulting container is used in the expression.
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
    inline const VectorCosh<C> cosh(const atl::VectorExpression<typename C::RET_TYPE, C>& c) {
        return VectorCosh<C>(c);
    }

    template<class C>
    struct VectorExp : atl::VectorExpression<typename C::RET_TYPE, VectorExp<C> > {
        const C& c_m;
        typedef typename C::RET_TYPE RET_TYPE;
        typedef typename C::BASE_TYPE BASE_TYPE;

        inline explicit VectorExp(const atl::VectorExpression<typename C::RET_TYPE, C>& c) : c_m(c.Cast()) {

        }

        inline const size_t Size(const int32_t & dimension) const {
            switch (dimension) {
                case 0:
                    return c_m.Size(0);
                default:
                    return 0;
            }
        }

        inline const RET_TYPE operator()(size_t i) const {
#ifdef ENABLE_BOUNDS_CHECKING
            assert(i < c_m.Size(0));
#endif   

            return ::exp(c_m(i));

        }

        /*
         *
         *Returns the first valid index for this container given dimension d. 
         */
        inline const size_t IndexMin(uint32_t d) const {
            return c_m.IndexMin(d);
        }

        /*
         *
         *Returns the last valid index for this container given dimension d. 
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
        inline const RET_TYPE AtRaw(const uint32_t& i) const {
            return ::exp(c_m.AtRaw(i));
        }

        /**
         * Used to determine if the resulting container is used in the expression.
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
    inline const VectorExp<C> exp(const atl::VectorExpression<typename C::RET_TYPE, C>& c) {
        return VectorExp<C>(c);
    }

    template<class C>
    struct VectorFabs : atl::VectorExpression<typename C::RET_TYPE, VectorFabs<C> > {
        const C& c_m;
        typedef typename C::RET_TYPE RET_TYPE;
        typedef typename C::BASE_TYPE BASE_TYPE;

        inline explicit VectorFabs(const atl::VectorExpression<typename C::RET_TYPE, C>& c) : c_m(c.Cast()) {

        }

        inline const size_t Size(const int32_t & dimension) const {
            switch (dimension) {
                case 0:
                    return c_m.Size(0);
                default:
                    return 0;
            }
        }

        inline const RET_TYPE operator()(size_t i) const {
#ifdef ENABLE_BOUNDS_CHECKING
            assert(i < c_m.Size(0));
#endif   

            return ::fabs(c_m(i));

        }

        /*
         *
         *Returns the first valid index for this container given dimension d. 
         */
        inline const size_t IndexMin(uint32_t d) const {
            return c_m.IndexMin(d);
        }

        /*
         *
         *Returns the last valid index for this container given dimension d. 
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
        inline const RET_TYPE AtRaw(const uint32_t& i) const {
            return ::fabs(c_m.AtRaw(i));
        }

        /**
         * Used to determine if the resulting container is used in the expression.
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
    inline const VectorFabs<C> fabs(const atl::VectorExpression<typename C::RET_TYPE, C>& c) {
        return VectorFabs<C>(c);
    }

    template<class C>
    struct VectorFloor : atl::VectorExpression<typename C::RET_TYPE, VectorFloor<C> > {
        const C& c_m;
        typedef typename C::RET_TYPE RET_TYPE;
        typedef typename C::BASE_TYPE BASE_TYPE;

        inline explicit VectorFloor(const atl::VectorExpression<typename C::RET_TYPE, C>& c) : c_m(c.Cast()) {

        }

        inline const size_t Size(const int32_t & dimension) const {
            switch (dimension) {
                case 0:
                    return c_m.Size(0);
                default:
                    return 0;
            }
        }

        inline const RET_TYPE operator()(size_t i) const {
#ifdef ENABLE_BOUNDS_CHECKING
            assert(i < c_m.Size(0));
#endif   

            return ::floor(c_m(i));

        }

        /*
         *
         *Returns the first valid index for this container given dimension d. 
         */
        inline const size_t IndexMin(uint32_t d) const {
            return c_m.IndexMin(d);
        }

        /*
         *
         *Returns the last valid index for this container given dimension d. 
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
        inline const RET_TYPE AtRaw(const uint32_t& i) const {
            return ::floor(c_m.AtRaw(i));
        }

        /**
         * Used to determine if the resulting container is used in the expression.
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
    inline const VectorFloor<C> floor(const atl::VectorExpression<typename C::RET_TYPE, C>& c) {
        return VectorFloor<C>(c);
    }

    template<class C>
    struct VectorLog : atl::VectorExpression<typename C::RET_TYPE, VectorLog<C> > {
        const C& c_m;
        typedef typename C::RET_TYPE RET_TYPE;
        typedef typename C::BASE_TYPE BASE_TYPE;

        inline explicit VectorLog(const atl::VectorExpression<typename C::RET_TYPE, C>& c) : c_m(c.Cast()) {

        }

        inline const size_t Size(const int32_t & dimension) const {
            switch (dimension) {
                case 0:
                    return c_m.Size(0);
                default:
                    return 0;
            }
        }

        inline const RET_TYPE operator()(size_t i) const {
#ifdef ENABLE_BOUNDS_CHECKING
            assert(i < c_m.Size(0));
#endif   

            return ::log(c_m(i));

        }

        /*
         *
         *Returns the first valid index for this container given dimension d. 
         */
        inline const size_t IndexMin(uint32_t d) const {
            return c_m.IndexMin(d);
        }

        /*
         *
         *Returns the last valid index for this container given dimension d. 
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
        inline const RET_TYPE AtRaw(const uint32_t& i) const {
            return ::log(c_m.AtRaw(i));
        }

        /**
         * Used to determine if the resulting container is used in the expression.
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
    inline const VectorLog<C> log(const atl::VectorExpression<typename C::RET_TYPE, C>& c) {
        return VectorLog<C>(c);
    }

    template<class C>
    struct VectorLog10 : atl::VectorExpression<typename C::RET_TYPE, VectorLog10<C> > {
        const C& c_m;
        typedef typename C::RET_TYPE RET_TYPE;
        typedef typename C::BASE_TYPE BASE_TYPE;

        inline explicit VectorLog10(const atl::VectorExpression<typename C::RET_TYPE, C>& c) : c_m(c.Cast()) {

        }

        inline const size_t Size(const int32_t & dimension) const {
            switch (dimension) {
                case 0:
                    return c_m.Size(0);
                default:
                    return 0;
            }
        }

        inline const RET_TYPE operator()(size_t i) const {
#ifdef ENABLE_BOUNDS_CHECKING
            assert(i < c_m.Size(0));
#endif   

            return ::log10(c_m(i));

        }

        /*
         *
         *Returns the first valid index for this container given dimension d. 
         */
        inline const size_t IndexMin(uint32_t d) const {
            return c_m.IndexMin(d);
        }

        /*
         *
         *Returns the last valid index for this container given dimension d. 
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
        inline const RET_TYPE AtRaw(const uint32_t& i) const {
            return ::log10(c_m.AtRaw(i));
        }

        /**
         * Used to determine if the resulting container is used in the expression.
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
    inline const VectorLog10<C> log10(const atl::VectorExpression<typename C::RET_TYPE, C>& c) {
        return VectorLog10<C>(c);
    }

    template<class C>
    struct VectorSin : atl::VectorExpression<typename C::RET_TYPE, VectorSin<C> > {
        const C& c_m;
        typedef typename C::RET_TYPE RET_TYPE;
        typedef typename C::BASE_TYPE BASE_TYPE;

        inline explicit VectorSin(const atl::VectorExpression<typename C::RET_TYPE, C>& c) : c_m(c.Cast()) {

        }

        inline const size_t Size(const int32_t & dimension) const {
            switch (dimension) {
                case 0:
                    return c_m.Size(0);
                default:
                    return 0;
            }
        }

        inline const RET_TYPE operator()(size_t i) const {
#ifdef ENABLE_BOUNDS_CHECKING
            assert(i < c_m.Size(0));
#endif   

            return ::sin(c_m(i));

        }

        /*
         *
         *Returns the first valid index for this container given dimension d. 
         */
        inline const size_t IndexMin(uint32_t d) const {
            return c_m.IndexMin(d);
        }

        /*
         *
         *Returns the last valid index for this container given dimension d. 
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
        inline const RET_TYPE AtRaw(const uint32_t& i) const {
            return ::sin(c_m.AtRaw(i));
        }

        /**
         * Used to determine if the resulting container is used in the expression.
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
    inline const VectorSin<C> sin(const atl::VectorExpression<typename C::RET_TYPE, C>& c) {
        return VectorSin<C>(c);
    }

    template<class C>
    struct VectorSinh : atl::VectorExpression<typename C::RET_TYPE, VectorSinh<C> > {
        const C& c_m;
        typedef typename C::RET_TYPE RET_TYPE;
        typedef typename C::BASE_TYPE BASE_TYPE;

        inline explicit VectorSinh(const atl::VectorExpression<typename C::RET_TYPE, C>& c) : c_m(c.Cast()) {

        }

        inline const size_t Size(const int32_t & dimension) const {
            switch (dimension) {
                case 0:
                    return c_m.Size(0);
                default:
                    return 0;
            }
        }

        inline const RET_TYPE operator()(size_t i) const {
#ifdef ENABLE_BOUNDS_CHECKING
            assert(i < c_m.Size(0));
#endif   

            return ::sinh(c_m(i));

        }

        /*
         *
         *Returns the first valid index for this container given dimension d. 
         */
        inline const size_t IndexMin(uint32_t d) const {
            return c_m.IndexMin(d);
        }

        /*
         *
         *Returns the last valid index for this container given dimension d. 
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
        inline const RET_TYPE AtRaw(const uint32_t& i) const {
            return ::sinh(c_m.AtRaw(i));
        }

        /**
         * Used to determine if the resulting container is used in the expression.
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
    inline const VectorSinh<C> sinh(const atl::VectorExpression<typename C::RET_TYPE, C>& c) {
        return VectorSinh<C>(c);
    }

    template<class C>
    struct VectorSqrt : atl::VectorExpression<typename C::RET_TYPE, VectorSqrt<C> > {
        const C& c_m;
        typedef typename C::RET_TYPE RET_TYPE;
        typedef typename C::BASE_TYPE BASE_TYPE;

        inline explicit VectorSqrt(const atl::VectorExpression<typename C::RET_TYPE, C>& c) : c_m(c.Cast()) {

        }

        inline const size_t Size(const int32_t & dimension) const {
            switch (dimension) {
                case 0:
                    return c_m.Size(0);
                default:
                    return 0;
            }
        }

        inline const RET_TYPE operator()(size_t i) const {
#ifdef ENABLE_BOUNDS_CHECKING
            assert(i < c_m.Size(0));
#endif   

            return ::sqrt(c_m(i));

        }

        /*
         *
         *Returns the first valid index for this container given dimension d. 
         */
        inline const size_t IndexMin(uint32_t d) const {
            return c_m.IndexMin(d);
        }

        /*
         *
         *Returns the last valid index for this container given dimension d. 
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
        inline const RET_TYPE AtRaw(const uint32_t& i) const {
            return ::sqrt(c_m.AtRaw(i));
        }

        /**
         * Used to determine if the resulting container is used in the expression.
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
    inline const VectorSqrt<C> sqrt(const atl::VectorExpression<typename C::RET_TYPE, C>& c) {
        return VectorSqrt<C>(c);
    }

    template<class C>
    struct VectorTan : atl::VectorExpression<typename C::RET_TYPE, VectorTan<C> > {
        const C& c_m;
        typedef typename C::RET_TYPE RET_TYPE;
        typedef typename C::BASE_TYPE BASE_TYPE;

        inline explicit VectorTan(const atl::VectorExpression<typename C::RET_TYPE, C>& c) : c_m(c.Cast()) {

        }

        inline const size_t Size(const int32_t & dimension) const {
            switch (dimension) {
                case 0:
                    return c_m.Size(0);
                default:
                    return 0;
            }
        }

        inline const RET_TYPE operator()(size_t i) const {
#ifdef ENABLE_BOUNDS_CHECKING
            assert(i < c_m.Size(0));
#endif   

            return ::tan(c_m(i));

        }

        /*
         *
         *Returns the first valid index for this container given dimension d. 
         */
        inline const size_t IndexMin(uint32_t d) const {
            return c_m.IndexMin(d);
        }

        /*
         *
         *Returns the last valid index for this container given dimension d. 
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
        inline const RET_TYPE AtRaw(const uint32_t& i) const {
            return ::tan(c_m.AtRaw(i));
        }

        /**
         * Used to determine if the resulting container is used in the expression.
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
    inline const VectorTan<C> tan(const atl::VectorExpression<typename C::RET_TYPE, C>& c) {
        return VectorTan<C>(c);
    }

    template<class C>
    struct VectorTanh : atl::VectorExpression<typename C::RET_TYPE, VectorTanh<C> > {
        const C& c_m;
        typedef typename C::RET_TYPE RET_TYPE;
        typedef typename C::BASE_TYPE BASE_TYPE;

        inline explicit VectorTanh(const atl::VectorExpression<typename C::RET_TYPE, C>& c) : c_m(c.Cast()) {

        }

        inline const size_t Size(const int32_t & dimension) const {
            switch (dimension) {
                case 0:
                    return c_m.Size(0);
                default:
                    return 0;
            }
        }

        inline const RET_TYPE operator()(size_t i) const {
#ifdef ENABLE_BOUNDS_CHECKING
            assert(i < c_m.Size(0));
#endif   

            return ::tanh(c_m(i));

        }

        /*
         *
         *Returns the first valid index for this container given dimension d. 
         */
        inline const size_t IndexMin(uint32_t d) const {
            return c_m.IndexMin(d);
        }

        /*
         *
         *Returns the last valid index for this container given dimension d. 
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
        inline const RET_TYPE AtRaw(const uint32_t& i) const {
            return ::tanh(c_m.AtRaw(i));
        }

        /**
         * Used to determine if the resulting container is used in the expression.
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
    inline const VectorTanh<C> tanh(const atl::VectorExpression<typename C::RET_TYPE, C>& c) {
        return VectorTanh<C>(c);
    }

    template< class LHS, class RHS>
    struct VectorPow : VectorExpression<typename atl::PromoteType<typename LHS::RET_TYPE, typename RHS::RET_TYPE >::return_type, VectorPow<LHS, RHS> > {
        typedef typename LHS::RET_TYPE RET_TYPEL;
        typedef typename RHS::RET_TYPE RET_TYPER;

        typedef typename atl::PromoteType<RET_TYPEL, RET_TYPER>::return_type RET_TYPE;
        typedef typename atl::PromoteType<typename LHS::BASE_TYPE, typename RHS::BASE_TYPE>::return_type BASE_TYPE;

        const LHS& lhs_m;
        const RHS& rhs_m;

        inline explicit VectorPow(const atl::VectorExpression<typename LHS::RET_TYPE, LHS>& lhs, const atl::VectorExpression<typename RHS::RET_TYPE, RHS>& rhs) : lhs_m(lhs.Cast()), rhs_m(rhs.Cast()) {



        }

        inline const size_t Size(const int32_t & dimension) const {
            return lhs_m.Size(dimension) < rhs_m.Size(dimension) ? lhs_m.Size(dimension) : rhs_m.Size(dimension);
        }

        inline const RET_TYPE operator()(const uint32_t& i) const {
            return ::pow(lhs_m(i), rhs_m(i));
        }

        /*
         *
         *Returns the first valid index for this container given dimension d. 
         */
        inline const size_t IndexMin(uint32_t d) const {
            return std::max(lhs_m.IndexMin(d), rhs_m.IndexMin(d));
        }

        /*
         *
         *Returns the last valid index for this container given dimension d. 
         */
        inline const size_t IndexMax(uint32_t d) const {
            return std::min(lhs_m.IndexMax(d), rhs_m.IndexMax(d));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i) const {
            return ::pow(lhs_m.AtRaw(i), rhs_m.AtRaw(i));
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
    struct VectorPowScalar : VectorExpression<typename PromoteType<typename LHS::RET_TYPE, T>::return_type, VectorPowScalar< LHS, T> > {
        const LHS& lhs_m;
        const T& rhs_m;

        typedef typename atl::PromoteType<typename LHS::RET_TYPE, T>::return_type RET_TYPE;
        typedef typename atl::PromoteType<typename LHS::BASE_TYPE, T>::return_type BASE_TYPE;

        inline explicit VectorPowScalar(const VectorExpression<typename LHS::RET_TYPE, LHS>& lhs, const T & rhs) : lhs_m(lhs.Cast()), rhs_m(rhs) {



        }

        inline const size_t Size(const int32_t & dimension) const {
            return lhs_m.Size(dimension);
        }

        inline const RET_TYPE operator()(const uint32_t& i) const {
            return ::pow(lhs_m(i), rhs_m);
        }

        /*
         *
         *Returns the first valid index for this container given dimension d. 
         */
        inline const size_t IndexMin(uint32_t d) const {
            return (lhs_m.IndexMin(d));
        }

        /*
         *
         *Returns the last valid index for this container given dimension d. 
         */
        inline const size_t IndexMax(uint32_t d) const {
            return (lhs_m.IndexMax(d));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i) const {
            return ::pow(lhs_m.AtRaw(i), rhs_m);
        }

        /**
         * Used to determine if the resulting matrix is used in the expression.
         * 
         * @param aliased
         * @param ptr
         */
        inline void IsAliased(bool& aliased, void* ptr) const {
            lhs_m.IsAliased(aliased, ptr);
        }

        inline void ExpressionLength(uint32_t& length)const {
            length++;
            lhs_m.ExpressionLength(length);
        }


    };

    template<class T, class RHS>
    struct ScalarVectorPow : VectorExpression<typename PromoteType<typename RHS::RET_TYPE, T>::return_type, ScalarVectorPow<T, RHS> > {
        typedef typename atl::PromoteType<T, typename RHS::RET_TYPE>::return_type RET_TYPE;
        typedef typename atl::PromoteType<T, typename RHS::BASE_TYPE>::return_type BASE_TYPE;

        const T& lhs_m;
        const RHS& rhs_m;

        inline explicit ScalarVectorPow(const T& lhs, const VectorExpression<typename RHS::RET_TYPE, RHS> & rhs) : lhs_m(lhs), rhs_m(rhs.Cast()) {



        }

        inline const size_t Size(const int32_t & dimension) const {
            return rhs_m.Size(dimension);
        }

        inline const RET_TYPE operator()(const uint32_t& i) const {
            return ::pow(lhs_m, rhs_m(i));
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

        inline const RET_TYPE AtRaw(const uint32_t& i) const {
            return ::pow(lhs_m, rhs_m.AtRaw(i));
        }

        /**
         * Used to determine if the resulting matrix is used in the expression.
         * 
         * @param aliased
         * @param ptr
         */
        inline void IsAliased(bool& aliased, void* ptr) const {
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
    inline const VectorPow< LHS, RHS> operator+(const MatrixExpression<typename LHS::RET_TYPE, LHS>& a,
            const MatrixExpression<typename RHS::RET_TYPE, RHS>& b) {
        return VectorPow< LHS, RHS > (a.Cast(), b.Cast());
    }



#define ATL_VECTOR_VectorPow_SCALAR(TYPE) \
    template< class LHS>      \
        inline const VectorPowScalar<LHS,TYPE> \
    pow(const VectorExpression<typename LHS::RET_TYPE, LHS>& a, const TYPE& b) {\
            return VectorPowScalar<LHS,TYPE > (a.Cast(), b);\
        } \
    

    ATL_VECTOR_VectorPow_SCALAR(short)
    ATL_VECTOR_VectorPow_SCALAR(unsigned short)
    ATL_VECTOR_VectorPow_SCALAR(int)
    ATL_VECTOR_VectorPow_SCALAR(unsigned int)
    ATL_VECTOR_VectorPow_SCALAR(long)
    ATL_VECTOR_VectorPow_SCALAR(unsigned long)
    ATL_VECTOR_VectorPow_SCALAR(float)
    ATL_VECTOR_VectorPow_SCALAR(double)
    ATL_VECTOR_VectorPow_SCALAR(long double)
    ATL_VECTOR_VectorPow_SCALAR(atl::Variable<float>)
    ATL_VECTOR_VectorPow_SCALAR(atl::Variable<double>)
    ATL_VECTOR_VectorPow_SCALAR(atl::Variable<long double>)


#define ATL_VectorPow_SCALAR_VECTOR(TYPE) \
    template< class RHS>      \
        inline const ScalarVectorPow<TYPE,RHS> \
    pow(const TYPE& a, const VectorExpression<typename RHS::RET_TYPE, RHS>& b ) {\
            return ScalarVectorPow<TYPE,RHS > (a, b.Cast());\
        } \
            
    ATL_VectorPow_SCALAR_VECTOR(short)
    ATL_VectorPow_SCALAR_VECTOR(unsigned short)
    ATL_VectorPow_SCALAR_VECTOR(int)
    ATL_VectorPow_SCALAR_VECTOR(unsigned int)
    ATL_VectorPow_SCALAR_VECTOR(long)
    ATL_VectorPow_SCALAR_VECTOR(unsigned long)
    ATL_VectorPow_SCALAR_VECTOR(float)
    ATL_VectorPow_SCALAR_VECTOR(double)
    ATL_VectorPow_SCALAR_VECTOR(long double)
    ATL_VectorPow_SCALAR_VECTOR(atl::Variable<float>)
    ATL_VectorPow_SCALAR_VECTOR(atl::Variable<double>)
    ATL_VectorPow_SCALAR_VECTOR(atl::Variable<long double>)



    template< class LHS, class RHS>
    struct VectorATan2 : VectorExpression<typename atl::PromoteType<typename LHS::RET_TYPE, typename RHS::RET_TYPE >::return_type, VectorATan2<LHS, RHS> > {
        typedef typename LHS::RET_TYPE RET_TYPEL;
        typedef typename RHS::RET_TYPE RET_TYPER;

        typedef typename atl::PromoteType<RET_TYPEL, RET_TYPER>::return_type RET_TYPE;
        typedef typename atl::PromoteType<typename LHS::BASE_TYPE, typename RHS::BASE_TYPE>::return_type BASE_TYPE;

        const LHS& lhs_m;
        const RHS& rhs_m;

        inline explicit VectorATan2(const atl::VectorExpression<typename LHS::RET_TYPE, LHS>& lhs, const atl::VectorExpression<typename RHS::RET_TYPE, RHS>& rhs) : lhs_m(lhs.Cast()), rhs_m(rhs.Cast()) {



        }

        inline const size_t Size(const int32_t & dimension) const {
            return lhs_m.Size(dimension) < rhs_m.Size(dimension) ? lhs_m.Size(dimension) : rhs_m.Size(dimension);
        }

        inline const RET_TYPE operator()(const uint32_t& i) const {
            return ::atan2(lhs_m(i), rhs_m(i));
        }

        /*
         *
         *Returns the first valid index for this container given dimension d. 
         */
        inline const size_t IndexMin(uint32_t d) const {
            return std::max(lhs_m.IndexMin(d), rhs_m.IndexMin(d));
        }

        /*
         *
         *Returns the last valid index for this container given dimension d. 
         */
        inline const size_t IndexMax(uint32_t d) const {
            return std::min(lhs_m.IndexMax(d), rhs_m.IndexMax(d));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i) const {
            return ::atan2(lhs_m.AtRaw(i), rhs_m.AtRaw(i));
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
    struct VectorATan2Scalar : VectorExpression<typename PromoteType<typename LHS::RET_TYPE, T>::return_type, VectorATan2Scalar< LHS, T> > {
        const LHS& lhs_m;
        const T& rhs_m;

        typedef typename atl::PromoteType<typename LHS::RET_TYPE, T>::return_type RET_TYPE;
        typedef typename atl::PromoteType<typename LHS::BASE_TYPE, T>::return_type BASE_TYPE;

        inline explicit VectorATan2Scalar(const VectorExpression<typename LHS::RET_TYPE, LHS>& lhs, const T & rhs) : lhs_m(lhs.Cast()), rhs_m(rhs) {



        }

        inline const size_t Size(const int32_t & dimension) const {
            return lhs_m.Size(dimension);
        }

        inline const RET_TYPE operator()(const uint32_t& i) const {
            return ::atan2(lhs_m(i), rhs_m);
        }

        /*
         *
         *Returns the first valid index for this container given dimension d. 
         */
        inline const size_t IndexMin(uint32_t d) const {
            return (lhs_m.IndexMin(d));
        }

        /*
         *
         *Returns the last valid index for this container given dimension d. 
         */
        inline const size_t IndexMax(uint32_t d) const {
            return (lhs_m.IndexMax(d));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i) const {
            return ::atan2(lhs_m.AtRaw(i), rhs_m);
        }

        /**
         * Used to determine if the resulting matrix is used in the expression.
         * 
         * @param aliased
         * @param ptr
         */
        inline void IsAliased(bool& aliased, void* ptr) const {
            lhs_m.IsAliased(aliased, ptr);
        }

        inline void ExpressionLength(uint32_t& length)const {
            length++;
            lhs_m.ExpressionLength(length);
        }


    };

    template<class T, class RHS>
    struct ScalarVectorATan2 : VectorExpression<typename PromoteType<typename RHS::RET_TYPE, T>::return_type, ScalarVectorATan2<T, RHS> > {
        typedef typename atl::PromoteType<T, typename RHS::RET_TYPE>::return_type RET_TYPE;
        typedef typename atl::PromoteType<T, typename RHS::BASE_TYPE>::return_type BASE_TYPE;

        const T& lhs_m;
        const RHS& rhs_m;

        inline explicit ScalarVectorATan2(const T& lhs, const VectorExpression<typename RHS::RET_TYPE, RHS> & rhs) : lhs_m(lhs), rhs_m(rhs.Cast()) {



        }

        inline const size_t Size(const int32_t & dimension) const {
            return rhs_m.Size(dimension);
        }

        inline const RET_TYPE operator()(const uint32_t& i) const {
            return ::atan2(lhs_m, rhs_m(i));
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

        inline const RET_TYPE AtRaw(const uint32_t& i) const {
            return ::atan2(lhs_m, rhs_m.AtRaw(i));
        }

        /**
         * Used to determine if the resulting matrix is used in the expression.
         * 
         * @param aliased
         * @param ptr
         */
        inline void IsAliased(bool& aliased, void* ptr) const {
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
    inline const VectorATan2< LHS, RHS> operator+(const MatrixExpression<typename LHS::RET_TYPE, LHS>& a,
            const MatrixExpression<typename RHS::RET_TYPE, RHS>& b) {
        return VectorATan2< LHS, RHS > (a.Cast(), b.Cast());
    }



#define ATL_VECTOR_VectorATan2_SCALAR(TYPE) \
    template< class LHS>      \
        inline const VectorATan2Scalar<LHS,TYPE> \
    atan2(const VectorExpression<typename LHS::RET_TYPE, LHS>& a, const TYPE& b) {\
            return VectorATan2Scalar<LHS,TYPE > (a.Cast(), b);\
        } \
    
    ATL_VECTOR_VectorATan2_SCALAR(short)
    ATL_VECTOR_VectorATan2_SCALAR(unsigned short)
    ATL_VECTOR_VectorATan2_SCALAR(int)
    ATL_VECTOR_VectorATan2_SCALAR(unsigned int)
    ATL_VECTOR_VectorATan2_SCALAR(long)
    ATL_VECTOR_VectorATan2_SCALAR(unsigned long)
    ATL_VECTOR_VectorATan2_SCALAR(float)
    ATL_VECTOR_VectorATan2_SCALAR(double)
    ATL_VECTOR_VectorATan2_SCALAR(long double)
    ATL_VECTOR_VectorATan2_SCALAR(atl::Variable<float>)
    ATL_VECTOR_VectorATan2_SCALAR(atl::Variable<double>)
    ATL_VECTOR_VectorATan2_SCALAR(atl::Variable<long double>)


#define ATL_VectorATan2_SCALAR_VECTOR(TYPE) \
    template< class RHS>      \
        inline const ScalarVectorATan2<TYPE,RHS> \
    atan2(const TYPE& a, const VectorExpression<typename RHS::RET_TYPE, RHS>& b ) {\
            return ScalarVectorATan2<TYPE,RHS > (a, b.Cast());\
        } \
            
    ATL_VectorATan2_SCALAR_VECTOR(short)
    ATL_VectorATan2_SCALAR_VECTOR(unsigned short)
    ATL_VectorATan2_SCALAR_VECTOR(int)
    ATL_VectorATan2_SCALAR_VECTOR(unsigned int)
    ATL_VectorATan2_SCALAR_VECTOR(long)
    ATL_VectorATan2_SCALAR_VECTOR(unsigned long)
    ATL_VectorATan2_SCALAR_VECTOR(float)
    ATL_VectorATan2_SCALAR_VECTOR(double)
    ATL_VectorATan2_SCALAR_VECTOR(long double)
    ATL_VectorATan2_SCALAR_VECTOR(atl::Variable<float>)
    ATL_VectorATan2_SCALAR_VECTOR(atl::Variable<double>)
    ATL_VectorATan2_SCALAR_VECTOR(atl::Variable<long double>)





}
#endif