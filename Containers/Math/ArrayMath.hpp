#ifndef ARRAYMATH_HPP
#define ARRAYMATH_HPP

#ifdef ENABLE_BOUNDS_CHECKING
#include <assert.h>
#endif




namespace atl {

    template<class C>
    struct ArrayACos : atl::ArrayExpression<typename C::RET_TYPE, ArrayACos<C> > {
        const C& c_m;
        typedef typename C::RET_TYPE RET_TYPE;
        typedef typename C::BASE_TYPE BASE_TYPE;

        inline explicit ArrayACos(const atl::ArrayExpression<typename C::RET_TYPE, C>& c) : c_m(c.Cast()) {

        }

        inline const size_t Dimensions() const {
            return c_m.Dimensions();
        }

        inline const size_t Size(const int32_t & d) const {
            return c_m.Size(d);
        }

        inline const RET_TYPE operator()(const uint32_t & i) const {
            return ::acos(c_m(i));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j) const {
            return ::acos(c_m(i, j));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k) const {
            return ::acos(c_m(i, j, k));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l)const {
            return ::acos(c_m(i, j, k, l));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m) const {
            return ::acos(c_m(i, j, k, l, m));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n) const {
            return ::acos(c_m(i, j, k, l, m, n));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n, const uint32_t & o) const {
            return ::acos(c_m(i, j, k, l, m, n, o));
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

        inline const RET_TYPE AtRaw(const uint32_t & i) const {
            return ::acos(c_m.AtRaw(i));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j) const {
            return ::acos(c_m.AtRaw(i, j));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k) const {
            return ::acos(c_m.AtRaw(i, j, k));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l)const {
            return ::acos(c_m.AtRaw(i, j, k, l));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m) const {
            return ::acos(c_m.AtRaw(i, j, k, l, m));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n) const {
            return ::acos(c_m.AtRaw(i, j, k, l, m, n));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n, const uint32_t & o) const {
            return ::acos(c_m.AtRaw(i, j, k, l, m, n, o));
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
    inline const ArrayACos<C> acos(const atl::ArrayExpression<typename C::RET_TYPE, C>& c) {
        return ArrayACos<C>(c);
    }

    template<class C>
    struct ArrayASin : atl::ArrayExpression<typename C::RET_TYPE, ArrayASin<C> > {
        const C& c_m;
        typedef typename C::RET_TYPE RET_TYPE;
        typedef typename C::BASE_TYPE BASE_TYPE;

        inline explicit ArrayASin(const atl::ArrayExpression<typename C::RET_TYPE, C>& c) : c_m(c.Cast()) {

        }

        inline const size_t Dimensions() const {
            return c_m.Dimensions();
        }

        inline const size_t Size(const int32_t & d) const {
            return c_m.Size(d);
        }

        inline const RET_TYPE operator()(const uint32_t & i) const {
            return ::asin(c_m(i));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j) const {
            return ::asin(c_m(i, j));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k) const {
            return ::asin(c_m(i, j, k));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l)const {
            return ::asin(c_m(i, j, k, l));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m) const {
            return ::asin(c_m(i, j, k, l, m));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n) const {
            return ::asin(c_m(i, j, k, l, m, n));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n, const uint32_t & o) const {
            return ::asin(c_m(i, j, k, l, m, n, o));
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

        inline const RET_TYPE AtRaw(const uint32_t & i) const {
            return ::asin(c_m.AtRaw(i));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j) const {
            return ::asin(c_m.AtRaw(i, j));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k) const {
            return ::asin(c_m.AtRaw(i, j, k));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l)const {
            return ::asin(c_m.AtRaw(i, j, k, l));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m) const {
            return ::asin(c_m.AtRaw(i, j, k, l, m));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n) const {
            return ::asin(c_m.AtRaw(i, j, k, l, m, n));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n, const uint32_t & o) const {
            return ::asin(c_m.AtRaw(i, j, k, l, m, n, o));
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
    inline const ArrayASin<C> asin(const atl::ArrayExpression<typename C::RET_TYPE, C>& c) {
        return ArrayASin<C>(c);
    }

    template<class C>
    struct ArrayATan : atl::ArrayExpression<typename C::RET_TYPE, ArrayATan<C> > {
        const C& c_m;
        typedef typename C::RET_TYPE RET_TYPE;
        typedef typename C::BASE_TYPE BASE_TYPE;

        inline explicit ArrayATan(const atl::ArrayExpression<typename C::RET_TYPE, C>& c) : c_m(c.Cast()) {

        }

        inline const size_t Dimensions() const {
            return c_m.Dimensions();
        }

        inline const size_t Size(const int32_t & d) const {
            return c_m.Size(d);
        }

        inline const RET_TYPE operator()(const uint32_t & i) const {
            return ::atan(c_m(i));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j) const {
            return ::atan(c_m(i, j));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k) const {
            return ::atan(c_m(i, j, k));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l)const {
            return ::atan(c_m(i, j, k, l));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m) const {
            return ::atan(c_m(i, j, k, l, m));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n) const {
            return ::atan(c_m(i, j, k, l, m, n));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n, const uint32_t & o) const {
            return ::atan(c_m(i, j, k, l, m, n, o));
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

        inline const RET_TYPE AtRaw(const uint32_t & i) const {
            return ::atan(c_m.AtRaw(i));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j) const {
            return ::atan(c_m.AtRaw(i, j));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k) const {
            return ::atan(c_m.AtRaw(i, j, k));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l)const {
            return ::atan(c_m.AtRaw(i, j, k, l));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m) const {
            return ::atan(c_m.AtRaw(i, j, k, l, m));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n) const {
            return ::atan(c_m.AtRaw(i, j, k, l, m, n));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n, const uint32_t & o) const {
            return ::atan(c_m.AtRaw(i, j, k, l, m, n, o));
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
    inline const ArrayATan<C> atan(const atl::ArrayExpression<typename C::RET_TYPE, C>& c) {
        return ArrayATan<C>(c);
    }

    template<class C>
    struct ArrayACeil : atl::ArrayExpression<typename C::RET_TYPE, ArrayACeil<C> > {
        const C& c_m;
        typedef typename C::RET_TYPE RET_TYPE;
        typedef typename C::BASE_TYPE BASE_TYPE;

        inline explicit ArrayACeil(const atl::ArrayExpression<typename C::RET_TYPE, C>& c) : c_m(c.Cast()) {

        }

        inline const size_t Dimensions() const {
            return c_m.Dimensions();
        }

        inline const size_t Size(const int32_t & d) const {
            return c_m.Size(d);
        }

        inline const RET_TYPE operator()(const uint32_t & i) const {
            return ::ceil(c_m(i));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j) const {
            return ::ceil(c_m(i, j));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k) const {
            return ::ceil(c_m(i, j, k));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l)const {
            return ::ceil(c_m(i, j, k, l));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m) const {
            return ::ceil(c_m(i, j, k, l, m));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n) const {
            return ::ceil(c_m(i, j, k, l, m, n));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n, const uint32_t & o) const {
            return ::ceil(c_m(i, j, k, l, m, n, o));
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

        inline const RET_TYPE AtRaw(const uint32_t & i) const {
            return ::ceil(c_m.AtRaw(i));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j) const {
            return ::ceil(c_m.AtRaw(i, j));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k) const {
            return ::ceil(c_m.AtRaw(i, j, k));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l)const {
            return ::ceil(c_m.AtRaw(i, j, k, l));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m) const {
            return ::ceil(c_m.AtRaw(i, j, k, l, m));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n) const {
            return ::ceil(c_m.AtRaw(i, j, k, l, m, n));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n, const uint32_t & o) const {
            return ::ceil(c_m.AtRaw(i, j, k, l, m, n, o));
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
    inline const ArrayACeil<C> ceil(const atl::ArrayExpression<typename C::RET_TYPE, C>& c) {
        return ArrayACeil<C>(c);
    }

    template<class C>
    struct ArrayCos : atl::ArrayExpression<typename C::RET_TYPE, ArrayCos<C> > {
        const C& c_m;
        typedef typename C::RET_TYPE RET_TYPE;
        typedef typename C::BASE_TYPE BASE_TYPE;

        inline explicit ArrayCos(const atl::ArrayExpression<typename C::RET_TYPE, C>& c) : c_m(c.Cast()) {

        }

        inline const size_t Dimensions() const {
            return c_m.Dimensions();
        }

        inline const size_t Size(const int32_t & d) const {
            return c_m.Size(d);
        }

        inline const RET_TYPE operator()(const uint32_t & i) const {
            return ::cos(c_m(i));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j) const {
            return ::cos(c_m(i, j));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k) const {
            return ::cos(c_m(i, j, k));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l)const {
            return ::cos(c_m(i, j, k, l));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m) const {
            return ::cos(c_m(i, j, k, l, m));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n) const {
            return ::cos(c_m(i, j, k, l, m, n));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n, const uint32_t & o) const {
            return ::cos(c_m(i, j, k, l, m, n, o));
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

        inline const RET_TYPE AtRaw(const uint32_t & i) const {
            return ::cos(c_m.AtRaw(i));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j) const {
            return ::cos(c_m.AtRaw(i, j));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k) const {
            return ::cos(c_m.AtRaw(i, j, k));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l)const {
            return ::cos(c_m.AtRaw(i, j, k, l));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m) const {
            return ::cos(c_m.AtRaw(i, j, k, l, m));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n) const {
            return ::cos(c_m.AtRaw(i, j, k, l, m, n));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n, const uint32_t & o) const {
            return ::cos(c_m.AtRaw(i, j, k, l, m, n, o));
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
    inline const ArrayCos<C> cos(const atl::ArrayExpression<typename C::RET_TYPE, C>& c) {
        return ArrayCos<C>(c);
    }

    template<class C>
    struct ArrayCosh : atl::ArrayExpression<typename C::RET_TYPE, ArrayCosh<C> > {
        const C& c_m;
        typedef typename C::RET_TYPE RET_TYPE;
        typedef typename C::BASE_TYPE BASE_TYPE;

        inline explicit ArrayCosh(const atl::ArrayExpression<typename C::RET_TYPE, C>& c) : c_m(c.Cast()) {

        }

        inline const size_t Dimensions() const {
            return c_m.Dimensions();
        }

        inline const size_t Size(const int32_t & d) const {
            return c_m.Size(d);
        }

        inline const RET_TYPE operator()(const uint32_t & i) const {
            return ::cosh(c_m(i));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j) const {
            return ::cosh(c_m(i, j));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k) const {
            return ::cosh(c_m(i, j, k));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l)const {
            return ::cosh(c_m(i, j, k, l));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m) const {
            return ::cosh(c_m(i, j, k, l, m));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n) const {
            return ::cosh(c_m(i, j, k, l, m, n));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n, const uint32_t & o) const {
            return ::cosh(c_m(i, j, k, l, m, n, o));
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

        inline const RET_TYPE AtRaw(const uint32_t & i) const {
            return ::cosh(c_m.AtRaw(i));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j) const {
            return ::cosh(c_m.AtRaw(i, j));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k) const {
            return ::cosh(c_m.AtRaw(i, j, k));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l)const {
            return ::cosh(c_m.AtRaw(i, j, k, l));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m) const {
            return ::cosh(c_m.AtRaw(i, j, k, l, m));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n) const {
            return ::cosh(c_m.AtRaw(i, j, k, l, m, n));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n, const uint32_t & o) const {
            return ::cosh(c_m.AtRaw(i, j, k, l, m, n, o));
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
    inline const ArrayCosh<C> cosh(const atl::ArrayExpression<typename C::RET_TYPE, C>& c) {
        return ArrayCosh<C>(c);
    }

    template<class C>
    struct ArrayExp : atl::ArrayExpression<typename C::RET_TYPE, ArrayExp<C> > {
        const C& c_m;
        typedef typename C::RET_TYPE RET_TYPE;
        typedef typename C::BASE_TYPE BASE_TYPE;

        inline explicit ArrayExp(const atl::ArrayExpression<typename C::RET_TYPE, C>& c) : c_m(c.Cast()) {

        }

        inline const size_t Dimensions() const {
            return c_m.Dimensions();
        }

        inline const size_t Size(const int32_t & d) const {
            return c_m.Size(d);
        }

        inline const RET_TYPE operator()(const uint32_t & i) const {
            return ::exp(c_m(i));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j) const {
            return ::exp(c_m(i, j));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k) const {
            return ::exp(c_m(i, j, k));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l)const {
            return ::exp(c_m(i, j, k, l));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m) const {
            return ::exp(c_m(i, j, k, l, m));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n) const {
            return ::exp(c_m(i, j, k, l, m, n));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n, const uint32_t & o) const {
            return ::exp(c_m(i, j, k, l, m, n, o));
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

        inline const RET_TYPE AtRaw(const uint32_t & i) const {
            return ::exp(c_m.AtRaw(i));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j) const {
            return ::exp(c_m.AtRaw(i, j));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k) const {
            return ::exp(c_m.AtRaw(i, j, k));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l)const {
            return ::exp(c_m.AtRaw(i, j, k, l));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m) const {
            return ::exp(c_m.AtRaw(i, j, k, l, m));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n) const {
            return ::exp(c_m.AtRaw(i, j, k, l, m, n));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n, const uint32_t & o) const {
            return ::exp(c_m.AtRaw(i, j, k, l, m, n, o));
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
    inline const ArrayExp<C> exp(const atl::ArrayExpression<typename C::RET_TYPE, C>& c) {
        return ArrayExp<C>(c);
    }

    template<class C>
    struct ArrayFabs : atl::ArrayExpression<typename C::RET_TYPE, ArrayFabs<C> > {
        const C& c_m;
        typedef typename C::RET_TYPE RET_TYPE;
        typedef typename C::BASE_TYPE BASE_TYPE;

        inline explicit ArrayFabs(const atl::ArrayExpression<typename C::RET_TYPE, C>& c) : c_m(c.Cast()) {

        }

        inline const size_t Dimensions() const {
            return c_m.Dimensions();
        }

        inline const size_t Size(const int32_t & d) const {
            return c_m.Size(d);
        }

        inline const RET_TYPE operator()(const uint32_t & i) const {
            return ::fabs(c_m(i));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j) const {
            return ::fabs(c_m(i, j));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k) const {
            return ::fabs(c_m(i, j, k));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l)const {
            return ::fabs(c_m(i, j, k, l));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m) const {
            return ::fabs(c_m(i, j, k, l, m));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n) const {
            return ::fabs(c_m(i, j, k, l, m, n));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n, const uint32_t & o) const {
            return ::fabs(c_m(i, j, k, l, m, n, o));
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

        inline const RET_TYPE AtRaw(const uint32_t & i) const {
            return ::fabs(c_m.AtRaw(i));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j) const {
            return ::fabs(c_m.AtRaw(i, j));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k) const {
            return ::fabs(c_m.AtRaw(i, j, k));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l)const {
            return ::fabs(c_m.AtRaw(i, j, k, l));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m) const {
            return ::fabs(c_m.AtRaw(i, j, k, l, m));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n) const {
            return ::fabs(c_m.AtRaw(i, j, k, l, m, n));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n, const uint32_t & o) const {
            return ::fabs(c_m.AtRaw(i, j, k, l, m, n, o));
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
    inline const ArrayFabs<C> fabs(const atl::ArrayExpression<typename C::RET_TYPE, C>& c) {
        return ArrayFabs<C>(c);
    }

    template<class C>
    struct ArrayFloor : atl::ArrayExpression<typename C::RET_TYPE, ArrayFloor<C> > {
        const C& c_m;
        typedef typename C::RET_TYPE RET_TYPE;
        typedef typename C::BASE_TYPE BASE_TYPE;

        inline explicit ArrayFloor(const atl::ArrayExpression<typename C::RET_TYPE, C>& c) : c_m(c.Cast()) {

        }

        inline const size_t Dimensions() const {
            return c_m.Dimensions();
        }

        inline const size_t Size(const int32_t & d) const {
            return c_m.Size(d);
        }

        inline const RET_TYPE operator()(const uint32_t & i) const {
            return ::floor(c_m(i));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j) const {
            return ::floor(c_m(i, j));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k) const {
            return ::floor(c_m(i, j, k));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l)const {
            return ::floor(c_m(i, j, k, l));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m) const {
            return ::floor(c_m(i, j, k, l, m));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n) const {
            return ::floor(c_m(i, j, k, l, m, n));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n, const uint32_t & o) const {
            return ::floor(c_m(i, j, k, l, m, n, o));
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

        inline const RET_TYPE AtRaw(const uint32_t & i) const {
            return ::floor(c_m.AtRaw(i));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j) const {
            return ::floor(c_m.AtRaw(i, j));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k) const {
            return ::floor(c_m.AtRaw(i, j, k));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l)const {
            return ::floor(c_m.AtRaw(i, j, k, l));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m) const {
            return ::floor(c_m.AtRaw(i, j, k, l, m));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n) const {
            return ::floor(c_m.AtRaw(i, j, k, l, m, n));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n, const uint32_t & o) const {
            return ::floor(c_m.AtRaw(i, j, k, l, m, n, o));
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
    inline const ArrayFloor<C> floor(const atl::ArrayExpression<typename C::RET_TYPE, C>& c) {
        return ArrayFloor<C>(c);
    }

    template<class C>
    struct ArrayLog : atl::ArrayExpression<typename C::RET_TYPE, ArrayLog<C> > {
        const C& c_m;
        typedef typename C::RET_TYPE RET_TYPE;
        typedef typename C::BASE_TYPE BASE_TYPE;

        inline explicit ArrayLog(const atl::ArrayExpression<typename C::RET_TYPE, C>& c) : c_m(c.Cast()) {

        }

        inline const size_t Dimensions() const {
            return c_m.Dimensions();
        }

        inline const size_t Size(const int32_t & d) const {
            return c_m.Size(d);
        }

        inline const RET_TYPE operator()(const uint32_t & i) const {
            return ::log(c_m(i));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j) const {
            return ::log(c_m(i, j));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k) const {
            return ::log(c_m(i, j, k));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l)const {
            return ::log(c_m(i, j, k, l));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m) const {
            return ::log(c_m(i, j, k, l, m));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n) const {
            return ::log(c_m(i, j, k, l, m, n));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n, const uint32_t & o) const {
            return ::log(c_m(i, j, k, l, m, n, o));
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

        inline const RET_TYPE AtRaw(const uint32_t & i) const {
            return ::log(c_m.AtRaw(i));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j) const {
            return ::log(c_m.AtRaw(i, j));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k) const {
            return ::log(c_m.AtRaw(i, j, k));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l)const {
            return ::log(c_m.AtRaw(i, j, k, l));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m) const {
            return ::log(c_m.AtRaw(i, j, k, l, m));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n) const {
            return ::log(c_m.AtRaw(i, j, k, l, m, n));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n, const uint32_t & o) const {
            return ::log(c_m.AtRaw(i, j, k, l, m, n, o));
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
    inline const ArrayLog<C> log(const atl::ArrayExpression<typename C::RET_TYPE, C>& c) {
        return ArrayLog<C>(c);
    }

    template<class C>
    struct ArrayLog10 : atl::ArrayExpression<typename C::RET_TYPE, ArrayLog10<C> > {
        const C& c_m;
        typedef typename C::RET_TYPE RET_TYPE;
        typedef typename C::BASE_TYPE BASE_TYPE;

        inline explicit ArrayLog10(const atl::ArrayExpression<typename C::RET_TYPE, C>& c) : c_m(c.Cast()) {

        }

        inline const size_t Dimensions() const {
            return c_m.Dimensions();
        }

        inline const size_t Size(const int32_t & d) const {
            return c_m.Size(d);
        }

        inline const RET_TYPE operator()(const uint32_t & i) const {
            return ::log10(c_m(i));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j) const {
            return ::log10(c_m(i, j));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k) const {
            return ::log10(c_m(i, j, k));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l)const {
            return ::log10(c_m(i, j, k, l));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m) const {
            return ::log10(c_m(i, j, k, l, m));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n) const {
            return ::log10(c_m(i, j, k, l, m, n));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n, const uint32_t & o) const {
            return ::log10(c_m(i, j, k, l, m, n, o));
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

        inline const RET_TYPE AtRaw(const uint32_t & i) const {
            return ::log10(c_m.AtRaw(i));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j) const {
            return ::log10(c_m.AtRaw(i, j));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k) const {
            return ::log10(c_m.AtRaw(i, j, k));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l)const {
            return ::log10(c_m.AtRaw(i, j, k, l));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m) const {
            return ::log10(c_m.AtRaw(i, j, k, l, m));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n) const {
            return ::log10(c_m.AtRaw(i, j, k, l, m, n));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n, const uint32_t & o) const {
            return ::log10(c_m.AtRaw(i, j, k, l, m, n, o));
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
    inline const ArrayLog10<C> log10(const atl::ArrayExpression<typename C::RET_TYPE, C>& c) {
        return ArrayLog10<C>(c);
    }

    template<class C>
    struct ArraySin : atl::ArrayExpression<typename C::RET_TYPE, ArraySin<C> > {
        const C& c_m;
        typedef typename C::RET_TYPE RET_TYPE;
        typedef typename C::BASE_TYPE BASE_TYPE;

        inline explicit ArraySin(const atl::ArrayExpression<typename C::RET_TYPE, C>& c) : c_m(c.Cast()) {

        }

        inline const size_t Dimensions() const {
            return c_m.Dimensions();
        }

        inline const size_t Size(const int32_t & d) const {
            return c_m.Size(d);
        }

        inline const RET_TYPE operator()(const uint32_t & i) const {
            return ::sin(c_m(i));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j) const {
            return ::sin(c_m(i, j));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k) const {
            return ::sin(c_m(i, j, k));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l)const {
            return ::sin(c_m(i, j, k, l));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m) const {
            return ::sin(c_m(i, j, k, l, m));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n) const {
            return ::sin(c_m(i, j, k, l, m, n));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n, const uint32_t & o) const {
            return ::sin(c_m(i, j, k, l, m, n, o));
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

        inline const RET_TYPE AtRaw(const uint32_t & i) const {
            return ::sin(c_m.AtRaw(i));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j) const {
            return ::sin(c_m.AtRaw(i, j));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k) const {
            return ::sin(c_m.AtRaw(i, j, k));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l)const {
            return ::sin(c_m.AtRaw(i, j, k, l));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m) const {
            return ::sin(c_m.AtRaw(i, j, k, l, m));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n) const {
            return ::sin(c_m.AtRaw(i, j, k, l, m, n));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n, const uint32_t & o) const {
            return ::sin(c_m.AtRaw(i, j, k, l, m, n, o));
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
    inline const ArraySin<C> sin(const atl::ArrayExpression<typename C::RET_TYPE, C>& c) {
        return ArraySin<C>(c);
    }

    template<class C>
    struct ArraySinh : atl::ArrayExpression<typename C::RET_TYPE, ArraySinh<C> > {
        const C& c_m;
        typedef typename C::RET_TYPE RET_TYPE;
        typedef typename C::BASE_TYPE BASE_TYPE;

        inline explicit ArraySinh(const atl::ArrayExpression<typename C::RET_TYPE, C>& c) : c_m(c.Cast()) {

        }

        inline const size_t Dimensions() const {
            return c_m.Dimensions();
        }

        inline const size_t Size(const int32_t & d) const {
            return c_m.Size(d);
        }

        inline const RET_TYPE operator()(const uint32_t & i) const {
            return ::sinh(c_m(i));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j) const {
            return ::sinh(c_m(i, j));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k) const {
            return ::sinh(c_m(i, j, k));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l)const {
            return ::sinh(c_m(i, j, k, l));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m) const {
            return ::sinh(c_m(i, j, k, l, m));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n) const {
            return ::sinh(c_m(i, j, k, l, m, n));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n, const uint32_t & o) const {
            return ::sinh(c_m(i, j, k, l, m, n, o));
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

        inline const RET_TYPE AtRaw(const uint32_t & i) const {
            return ::sinh(c_m.AtRaw(i));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j) const {
            return ::sinh(c_m.AtRaw(i, j));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k) const {
            return ::sinh(c_m.AtRaw(i, j, k));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l)const {
            return ::sinh(c_m.AtRaw(i, j, k, l));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m) const {
            return ::sinh(c_m.AtRaw(i, j, k, l, m));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n) const {
            return ::sinh(c_m.AtRaw(i, j, k, l, m, n));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n, const uint32_t & o) const {
            return ::sinh(c_m.AtRaw(i, j, k, l, m, n, o));
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
    inline const ArraySinh<C> sinh(const atl::ArrayExpression<typename C::RET_TYPE, C>& c) {
        return ArraySinh<C>(c);
    }

    template<class C>
    struct ArraySqrt : atl::ArrayExpression<typename C::RET_TYPE, ArraySqrt<C> > {
        const C& c_m;
        typedef typename C::RET_TYPE RET_TYPE;
        typedef typename C::BASE_TYPE BASE_TYPE;

        inline explicit ArraySqrt(const atl::ArrayExpression<typename C::RET_TYPE, C>& c) : c_m(c.Cast()) {

        }

        inline const size_t Dimensions() const {
            return c_m.Dimensions();
        }

        inline const size_t Size(const int32_t & d) const {
            return c_m.Size(d);
        }

        inline const RET_TYPE operator()(const uint32_t & i) const {
            return ::sqrt(c_m(i));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j) const {
            return ::sqrt(c_m(i, j));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k) const {
            return ::sqrt(c_m(i, j, k));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l)const {
            return ::sqrt(c_m(i, j, k, l));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m) const {
            return ::sqrt(c_m(i, j, k, l, m));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n) const {
            return ::sqrt(c_m(i, j, k, l, m, n));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n, const uint32_t & o) const {
            return ::sqrt(c_m(i, j, k, l, m, n, o));
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

        inline const RET_TYPE AtRaw(const uint32_t & i) const {
            return ::sqrt(c_m.AtRaw(i));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j) const {
            return ::sqrt(c_m.AtRaw(i, j));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k) const {
            return ::sqrt(c_m.AtRaw(i, j, k));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l)const {
            return ::sqrt(c_m.AtRaw(i, j, k, l));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m) const {
            return ::sqrt(c_m.AtRaw(i, j, k, l, m));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n) const {
            return ::sqrt(c_m.AtRaw(i, j, k, l, m, n));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n, const uint32_t & o) const {
            return ::sqrt(c_m.AtRaw(i, j, k, l, m, n, o));
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
    inline const ArraySqrt<C> sqrt(const atl::ArrayExpression<typename C::RET_TYPE, C>& c) {
        return ArraySqrt<C>(c);
    }

    template<class C>
    struct ArrayTan : atl::ArrayExpression<typename C::RET_TYPE, ArrayTan<C> > {
        const C& c_m;
        typedef typename C::RET_TYPE RET_TYPE;
        typedef typename C::BASE_TYPE BASE_TYPE;

        inline explicit ArrayTan(const atl::ArrayExpression<typename C::RET_TYPE, C>& c) : c_m(c.Cast()) {

        }

        inline const size_t Dimensions() const {
            return c_m.Dimensions();
        }

        inline const size_t Size(const int32_t & d) const {
            return c_m.Size(d);
        }

        inline const RET_TYPE operator()(const uint32_t & i) const {
            return ::tan(c_m(i));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j) const {
            return ::tan(c_m(i, j));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k) const {
            return ::tan(c_m(i, j, k));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l)const {
            return ::tan(c_m(i, j, k, l));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m) const {
            return ::tan(c_m(i, j, k, l, m));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n) const {
            return ::tan(c_m(i, j, k, l, m, n));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n, const uint32_t & o) const {
            return ::tan(c_m(i, j, k, l, m, n, o));
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

        inline const RET_TYPE AtRaw(const uint32_t & i) const {
            return ::tan(c_m.AtRaw(i));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j) const {
            return ::tan(c_m.AtRaw(i, j));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k) const {
            return ::tan(c_m.AtRaw(i, j, k));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l)const {
            return ::tan(c_m.AtRaw(i, j, k, l));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m) const {
            return ::tan(c_m.AtRaw(i, j, k, l, m));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n) const {
            return ::tan(c_m.AtRaw(i, j, k, l, m, n));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n, const uint32_t & o) const {
            return ::tan(c_m.AtRaw(i, j, k, l, m, n, o));
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
    inline const ArrayTan<C> tan(const atl::ArrayExpression<typename C::RET_TYPE, C>& c) {
        return ArrayTan<C>(c);
    }

    template<class C>
    struct ArrayTanh : atl::ArrayExpression<typename C::RET_TYPE, ArrayTanh<C> > {
        const C& c_m;
        typedef typename C::RET_TYPE RET_TYPE;
        typedef typename C::BASE_TYPE BASE_TYPE;

        inline explicit ArrayTanh(const atl::ArrayExpression<typename C::RET_TYPE, C>& c) : c_m(c.Cast()) {

        }

        inline const size_t Dimensions() const {
            return c_m.Dimensions();
        }

        inline const size_t Size(const int32_t & d) const {
            return c_m.Size(d);
        }

        inline const RET_TYPE operator()(const uint32_t & i) const {
            return ::tanh(c_m(i));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j) const {
            return ::tanh(c_m(i, j));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k) const {
            return ::tanh(c_m(i, j, k));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l)const {
            return ::tanh(c_m(i, j, k, l));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m) const {
            return ::tanh(c_m(i, j, k, l, m));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n) const {
            return ::tanh(c_m(i, j, k, l, m, n));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n, const uint32_t & o) const {
            return ::tanh(c_m(i, j, k, l, m, n, o));
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

        inline const RET_TYPE AtRaw(const uint32_t & i) const {
            return ::tanh(c_m.AtRaw(i));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j) const {
            return ::tanh(c_m.AtRaw(i, j));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k) const {
            return ::tanh(c_m.AtRaw(i, j, k));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l)const {
            return ::tanh(c_m.AtRaw(i, j, k, l));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m) const {
            return ::tanh(c_m.AtRaw(i, j, k, l, m));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n) const {
            return ::tanh(c_m.AtRaw(i, j, k, l, m, n));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n, const uint32_t & o) const {
            return ::tanh(c_m.AtRaw(i, j, k, l, m, n, o));
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
    inline const ArrayTanh<C> tanh(const atl::ArrayExpression<typename C::RET_TYPE, C>& c) {
        return ArrayTanh<C>(c);
    }

    template< class LHS, class RHS>
    struct ArrayPow : ArrayExpression<typename atl::PromoteType<typename LHS::RET_TYPE, typename RHS::RET_TYPE >::return_type, ArrayPow<LHS, RHS> > {
        typedef typename LHS::RET_TYPE RET_TYPEL;
        typedef typename RHS::RET_TYPE RET_TYPER;

        typedef typename atl::PromoteType<RET_TYPEL, RET_TYPER>::return_type RET_TYPE;
        typedef typename atl::PromoteType<typename LHS::BASE_TYPE, typename RHS::BASE_TYPE>::return_type BASE_TYPE;

        const LHS& lhs_m;
        const RHS& rhs_m;

        inline explicit ArrayPow(const atl::ArrayExpression<typename LHS::RET_TYPE, LHS>& lhs, const atl::ArrayExpression<typename RHS::RET_TYPE, RHS>& rhs) : lhs_m(lhs.Cast()), rhs_m(rhs.Cast()) {



        }

        inline const size_t Size(const int32_t & dimension) const {
            return lhs_m.Size(dimension) < rhs_m.Size(dimension) ? lhs_m.Size(dimension) : rhs_m.Size(dimension);
        }

        inline const size_t Dimensions() const {
            return lhs_m.Dimensions() < rhs_m.Dimensions() ? lhs_m.Dimensions() : rhs_m.Dimensions();
        }

        inline const RET_TYPE operator()(const uint32_t & i) const {
            return ::pow(lhs_m(i), rhs_m(i));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j) const {
            return ::pow(lhs_m(i, j), rhs_m(i, j));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k) const {
            return ::pow(lhs_m(i, j, k), rhs_m(i, j, k));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l)const {
            return ::pow(lhs_m(i, j, k, l), rhs_m(i, j, k, l));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m) const {
            return ::pow(lhs_m(i, j, k, l, m), rhs_m(i, j, k, l, m));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n) const {
            return ::pow(lhs_m(i, j, k, l, m, n), rhs_m(i, j, k, l, m, n));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n, const uint32_t & o) const {
            return ::pow(lhs_m(i, j, k, l, m, n, o), rhs_m(i, j, k, l, m, n, o));
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

        inline const RET_TYPE AtRaw(const uint32_t & i) const {
            return ::pow(lhs_m.AtRaw(i), rhs_m.AtRaw(i));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j) const {
            return ::pow(lhs_m.AtRaw(i, j), rhs_m.AtRaw(i, j));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k) const {
            return ::pow(lhs_m.AtRaw(i, j, k), rhs_m.AtRaw(i, j, k));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l)const {
            return ::pow(lhs_m.AtRaw(i, j, k, l), rhs_m.AtRaw(i, j, k, l));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m) const {
            return ::pow(lhs_m.AtRaw(i, j, k, l, m), rhs_m.AtRaw(i, j, k, l, m));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n) const {
            return ::pow(lhs_m.AtRaw(i, j, k, l, m, n), rhs_m.AtRaw(i, j, k, l, m, n));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n, const uint32_t & o) const {
            return ::pow(lhs_m.AtRaw(i, j, k, l, m, n, o), rhs_m.AtRaw(i, j, k, l, m, n, o));
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
    struct ArrayPowScalar : ArrayExpression<typename PromoteType<typename LHS::RET_TYPE, T>::return_type, ArrayPowScalar< LHS, T> > {
        const LHS& lhs_m;
        const T& rhs_m;

        typedef typename atl::PromoteType<typename LHS::RET_TYPE, T>::return_type RET_TYPE;
        typedef typename atl::PromoteType<typename LHS::BASE_TYPE, T>::return_type BASE_TYPE;

        inline explicit ArrayPowScalar(const ArrayExpression<typename LHS::RET_TYPE, LHS>& lhs, const T & rhs) : lhs_m(lhs.Cast()), rhs_m(rhs) {



        }

        inline const size_t Size(const int32_t & dimension) const {
            return lhs_m.Size(dimension);
        }

        inline const size_t Dimensions() const {
            return lhs_m.Dimensions();
        }

        inline const RET_TYPE operator()(const uint32_t & i) const {
            return ::pow(lhs_m(i), rhs_m);
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j) const {
            return ::pow(lhs_m(i, j), rhs_m);
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k) const {
            return ::pow(lhs_m(i, j, k), rhs_m);
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l)const {
            return ::pow(lhs_m(i, j, k, l), rhs_m);
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m) const {
            return ::pow(lhs_m(i, j, k, l, m), rhs_m);
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n) const {
            return ::pow(lhs_m(i, j, k, l, m, n), rhs_m);
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n, const uint32_t & o) const {
            return ::pow(lhs_m(i, j, k, l, m, n, o), rhs_m);
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

        inline const RET_TYPE AtRaw(const uint32_t & i) const {
            return ::pow(lhs_m.AtRaw(i), rhs_m);
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j) const {
            return ::pow(lhs_m.AtRaw(i, j), rhs_m);
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k) const {
            return ::pow(lhs_m.AtRaw(i, j, k), rhs_m);
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l)const {
            return ::pow(lhs_m.AtRaw(i, j, k, l), rhs_m);
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m) const {
            return ::pow(lhs_m.AtRaw(i, j, k, l, m), rhs_m);
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n) const {
            return ::pow(lhs_m.AtRaw(i, j, k, l, m, n), rhs_m);
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n, const uint32_t & o) const {
            return ::pow(lhs_m.AtRaw(i, j, k, l, m, n, o), rhs_m);
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
    struct ScalarArrayPow : ArrayExpression<typename PromoteType<typename RHS::RET_TYPE, T>::return_type, ScalarArrayPow<T, RHS> > {
        typedef typename atl::PromoteType<T, typename RHS::RET_TYPE>::return_type RET_TYPE;
        typedef typename atl::PromoteType<T, typename RHS::BASE_TYPE>::return_type BASE_TYPE;

        const T& lhs_m;
        const RHS& rhs_m;

        inline explicit ScalarArrayPow(const T& lhs, const ArrayExpression<typename RHS::RET_TYPE, RHS> & rhs) : lhs_m(lhs), rhs_m(rhs.Cast()) {



        }

        inline const size_t Size(const int32_t & dimension) const {
            return rhs_m.Size(dimension);
        }

        inline const size_t Dimensions() const {
            return rhs_m.Dimensions();
        }

        inline const RET_TYPE operator()(const uint32_t & i) const {
            return ::pow(lhs_m, rhs_m(i));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j) const {
            return ::pow(lhs_m, rhs_m(i, j));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k) const {
            return ::pow(lhs_m, rhs_m(i, j, k));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l)const {
            return ::pow(lhs_m, rhs_m(i, j, k, l));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m) const {
            return ::pow(lhs_m, rhs_m(i, j, k, l, m));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n) const {
            return ::pow(lhs_m, rhs_m(i, j, k, l, m, n));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n, const uint32_t & o) const {
            return ::pow(lhs_m, rhs_m(i, j, k, l, m, n, o));
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

        inline const RET_TYPE AtRaw(const uint32_t & i) const {
            return ::pow(lhs_m, rhs_m.AtRaw(i));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j) const {
            return ::pow(lhs_m, rhs_m.AtRaw(i, j));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k) const {
            return ::pow(lhs_m, rhs_m.AtRaw(i, j, k));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l)const {
            return ::pow(lhs_m, rhs_m.AtRaw(i, j, k, l));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m) const {
            return ::pow(lhs_m, rhs_m.AtRaw(i, j, k, l, m));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n) const {
            return ::pow(lhs_m, rhs_m.AtRaw(i, j, k, l, m, n));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n, const uint32_t & o) const {
            return ::pow(lhs_m, rhs_m.AtRaw(i, j, k, l, m, n, o));
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
    inline const ArrayPow< LHS, RHS> operator+(const MatrixExpression<typename LHS::RET_TYPE, LHS>& a,
            const MatrixExpression<typename RHS::RET_TYPE, RHS>& b) {
        return ArrayPow< LHS, RHS > (a.Cast(), b.Cast());
    }



#define ATL_ARRAY_ArrayPow_SCALAR(TYPE) \
    template< class LHS>      \
        inline const ArrayPowScalar<LHS,TYPE> \
    pow(const ArrayExpression<typename LHS::RET_TYPE, LHS>& a, const TYPE& b) {\
            return ArrayPowScalar<LHS,TYPE > (a.Cast(), b);\
        } \
    

    ATL_ARRAY_ArrayPow_SCALAR(short)
    ATL_ARRAY_ArrayPow_SCALAR(unsigned short)
    ATL_ARRAY_ArrayPow_SCALAR(int)
    ATL_ARRAY_ArrayPow_SCALAR(unsigned int)
    ATL_ARRAY_ArrayPow_SCALAR(long)
    ATL_ARRAY_ArrayPow_SCALAR(unsigned long)
    ATL_ARRAY_ArrayPow_SCALAR(float)
    ATL_ARRAY_ArrayPow_SCALAR(double)
    ATL_ARRAY_ArrayPow_SCALAR(long double)
    ATL_ARRAY_ArrayPow_SCALAR(atl::Variable<float>)
    ATL_ARRAY_ArrayPow_SCALAR(atl::Variable<double>)
    ATL_ARRAY_ArrayPow_SCALAR(atl::Variable<long double>)


#define ATL_ArrayPow_SCALAR_ARRAY(TYPE) \
    template< class RHS>      \
        inline const ScalarArrayPow<TYPE,RHS> \
    pow(const TYPE& a, const ArrayExpression<typename RHS::RET_TYPE, RHS>& b ) {\
            return ScalarArrayPow<TYPE,RHS > (a, b.Cast());\
        } \
            
    ATL_ArrayPow_SCALAR_ARRAY(short)
    ATL_ArrayPow_SCALAR_ARRAY(unsigned short)
    ATL_ArrayPow_SCALAR_ARRAY(int)
    ATL_ArrayPow_SCALAR_ARRAY(unsigned int)
    ATL_ArrayPow_SCALAR_ARRAY(long)
    ATL_ArrayPow_SCALAR_ARRAY(unsigned long)
    ATL_ArrayPow_SCALAR_ARRAY(float)
    ATL_ArrayPow_SCALAR_ARRAY(double)
    ATL_ArrayPow_SCALAR_ARRAY(long double)
    ATL_ArrayPow_SCALAR_ARRAY(atl::Variable<float>)
    ATL_ArrayPow_SCALAR_ARRAY(atl::Variable<double>)
    ATL_ArrayPow_SCALAR_ARRAY(atl::Variable<long double>)



    template< class LHS, class RHS>
    struct ArrayATan2 : ArrayExpression<typename atl::PromoteType<typename LHS::RET_TYPE, typename RHS::RET_TYPE >::return_type, ArrayATan2<LHS, RHS> > {
        typedef typename LHS::RET_TYPE RET_TYPEL;
        typedef typename RHS::RET_TYPE RET_TYPER;

        typedef typename atl::PromoteType<RET_TYPEL, RET_TYPER>::return_type RET_TYPE;
        typedef typename atl::PromoteType<typename LHS::BASE_TYPE, typename RHS::BASE_TYPE>::return_type BASE_TYPE;

        const LHS& lhs_m;
        const RHS& rhs_m;

        inline explicit ArrayATan2(const atl::ArrayExpression<typename LHS::RET_TYPE, LHS>& lhs, const atl::ArrayExpression<typename RHS::RET_TYPE, RHS>& rhs) : lhs_m(lhs.Cast()), rhs_m(rhs.Cast()) {



        }

        inline const size_t Size(const int32_t & dimension) const {
            return lhs_m.Size(dimension) < rhs_m.Size(dimension) ? lhs_m.Size(dimension) : rhs_m.Size(dimension);
        }

        inline const size_t Dimensions() const {
            return lhs_m.Dimensions() < rhs_m.Dimensions() ? lhs_m.Dimensions() : rhs_m.Dimensions();
        }

        inline const RET_TYPE operator()(const uint32_t & i) const {
            return ::atan2(lhs_m(i), rhs_m(i));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j) const {
            return ::atan2(lhs_m(i, j), rhs_m(i, j));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k) const {
            return ::atan2(lhs_m(i, j, k), rhs_m(i, j, k));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l)const {
            return ::atan2(lhs_m(i, j, k, l), rhs_m(i, j, k, l));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m) const {
            return ::atan2(lhs_m(i, j, k, l, m), rhs_m(i, j, k, l, m));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n) const {
            return ::atan2(lhs_m(i, j, k, l, m, n), rhs_m(i, j, k, l, m, n));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n, const uint32_t & o) const {
            return ::atan2(lhs_m(i, j, k, l, m, n, o), rhs_m(i, j, k, l, m, n, o));
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

        inline const RET_TYPE AtRaw(const uint32_t & i) const {
            return ::atan2(lhs_m.AtRaw(i), rhs_m.AtRaw(i));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j) const {
            return ::atan2(lhs_m.AtRaw(i, j), rhs_m.AtRaw(i, j));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k) const {
            return ::atan2(lhs_m.AtRaw(i, j, k), rhs_m.AtRaw(i, j, k));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l)const {
            return ::atan2(lhs_m.AtRaw(i, j, k, l), rhs_m.AtRaw(i, j, k, l));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m) const {
            return ::atan2(lhs_m.AtRaw(i, j, k, l, m), rhs_m.AtRaw(i, j, k, l, m));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n) const {
            return ::atan2(lhs_m.AtRaw(i, j, k, l, m, n), rhs_m.AtRaw(i, j, k, l, m, n));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n, const uint32_t & o) const {
            return ::atan2(lhs_m.AtRaw(i, j, k, l, m, n, o), rhs_m.AtRaw(i, j, k, l, m, n, o));
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
    struct ArrayATan2Scalar : ArrayExpression<typename PromoteType<typename LHS::RET_TYPE, T>::return_type, ArrayATan2Scalar< LHS, T> > {
        const LHS& lhs_m;
        const T& rhs_m;

        typedef typename atl::PromoteType<typename LHS::RET_TYPE, T>::return_type RET_TYPE;
        typedef typename atl::PromoteType<typename LHS::BASE_TYPE, T>::return_type BASE_TYPE;

        inline explicit ArrayATan2Scalar(const ArrayExpression<typename LHS::RET_TYPE, LHS>& lhs, const T & rhs) : lhs_m(lhs.Cast()), rhs_m(rhs) {



        }

        inline const size_t Size(const int32_t & dimension) const {
            return lhs_m.Size(dimension);
        }

        inline const size_t Dimensions() const {
            return lhs_m.Dimensions();
        }

        inline const RET_TYPE operator()(const uint32_t & i) const {
            return ::atan2(lhs_m(i), rhs_m);
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j) const {
            return ::atan2(lhs_m(i, j), rhs_m);
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k) const {
            return ::atan2(lhs_m(i, j, k), rhs_m);
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l)const {
            return ::atan2(lhs_m(i, j, k, l), rhs_m);
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m) const {
            return ::atan2(lhs_m(i, j, k, l, m), rhs_m);
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n) const {
            return ::atan2(lhs_m(i, j, k, l, m, n), rhs_m);
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n, const uint32_t & o) const {
            return ::atan2(lhs_m(i, j, k, l, m, n, o), rhs_m);
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

        inline const RET_TYPE AtRaw(const uint32_t & i) const {
            return ::atan2(lhs_m.AtRaw(i), rhs_m);
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j) const {
            return ::atan2(lhs_m.AtRaw(i, j), rhs_m);
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k) const {
            return ::atan2(lhs_m.AtRaw(i, j, k), rhs_m);
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l)const {
            return ::atan2(lhs_m.AtRaw(i, j, k, l), rhs_m);
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m) const {
            return ::atan2(lhs_m.AtRaw(i, j, k, l, m), rhs_m);
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n) const {
            return ::atan2(lhs_m.AtRaw(i, j, k, l, m, n), rhs_m);
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n, const uint32_t & o) const {
            return ::atan2(lhs_m.AtRaw(i, j, k, l, m, n, o), rhs_m);
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
    struct ScalarArrayATan2 : ArrayExpression<typename PromoteType<typename RHS::RET_TYPE, T>::return_type, ScalarArrayATan2<T, RHS> > {
        typedef typename atl::PromoteType<T, typename RHS::RET_TYPE>::return_type RET_TYPE;
        typedef typename atl::PromoteType<T, typename RHS::BASE_TYPE>::return_type BASE_TYPE;

        const T& lhs_m;
        const RHS& rhs_m;

        inline explicit ScalarArrayATan2(const T& lhs, const ArrayExpression<typename RHS::RET_TYPE, RHS> & rhs) : lhs_m(lhs), rhs_m(rhs.Cast()) {



        }

        inline const size_t Size(const int32_t & dimension) const {
            return rhs_m.Size(dimension);
        }

        inline const size_t Dimensions() const {
            return rhs_m.Dimensions();
        }

        inline const RET_TYPE operator()(const uint32_t & i) const {
            return ::atan2(lhs_m, rhs_m(i));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j) const {
            return ::atan2(lhs_m, rhs_m(i, j));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k) const {
            return ::atan2(lhs_m, rhs_m(i, j, k));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l)const {
            return ::atan2(lhs_m, rhs_m(i, j, k, l));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m) const {
            return ::atan2(lhs_m, rhs_m(i, j, k, l, m));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n) const {
            return ::atan2(lhs_m, rhs_m(i, j, k, l, m, n));
        }

        inline const RET_TYPE operator()(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n, const uint32_t & o) const {
            return ::atan2(lhs_m, rhs_m(i, j, k, l, m, n, o));
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

        inline const RET_TYPE AtRaw(const uint32_t & i) const {
            return ::atan2(lhs_m, rhs_m.AtRaw(i));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j) const {
            return ::atan2(lhs_m, rhs_m.AtRaw(i, j));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k) const {
            return ::atan2(lhs_m, rhs_m.AtRaw(i, j, k));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l)const {
            return ::atan2(lhs_m, rhs_m.AtRaw(i, j, k, l));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m) const {
            return ::atan2(lhs_m, rhs_m.AtRaw(i, j, k, l, m));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n) const {
            return ::atan2(lhs_m, rhs_m.AtRaw(i, j, k, l, m, n));
        }

        inline const RET_TYPE AtRaw(const uint32_t& i, const uint32_t & j, const uint32_t & k, const uint32_t & l, const uint32_t & m, const uint32_t & n, const uint32_t & o) const {
            return ::atan2(lhs_m, rhs_m.AtRaw(i, j, k, l, m, n, o));
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
    inline const ArrayATan2< LHS, RHS> operator+(const MatrixExpression<typename LHS::RET_TYPE, LHS>& a,
            const MatrixExpression<typename RHS::RET_TYPE, RHS>& b) {
        return ArrayATan2< LHS, RHS > (a.Cast(), b.Cast());
    }



#define ATL_ARRAY_ArrayATan2_SCALAR(TYPE) \
    template< class LHS>      \
        inline const ArrayATan2Scalar<LHS,TYPE> \
    atan2(const ArrayExpression<typename LHS::RET_TYPE, LHS>& a, const TYPE& b) {\
            return ArrayATan2Scalar<LHS,TYPE > (a.Cast(), b);\
        } \
    
    ATL_ARRAY_ArrayATan2_SCALAR(short)
    ATL_ARRAY_ArrayATan2_SCALAR(unsigned short)
    ATL_ARRAY_ArrayATan2_SCALAR(int)
    ATL_ARRAY_ArrayATan2_SCALAR(unsigned int)
    ATL_ARRAY_ArrayATan2_SCALAR(long)
    ATL_ARRAY_ArrayATan2_SCALAR(unsigned long)
    ATL_ARRAY_ArrayATan2_SCALAR(float)
    ATL_ARRAY_ArrayATan2_SCALAR(double)
    ATL_ARRAY_ArrayATan2_SCALAR(long double)
    ATL_ARRAY_ArrayATan2_SCALAR(atl::Variable<float>)
    ATL_ARRAY_ArrayATan2_SCALAR(atl::Variable<double>)
    ATL_ARRAY_ArrayATan2_SCALAR(atl::Variable<long double>)


#define ATL_ArrayATan2_SCALAR_ARRAY(TYPE) \
    template< class RHS>      \
        inline const ScalarArrayATan2<TYPE,RHS> \
    atan2(const TYPE& a, const ArrayExpression<typename RHS::RET_TYPE, RHS>& b ) {\
            return ScalarArrayATan2<TYPE,RHS > (a, b.Cast());\
        } \
            
    ATL_ArrayATan2_SCALAR_ARRAY(short)
    ATL_ArrayATan2_SCALAR_ARRAY(unsigned short)
    ATL_ArrayATan2_SCALAR_ARRAY(int)
    ATL_ArrayATan2_SCALAR_ARRAY(unsigned int)
    ATL_ArrayATan2_SCALAR_ARRAY(long)
    ATL_ArrayATan2_SCALAR_ARRAY(unsigned long)
    ATL_ArrayATan2_SCALAR_ARRAY(float)
    ATL_ArrayATan2_SCALAR_ARRAY(double)
    ATL_ArrayATan2_SCALAR_ARRAY(long double)
    ATL_ArrayATan2_SCALAR_ARRAY(atl::Variable<float>)
    ATL_ArrayATan2_SCALAR_ARRAY(atl::Variable<double>)
    ATL_ArrayATan2_SCALAR_ARRAY(atl::Variable<long double>)





}
#endif