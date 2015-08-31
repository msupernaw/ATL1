/* 
 * File:   VectorFlip.hpp
 * Author: matthewsupernaw
 *
 * Created on December 15, 2014, 8:58 AM
 */

#ifndef VECTORFLIP_HPP
#define	VECTORFLIP_HPP


namespace atl {

    template<class V>
    struct VectorFlip :
    VectorExpression<V::RET_TYPE, VectorFlip<V> > {
        uint32_t size;
        const V& vec;
    };


}


#endif	/* VECTORFLIP_HPP */

