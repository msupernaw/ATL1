

#ifndef FLAT_SET_HPP
#define	FLAT_SET_HPP

template <class T, class Compare = std::less<T>, class Allocator = std::allocator<T> >
class flat_set {
    std::vector<T, Allocator> data_m;
    Compare cmp;
public:
    typedef typename std::vector<T>::iterator iterator;
    typedef typename std::vector<T>::const_iterator const_iterator;
    typedef typename std::vector<T>::reverse_iterator reverse_iterator;
    typedef typename std::vector<T>::const_reverse_iterator const_reverse_iterator;
    
    flat_set(const Compare& c = Compare())
    : data_m(), cmp(c) {
        //        data_m.reserve(32);
    }
    
    template <class InputIterator>
    flat_set(InputIterator first, InputIterator last,
             const Compare& c = Compare())
    : data_m(first, last), cmp(c) {
        std::sort(begin(), end(), cmp);
    }
    
    inline iterator erase(const T& t) {
        iterator i = this->find(t);
        if(i != this->end()){
            return data_m.erase(i);
        }else{
            return this->end();
        }
    }
    
    inline iterator insert(const T& t) {
        iterator i = std::lower_bound(begin(), end(), t, cmp);
        if (i == end() || cmp(t, *i))
            data_m.insert(i, t);
        return i;
    }
    
    inline const_iterator find(const T& t) const {
        const_iterator i = std::lower_bound(begin(), end(), t, cmp);
        return i == end() || cmp(t, *i) ? end() : i;
    }
    
    inline iterator find(const T& t)  {
        iterator i = std::lower_bound(begin(), end(), t, cmp);
        return i == end() || cmp(t, *i) ? end() : i;
    }
    
    inline iterator begin() {
        return data_m.begin();
    }
    
    inline iterator end() {
        return data_m.end();
    }
    
    inline const_iterator begin() const noexcept {
        return data_m.begin();
    }
    
    inline const_iterator end() const noexcept {
        return data_m.end();
    }
    
    inline reverse_iterator rbegin() {
        return data_m.rbegin();
    }
    
    inline reverse_iterator rend() {
        return data_m.rend();
    }
    
    inline const_reverse_iterator rbegin() const {
        return data_m.rbegin();
    }
    
    inline const_reverse_iterator rend() const {
        return data_m.rend();
    }
    
    inline size_t size() const{
        return data_m.size();
    }
    
    inline bool contains(const T& t){
        return std::binary_search(this->data_m.begin(), this->data_m.end(),t);
    }
    
    inline void clear() {
        data_m.clear();
    }
    
    
    
};




#endif	/* FLAT_SET_HPP */

