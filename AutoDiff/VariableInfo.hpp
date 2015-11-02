/* 
 * File:   DerivativePointer.hpp
 * Author: matthewsupernaw
 *
 * Created on April 8, 2015, 12:12 PM
 */

#ifndef DERIVATIVEPOINTER_HPP
#define	DERIVATIVEPOINTER_HPP
#include <atomic>
#include <mutex>
#include <stack>
#include <memory>
#include "../Utilities/Platform.hpp"
#include "../Utilities/flat_set.hpp"


namespace atl {

    /*!
     * Creates a unique identifier for variables. Identifiers are recyclable.
     * @return 
     */
    class VariableIdGenerator {
        std::stack<uint32_t> available;
        std::atomic<uint32_t> available_size;
        static std::mutex mutex_g;

    public:
        static std::shared_ptr<VariableIdGenerator> instance();

        const uint32_t next() {
            mutex_g.lock();
            uint32_t ret;
            if (!available.empty() > 0) {
                ret = available.top();
                available.pop();

            } else {
                ret = ++_id;
            }



            mutex_g.unlock();
            return ret; //(++_id);
        }

        void release(const uint32_t& id) {
            mutex_g.lock();
            available.push(id);
            available_size++;
            mutex_g.unlock();
        }

        const uint32_t current() {
            return _id;
        }

        //    private:

        VariableIdGenerator() : _id(0), available_size(0) {
        }

        std::atomic<uint32_t> _id;
    };

    std::mutex VariableIdGenerator::mutex_g;
    static std::shared_ptr<VariableIdGenerator> only_copy;

    inline std::shared_ptr<VariableIdGenerator>
    VariableIdGenerator::instance() {

        if (!only_copy) {
            only_copy = std::make_shared<VariableIdGenerator>();
        }

        return only_copy;
    }

    /**
     *quick and dirty hash table to hold computed Hessian info.
     */
    template<typename T, uint32_t BUCKET_SIZE = 100 >
    class hash_table {
    public:
        uint32_t default_number_of_buckets;
        uint32_t current_number_of_buckets;
        std::vector<std::vector<T > > buckets;

        typedef T value_type;
        typedef T* iterator;
        typedef const T* const_iterator;

        hash_table(uint32_t default_number_of_buckets = 1) : default_number_of_buckets(default_number_of_buckets) {
            buckets.resize(default_number_of_buckets);
            current_number_of_buckets = default_number_of_buckets;
        }

        inline void insert(uint32_t key, const value_type& value) {
            uint32_t bucket = key / BUCKET_SIZE;
            if (unsigned(bucket) >= unsigned(current_number_of_buckets)) {
                buckets.resize(bucket + 1);
                current_number_of_buckets = bucket + 1;
            }

            uint32_t local_key = key - (BUCKET_SIZE * bucket);
            if (buckets[bucket].empty()) {
                buckets[bucket].resize(BUCKET_SIZE);
            }

            buckets[bucket][local_key] = value;
        }

        inline T& operator[](const uint32_t& key) {

            uint32_t bucket = key / BUCKET_SIZE;

            if (unsigned(bucket) >= unsigned(current_number_of_buckets)) {
                buckets.resize(bucket + 1);
                current_number_of_buckets = bucket + 1;
            }

            uint32_t local_key = key - (BUCKET_SIZE * bucket);
            if (buckets[bucket].empty()) {
                buckets[bucket].resize(BUCKET_SIZE);
            }
            return buckets[bucket][local_key];
        }

        inline const_iterator find(const uint32_t& key) const {
            uint32_t bucket = key / BUCKET_SIZE;
            if (unsigned(bucket) >= unsigned(current_number_of_buckets)) {
                return this->end();
            }

            if (buckets[bucket].empty()) {
                return this->end();
            }
            uint32_t local_key = key - (BUCKET_SIZE * bucket);
            return &(*(buckets[bucket].begin() + local_key));

        }

        inline iterator find(const uint32_t& key) {
            uint32_t bucket = unsigned(key) / BUCKET_SIZE;
            if (unsigned(bucket) >= unsigned(current_number_of_buckets)) {
                return this->end();
            }

            if (buckets[bucket].empty()) {
                return this->end();
            }
            uint32_t local_key = key - (BUCKET_SIZE * bucket);
            return &(*(buckets[bucket].begin() + local_key));

        }

        inline const_iterator end() const {
            return NULL; //buckets[buckets.size() - 1].end();
        }

        inline iterator end() {
            return NULL; //*(buckets[buckets.size() - 1].end());
        }

        inline void clear() {
            buckets.resize(default_number_of_buckets);
            buckets[0].clear();
            current_number_of_buckets = default_number_of_buckets;
        }


    };

    template<typename REAL_T>
    class VariableInfo {
    public:
        static std::mutex mutex_g;
        static std::mutex vinfo_mutex_g;
        static std::vector<VariableInfo<REAL_T>* > freed;
        REAL_T dvalue;
        REAL_T vvalue;
        std::atomic<uint32_t> count;
        typedef hash_table<REAL_T> HessianInfo;
        HessianInfo hessian_row;
        typedef typename HessianInfo::iterator row_iterator;
        uint32_t id;

        VariableInfo() : vvalue(0.0), dvalue(0.0), count(1), id(VariableIdGenerator::instance()->next()) {
            //            hessian_row.set_empty_key(0);

        }

        inline void Aquire() {
            count++;
        }

        inline REAL_T GetHessianRowValue(VariableInfo<REAL_T>* var) {
            row_iterator it = hessian_row.find(var);
            if (it != hessian_row.end()) {
                return (*it).second;
            } else {
                return static_cast<REAL_T> (0.0);
            }
        }

        inline void Release() {
            count--;
            if ((count) == 0) {
                VariableInfo<REAL_T>::vinfo_mutex_g.lock();
                freed.push_back(this);
                VariableInfo<REAL_T>::vinfo_mutex_g.unlock();
            }
        }

        static void FreeAll() {
            for (int i = 0; i < freed.size(); i++) {
                VariableIdGenerator::instance()->release(freed[i]->id);
                delete freed[i];
            }
            freed.resize(0);
        }

        inline void Reset() {
            this->hessian_row.clear(); // = hash_table<REAL_T>();
            this->dvalue = 0;
        }

        void* operator new(size_t size) {
            return malloc(size);
        }

        void operator delete(void* ptr) {
            free(ptr);
        }
    };


    template<typename REAL_T>
    std::vector<VariableInfo<REAL_T>* > VariableInfo<REAL_T>::freed;

    template<typename REAL_T>
    std::mutex VariableInfo<REAL_T>::vinfo_mutex_g;

    template<typename REAL_T>
    std::mutex VariableInfo<REAL_T>::mutex_g;


}


#endif	/* DERIVATIVEPOINTER_HPP */

