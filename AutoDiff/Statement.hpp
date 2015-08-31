/* 
 * File:   Statement.hpp
 * Author: matthewsupernaw
 *
 * Created on June 18, 2014, 1:16 PM
 */

#ifndef ET4AD_STATEMENT_HPP
#define	ET4AD_STATEMENT_HPP
#include <iostream>

/**
 * Operation values used for recording expressions into a post-order 
 * expression tree. These operations are used primarily for supporting 
 * operations such as arbitray order derivatives,uncertainty calculations
 * and expression string building.
 */
enum Operation {
    CONSTANT = 0,
    VARIABLE,
    PLUS,
    MINUS,
    MULTIPLY,
    DIVIDE,
    SIN,
    COS,
    TAN,
    ASIN,
    ACOS,
    ATAN,
    ATAN2,
    ATAN3,
    ATAN4,
    SQRT,
    POW,
    POW1,
    POW2,
    LOG,
    LOG10,
    EXP,
    MFEXP,
    SINH,
    COSH,
    TANH,
    ABS,
    FABS,
    FLOOR,
    CEIL,

    PLUS_EQUALS,
    MINUS_EQUALS,
    TIMES_EQUALS,
    DIVIDE_EQUALS,
    NONE
};

/**
 * Statement class is used store the expression in a post-order vector. 
 * Ultimately used for higher order derivatives and expression string 
 * building.
 * 
 * @param op
 */
template<class REAL_T>
struct Statement {
    typedef REAL_T BASE_TYPE;

    Statement() : op_m(CONSTANT), value_m(0), id_m(0) {

    }

    Statement(const Operation &op) : op_m(op), value_m(0), id_m(0) {

    }

    Statement(const Operation &op, const REAL_T &value) : op_m(op), value_m(value), id_m(0) {

    }

    Statement(const Operation &op, const REAL_T &value, const uint32_t &id) : op_m(op), value_m(value), id_m(id) {

    }

    Statement(const Statement &orig) : op_m(orig.op_m), value_m(orig.value_m), id_m(orig.id_m) {

    }

    uint32_t GetId() const {
        return id_m;
    }

    void SetId(uint32_t id) {
        this->id_m = id;
    }

    Operation GetOp() const {
        return op_m;
    }

    void SetOp(Operation op) {
        this->op_m = op;
    }

    REAL_T GetValue() const {
        return value_m;
    }

    void SetValue(REAL_T value) {
        this->value_m = value;
    }



    //private:
    Operation op_m;
    REAL_T value_m;
    uint32_t id_m;

};

template<class REAL_T>
std::ostream& operator<<(std::ostream& out, const Statement<REAL_T> & s) {
    out << "statement[" << s.op_m << ", " << s.value_m << ", " << s.id_m << "]\n";
    return out;
}

#endif	/* STATEMENT_HPP */

