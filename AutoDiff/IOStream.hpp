/* 
 * File:   IOStream.hpp
 * Author: matthewsupernaw
 *
 * Created on April 1, 2014, 11:06 AM
 */

#ifndef IOSTREAM_HPP
#define	IOSTREAM_HPP

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>


namespace util {

    /* Convert double to string with specified number of places after the decimal. */
    std::string prd(const double &x, const int decDigits) {
        std::stringstream ss;
        ss << std::fixed;
        ss.precision(decDigits); // set # places after decimal
        ss << x;
        return ss.str();
    }

    /* Convert double to string with specified number of places after the decimal
       and left padding. */
    std::string prd(const double &x, const int decDigits, const int width) {
        std::stringstream ss;
        ss << std::fixed << std::right;
        ss.fill(' '); // fill space around displayed #
        ss.width(width); // set  width around displayed #
        ss.precision(decDigits); // set # places after decimal
        ss << x;
        return ss.str();
    }

    /*! Center-aligns string within a field of width w. Pads with blank spaces
        to enforce alignment. */
    std::string center(const std::string &s, const int w) {
        std::stringstream ss, spaces;
        int padding = w - s.size(); // count excess room to pad
        for (int i = 0; i < padding / 2; ++i)
            spaces << " ";
        ss << spaces.str() << s << spaces.str(); // format with padding
        if (padding > 0 && padding % 2 != 0) // if odd #, add 1 space
            ss << " ";
        return ss.str();
    }

    /**
     * Trim the left side.
     * 
     * @param s
     * @return 
     */
    static inline std::string& LeftTrim(std::string &s) {
        s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
        return s;
    }

    /**
     * Trim the right side.
     * 
     * @param s
     * @return 
     */
    static inline std::string& RightTrim(std::string &s) {
        s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
        return s;
    }

    /**
     * Trim both sides.
     * 
     * @param s
     * @return 
     */
    static inline std::string& Trim(std::string &s) {
        return LeftTrim(RightTrim(s));
    }

    static void Tokenize(const std::string& str, std::vector<std::string>& tokens,
            const std::string& delimiters = " ", const bool trimEmpty = true) {
        std::string::size_type pos, lastPos = 0;
        while (true) {
            pos = str.find_first_of(delimiters, lastPos);
            if (pos == std::string::npos) {
                pos = str.length();

                if (pos != lastPos || !trimEmpty)
                    tokens.push_back(std::vector<std::string>::value_type(str.data() + lastPos,
                        (std::vector<std::string>::value_type::size_type)pos - lastPos));

                break;
            } else {
                if (pos != lastPos || !trimEmpty)
                    tokens.push_back(std::vector<std::string>::value_type(str.data() + lastPos,
                        (std::vector<std::string>::value_type::size_type)pos - lastPos));
            }

            lastPos = pos + 1;
        }
    };

    static bool StartsWith(const std::string &value1, const std::string &value2) {
        return value1.find(value2) == 0;
    }

    template <typename T>
    T StringToNumber(const std::string &Text) {
        std::istringstream ss(Text);
        T result;
        return ss >> result ? result : 0;
    }
}

/**
 * Simple class for reading space and tab delimited 
 * streamed input.
 */
template<class T>
class StreamedDataFile {
public:

    StreamedDataFile() {
    }

    StreamedDataFile(const std::string& path) {
        this->Parse(path);
    }

    StreamedDataFile(const StreamedDataFile<T>& orig) {
    }

    ~StreamedDataFile() {
    }

    void Parse(std::string path) {
        this->path_m = path;
        this->data_index_m = 0;
        this->data_m.clear();

        std::ifstream input;
        input.open(path.c_str());
        std::string line;
        this->data_m.reserve(10000);
        while (input.good()) {
            std::getline(input, line);

            if (!util::StartsWith(line, "#")) {
                std::vector<std::string> tokens;

                line = util::LeftTrim(line);
                line = util::RightTrim(line);
                util::Tokenize(line, tokens, " \t", true);

                for (int i = 0; i < tokens.size(); i++) {
                    if (util::StartsWith(tokens.at(i), "#"))break;
                    this->data_m.push_back(util::StringToNumber<T > (tokens.at(i)));
                }


            }

        }
    }

    /**
     * Returns the next value in the stream.
     * 
     * @return 
     */
    T Next() {
        if (this->data_index_m< this->data_m.size()) {
            this->data_index_m++;
            return this->data_m.at(this->data_index_m - 1);
        }
        return -999;
    }

    /**
     * True if there exist more values in 
     * the data stream, else false.
     * 
     * @return 
     */
    bool HasMore() {
        return (this->data_index_m < data_m.size());
    }

    std::string GetPath() const {
        return path_m;
    }



private:
    std::string path_m;
    std::vector<T> data_m;
    size_t data_index_m;
};


#endif	/* IOSTREAM_HPP */

