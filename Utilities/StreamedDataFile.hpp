#ifndef STREAMEDDATAFILE_HPP
#define STREAMEDDATAFILE
#include <iostream>
#include <sstream>
#include <limits>
#include <fstream>


namespace atl{

class StreamedDataFile {
    std::ifstream input;
    
public:
    
    StreamedDataFile() {
        
    }
    
    StreamedDataFile(const char* path) {
        input.open(path);
    }
    
    void open(const char* path) {
        input.open(path);
    }
    
    bool good() {
        return input.good();
    }
    
    void operator>>(char* v) {
        for (;;) {
            
            input >> v;
            
            if (input.eof() || input.bad()) {
                break;
            } else if (input.fail()) {
                input.clear(); // unset failbit
                input.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // skip bad input
            } else {
                break;
            }
        }
    }
    
    void operator>>(unsigned char* v) {
        return operator>>((char*) v);
    }
    
    void operator>>(signed char*v) {
        return operator>>((char*) v);
    }
    
    void operator>>(char& v) {
        for (;;) {
            
            input >> v;
            
            if (input.eof() || input.bad()) {
                break;
            } else if (input.fail()) {
                input.clear(); // unset failbit
                input.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // skip bad input
            } else {
                break;
            }
        }
    }
    
    void operator>>(unsigned char& v) {
        return operator>>((char&) v);
    }
    
    void operator>>(signed char& v) {
        return operator>>((char&) v);
    }
    
    void operator>>(int& v) {
        for (;;) {
            
            input >> v;
            
            if (input.eof() || input.bad()) {
                break;
            } else if (input.fail()) {
                input.clear(); // unset failbit
                input.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // skip bad input
            } else {
                break;
            }
        }
    }
    
    void operator>>(long& v) {
        for (;;) {
            
            input >> v;
            
            if (input.eof() || input.bad()) {
                break;
            } else if (input.fail()) {
                input.clear(); // unset failbit
                input.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // skip bad input
            } else {
                break;
            }
        }
    }
    
    void operator>>(short& v) {
        for (;;) {
            
            input >> v;
            
            if (input.eof() || input.bad()) {
                break;
            } else if (input.fail()) {
                input.clear(); // unset failbit
                input.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // skip bad input
            } else {
                break;
            }
        }
    }
    
    void operator>>(unsigned int& v) {
        for (;;) {
            
            input >> v;
            
            if (input.eof() || input.bad()) {
                break;
            } else if (input.fail()) {
                input.clear(); // unset failbit
                input.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // skip bad input
            } else {
                break;
            }
        }
    }
    
    void operator>>(unsigned long& v) {
        for (;;) {
            
            input >> v;
            
            if (input.eof() || input.bad()) {
                break;
            } else if (input.fail()) {
                input.clear(); // unset failbit
                input.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // skip bad input
            } else {
                break;
            }
        }
    }
    
    void operator>>(unsigned short& v) {
        for (;;) {
            
            input >> v;
            
            if (input.eof() || input.bad()) {
                break;
            } else if (input.fail()) {
                input.clear(); // unset failbit
                input.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // skip bad input
            } else {
                break;
            }
        }
    }
    
    void operator>>(bool& v) {
        for (;;) {
            
            input >> v;
            
            if (input.eof() || input.bad()) {
                break;
            } else if (input.fail()) {
                input.clear(); // unset failbit
                input.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // skip bad input
            } else {
                break;
            }
        }
    }
    
    void operator>>(float& v) {
        for (;;) {
            
            input >> v;
            
            if (input.eof() || input.bad()) {
                break;
            } else if (input.fail()) {
                input.clear(); // unset failbit
                input.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // skip bad input
            } else {
                break;
            }
        }
    }
    
    void operator>>(double& v) {
        for (;;) {
            
            input >> v;
            
            if (input.eof() || input.bad()) {
                break;
            } else if (input.fail()) {
                input.clear(); // unset failbit
                input.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // skip bad input
            } else {
                break;
            }
        }
    }
    
    void operator>>(long double& v) {
        for (;;) {
            
            input >> v;
            
            if (input.eof() || input.bad()) {
                break;
            } else if (input.fail()) {
                input.clear(); // unset failbit
                input.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // skip bad input
            } else {
                break;
            }
        }
    }
private:
    
};

}


#endif