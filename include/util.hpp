#ifndef UTIL_HPP
#define UTIL_HPP

#include <iostream>
#include <cstdlib>
#include <string>

inline void error_exit(const std::string& msg) { std::cerr << "Error: " << msg << std::endl; std::exit(1); }
inline void warning(const std::string& msg) { std::cerr << "Warning: " << msg << std::endl; }
inline void note(const std::string& msg) { std::cerr << "Note: " << msg << std::endl; }

#endif