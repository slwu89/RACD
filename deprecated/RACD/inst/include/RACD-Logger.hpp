/*
 #      ____  ___   __________
 #     / __ \/   | / ____/ __ \
 #    / /_/ / /| |/ /   / / / /
 #   / _, _/ ___ / /___/ /_/ /
 #  /_/ |_/_/  |_\____/_____/
 #
 #  Sean Wu & John M. Marshall
 #  December 2017
 #
 #  Logger Class Declaration
*/

/* ######################################################################
 # includes and foward declarations
###################################################################### */

#ifndef RACD_LOGGER
#define RACD_LOGGER

/* C++ includes */
#include <fstream>
#include <string>
#include <iostream>

// #include "DEBUG.hpp"


/* ######################################################################
 # class declaration
###################################################################### */

class logger {
public:

  /* constructor & destructor */
  logger();
  ~logger();

  /* delete all copy & move semantics */
  logger(const logger&) = delete;
  logger& operator=(const logger&) = delete;
  logger(logger&&) = delete;
  logger& operator=(logger&&) = delete;

  /* interface methods */
  void open_log(const std::string& outfile_);
  std::ofstream& get_log(){return outfile;};
  void close_log();

private:
  std::ofstream outfile;
};

#endif
