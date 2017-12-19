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
 #  Logger Class Definition & Implementation
*/

#ifndef _RACD_LOGGER_
#define _RACD_LOGGER_

#include <fstream>
#include <string>
#include <iostream>

/* simple logging singleton
 * this logger must be initialized with logger::open_log()
 * before any other instances are made.
 */
class logger {
public:

  void open_log(const std::string& _out_trans){
    out_trans.open(_out_trans);
  };

  void log_trans(const std::string& trans){
    out_trans << trans << std::endl;
  };

  void suicide(){
    #ifdef DEBUG_HPP
    std::cout << "logger suiciding at " << this << std::endl;
    #endif
    out_trans.close();
  };

  static logger* instance(){
    if(!l_instance){
      l_instance = new logger;
    }
    return l_instance;
  };

private:
  std::ofstream out_trans; /* all state transitions */
  std::ofstream out_pop; /* population size */
  std::ofstream out_clinic; /* clinical incidence */
  std::ofstream out_slide; /* slide positivity */

  static logger* l_instance;

  /* private constructor */
  logger(){
    #ifdef DEBUG_HPP
    std::cout << "logger being born at " << this << std::endl;
    #endif
  };

  /* private destructor */
  ~logger(){
    #ifdef DEBUG_HPP
    std::cout << "logger being killed at " << this << std::endl;
    #endif
  };
};

logger *logger::l_instance = nullptr;

#endif
