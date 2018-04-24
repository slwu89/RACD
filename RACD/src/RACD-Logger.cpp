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

#include "RACD-Logger.hpp"

/* constructor & destructor */
logger::logger(){
  #ifdef DEBUG_RACD
  std::cout << "logger being born at " << this << std::endl;
  #endif
};

logger::~logger(){
  #ifdef DEBUG_RACD
  std::cout << "logger being killed at " << this << std::endl;
  #endif
};

/* utility methods */
logger& logger::instance(){
    static logger instance;
    return instance;
};

/* open logging file */
void logger::open_log(const std::string &_out_trans){
  out_trans.open(_out_trans);
}

/* log data */
void logger::log_trans(const std::string &trans){
  out_trans << trans << "\n";
}

/* close open streams */
void logger::close_log(){
  out_trans.close();
}
