// /*
//  #      ____  ___   __________
//  #     / __ \/   | / ____/ __ \
//  #    / /_/ / /| |/ /   / / / /
//  #   / _, _/ ___ / /___/ /_/ /
//  #  /_/ |_/_/  |_\____/_____/
//  #
//  #  Sean Wu & John M. Marshall
//  #  December 2017
//  #
//  #  Logger Class Definition & Implementation
// */
//
// #include "RACD-Logger.hpp"
//
// /* constructor & destructor */
// logger::logger(){
//   #ifdef DEBUG_RACD
//   std::cout << "logger being born at " << this << std::endl;
//   #endif
// };
//
// logger::~logger(){
//   close_log();
//   #ifdef DEBUG_RACD
//   std::cout << "logger being killed at " << this << std::endl;
//   #endif
// };
//
// /* open logging file */
// void logger::open_log(const std::string &outfile_){
//   outfile.open(outfile_);
// }
//
// /* close open streams */
// void logger::close_log(){
//   outfile.close();
// }
