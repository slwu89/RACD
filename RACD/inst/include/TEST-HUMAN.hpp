#ifndef _RACD_TESTHUMAN_
#define _RACD_TESTHUMAN_

#include <string>
#include <iostream>
#include "RACD-Parameters.hpp"

class test_human {
public:
  test_human(const std::string &_name) : name(_name) {
    std::cout << name << " is being born at " << this << std::endl;
  };
  ~test_human(){
    std::cout << name << " is being killed at " << this << std::endl;
  };

  void get_kappaD(){
    std::cout << name << " is getting kappaD: " << RACD_Parameters::instance()->get_kappaD() << std::endl;
  };

  void suicide(){
    std::cout << name << " is suiciding" << std::endl;
    delete this;
  };

  friend std::ostream& operator<<(std::ostream& os, test_human& h) {
    return os << "test_human " << h.name << " says hi"  << std::endl;
}

private:
  std::string           name;
};

#endif
