/*
 #      ____  ___   __________
 #     / __ \/   | / ____/ __ \
 #    / /_/ / /| |/ /   / / / /
 #   / _, _/ ___ / /___/ /_/ /
 #  /_/ |_/_/  |_\____/_____/
 #
 #  Sean Wu
 #  slwu89@berkeley.edu
 #  September 2019
 #
 #  stuff the whole program needs to see
*/

#include "output.hpp"


/* constructor & destructor */
output::output(){};
output::~output(){};

/* utility methods */
output& output::instance(){
    static output instance;
    return instance;
};
