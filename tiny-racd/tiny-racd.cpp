/*
 #      ____  ___   __________
 #     / __ \/   | / ____/ __ \
 #    / /_/ / /| |/ /   / / / /
 #   / _, _/ ___ / /___/ /_/ /
 #  /_/ |_/_/  |_\____/_____/
 #
 #  Sean Wu & John M. Marshall
 #  April 2019
 #
 #  the whole enchilada
*/

/* ################################################################################
#   Includes & auxiliary
################################################################################ */

// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::depends(RcppProgress)]]

// for Rcpp and R's RNG because we're lazy
#include <Rcpp.h>
#include <Rmath.h>
#include <progress.hpp>

// the best debugger
#include <iostream>

// for holding the people
#include <list>

// name one C++ project that doesn't have vectors in it...go!
#include <vector>

// for smart pointers
#include <memory>

// for holding the functions
#include <unordered_map>
#include <functional>

// for strings of course
#include <string>

// to pretend we know CS
#include <algorithm>

// for da maths
#include <math.h>


/* ################################################################################
#   output and parameters/global stuff
################################################################################ */

// current simulation time
static size_t tnow = 0;

// parameters
static std::unordered_map<std::string,double> parameters;


/* ################################################################################
#   House
################################################################################ */

// a house
// all expectations are averaging out heterogeneity in the people at this house
typedef struct house {

  // data members
  size_t            id;
  double            psi;  // P(a bite going anywhere goes here)
  double            w;    // E[P(feed and survive)]
  double            y;    // E[P(feed)]
  double            z;    // E[P(repelled without feeding)]
  double            c;    // E[P(a feed here will result in a mosquito infection)]
  size_t            n;    // number of people here
  bool              IRS;  // does my house have IRS

  house(const size_t id_, const double psi_, const double w_, const double y_, const double z_, const double c_,
        const size_t n_
  );
  ~house();
} house;

// constructor
house::house(const size_t id_, const double psi_, const double w_, const double y_, const double z_, const double c_,
             const size_t n_
) :
  id(id_), psi(psi_), w(w_), y(y_), z(z_), c(c_), n(n_), IRS(false)
{};

// destructor
house::~house(){};

// we're men of taste, who use smart pointers
using house_ptr = std::unique_ptr<house>;

// the houses!
static std::vector<house_ptr> houses;
