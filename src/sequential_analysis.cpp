#include <boost/math/distributions/beta.hpp>
#include <math.h>
#include <numeric>
#include "sequential_analysis.hpp"

ConnectivityMatrix::ConnectivityMatrix(
  const std::vector<std::string> origins,
  const std::vector<std::string> destinations) {
  /* Save the origins and destinations. */
  this->origins = origins;
  this->destinations = destinations;
  /* Resize the count matrix. */
  this->counts.resize(origins.size());
  for(uint i = 0; i < origins.size(); ++i)
    this->counts[i].resize(destinations.size());
  /* Set delta and pi. */
  delta = 0.005;
  pi = 0.05;
}

void ConnectivityMatrix::update(const std::vector< std::vector<int> > counts) {
  for(uint i = 0; i < counts.size(); ++i)
    for(uint j = 0; j < counts[i].size(); ++j)
      this->counts[i][j] += counts[i][j];
}

static double mu_beta(const double alpha, const double beta) {
  return alpha / (alpha + beta);
}

static double sd_beta(const double alpha, const double beta) {
  return sqrt((alpha * beta) / (pow(alpha + beta, 2) * (alpha + beta + 1)));
}

static double cv_beta(const double alpha, const double beta) {
  return sd_beta(alpha, beta) / mu_beta(alpha, beta);
}

static double compute_prob(const double delta, const double alpha,
			   const double beta) {
  /* Construct the beta distribution. */
  boost::math::beta_distribution<> dist(alpha, beta);
  /* Compute the probability from the CDF of the beta distribution. */
  return 1. - cdf(dist, delta);
}

double ConnectivityMatrix::obj_fn_cv() {
  double max_cv = 0.;
  for(uint i = 0; i < origins.size(); ++i) {
    const double alpha_sum = (double)
      (std::accumulate(counts[i].begin(), counts[i].end(), 0) +
       counts[i].size());
    for(uint j = 0; j < destinations.size(); ++j) {
      const double alpha = (double) counts[i][j] + 1;
      if(((alpha / alpha_sum) >= delta) &&
	 (compute_prob(delta, alpha, alpha_sum - alpha) >= pi)) {
	   max_cv = std::max(max_cv, cv_beta(alpha, alpha_sum - alpha));
      }
    }
  }
  return max_cv;
}

void ConnectivityMatrix::set_obj_fn_cv_args(const double delta,
					    const double pi) {
  this->delta = delta;
  this->pi = pi;
}
