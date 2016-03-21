#include <boost/math/distributions/beta.hpp>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
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

ConnectivityMatrix::ConnectivityMatrix(const ConnectivityMatrix & conn_mat) {
  /* Save the origins and destinations. */
  this->origins = conn_mat.get_origins();
  this->destinations = conn_mat.get_destinations();
  /* Set the counts. */
  this->counts.resize(origins.size());
  for(uint i = 0; i < origins.size(); ++i)
    this->counts[i].resize(destinations.size());
  this->update(conn_mat.get_counts());
  /* Set delta and pi. */
  delta = 0.005;
  pi = 0.05;
}

std::vector<uint> ConnectivityMatrix::allocate(const uint n) {
  /* Compute the total number of particles thus far. */
  int total_particles = 0;
  for(uint i = 0; i < origins.size(); ++i)
    total_particles += std::accumulate(counts[i].begin(), counts[i].end(), 0);
  return total_particles == 0 ? this->allocate_uniform(n) :
    this->allocate_optimized(n);
}

std::vector<uint> ConnectivityMatrix::allocate_uniform(const uint n) {
  /* Create a vector to store the allocation. */
  std::vector<uint> allocation(origins.size(), 0);
  /* Assign the rounded down number of particles to each origin. */
  allocation.assign(allocation.size(), n / allocation.size());
  /* Fill in the remaining particles one at a time. */
  uint total_assigned = (n / allocation.size()) * allocation.size();
  uint i = 0;
  while(total_assigned < n) {
    ++allocation[i++];
    ++total_assigned;
  }
  /* Return the allocation. */
  return allocation;
}

std::vector<uint> ConnectivityMatrix::allocate_optimized(const uint n) {
  /* Create a vector to store the allocation. */
  std::vector<uint> allocation(origins.size(), 0);
  /* Create vectors to store the prior costs and the expected posterior costs.*/
  std::vector<double> prior_costs(origins.size(), 0.);
  std::vector<double> exp_costs(origins.size(), 0.);
  for(uint i = 0; i < prior_costs.size(); ++i) {
    prior_costs[i] = this->obj_fn_cv(i);
    exp_costs[i] = this->expected_obj_fn_cv(i, 1);
  }
  /* Iterate through the particles. */
  int release_site = -1;
  for(uint i = 0; i < n; ++i) {
    /* Update the costs if necessary. */
    if(release_site != -1) {
      prior_costs[release_site] = exp_costs[release_site];
      exp_costs[release_site] = this->expected_obj_fn_cv(release_site,
        allocation[release_site] + 1);
    }
    /* Compute the larger of the prior or expected cost for each site. This 
     * corrects for objective functions that do not monotonically increase. */
    std::vector<double> max_costs(exp_costs.size(), 0.);
    for(uint j = 0; j < max_costs.size(); ++j)
      max_costs[j] = std::max(prior_costs[j], exp_costs[j]);
    /* Compute the expected cost after releasing a particle from each site. */
    std::vector<double> release_costs(exp_costs.size(), 0.);
    for(uint j = 0; j < release_costs.size(); ++j) {
      /* Use the expected cost for the release site. */
      double exp_cost = exp_costs[j];
      std::swap(max_costs[j], exp_cost);
      /* Compute the expected posterior costs. */
      release_costs[j] = *std::max_element(max_costs.begin(), max_costs.end());
      /* Restore the max cost. */
      std::swap(max_costs[j], exp_cost);
    }
    const double min_cost = *std::min_element(release_costs.begin(),
					      release_costs.end());
    uint n_matches = 0;
    for(uint j = 0; j < release_costs.size(); ++j) {
      if(release_costs[j] == min_cost) {
	release_site = j;
	++n_matches;
      }
    }
    /* If there is more than 1 match, then we need to make sure that the 
     * particle is released from one of the sites with the maximum 
     * prior cost. */
    if(n_matches > 1) {
      release_site = std::distance(prior_costs.begin(),
        std::max_element(prior_costs.begin(), prior_costs.end()));
    }
    /* Release the particle from the site that minimizes the expected cost. */
    ++allocation[release_site];
  }
  /* Return the allocation. */
  return allocation;  
}

void ConnectivityMatrix::update(const std::vector< std::vector<uint> > counts) {
  for(uint i = 0; i < counts.size(); ++i)
    for(uint j = 0; j < counts[i].size(); ++j)
      this->counts[i][j] += counts[i][j];
}

double ConnectivityMatrix::expected_obj_fn_cv(const uint i, const uint n) {
  /* Create a random number generator. Use the Mersenne-Twister (default in R).
   */
  gsl_rng * rng = gsl_rng_alloc(gsl_rng_mt19937);
  const uint n_reps = 250;
  /* Compute the C.V. for each p_{ij}. */
  std::vector<double> costs(n_reps, 0.);
  std::vector<double> p_sum(5, 0);
  /* Convert the counts to the Dirichlet parameters. */
  std::vector<double> alpha(counts[i].size(), 0.);
  for(uint j = 0; j < counts[i].size(); ++j)
    alpha[j] = (double) counts[i][j] + 1;
  for(uint r = 0; r < n_reps; ++r) {
    /* Create a connectivity matrix using the original counts. */
    ConnectivityMatrix conn_mat(*this);
    /* Generate a probability vector by drawing from a Dirichlet distribution.*/
    std::vector<double> p(destinations.size(), 0.);
    gsl_ran_dirichlet(rng, destinations.size(), alpha.data(), p.data());
    /* Generate a random sample for where these particles go. */
    std::vector< std::vector<uint> >new_counts(origins.size(),
      std::vector<uint>(destinations.size(), 0));
    gsl_ran_multinomial(rng, destinations.size(), n, p.data(),
    			new_counts[i].data());
    /* Update the counts. */
    conn_mat.update(new_counts);
    /* Compute the cost with the updated counts. */
    costs[r] = conn_mat.obj_fn_cv(i);
  }
  /* Cleanup the random number generator. */
  gsl_rng_free(rng);
  return std::accumulate(costs.begin(), costs.end(), 0.) / (double) n_reps;
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
  for(uint i = 0; i < origins.size(); ++i)
    max_cv = std::max(max_cv, this->obj_fn_cv(i));
  return max_cv;
}

double ConnectivityMatrix::obj_fn_cv(const uint i) {
  double max_cv = 0.;
  const double alpha_sum = (double)
    (std::accumulate(counts[i].begin(), counts[i].end(), 0) + counts[i].size());
  for(uint j = 0; j < destinations.size(); ++j) {
    const double alpha = (double) counts[i][j] + 1;
    if(((alpha / alpha_sum) >= delta) &&
       (compute_prob(delta, alpha, alpha_sum - alpha) >= pi)) {
      max_cv = std::max(max_cv, cv_beta(alpha, alpha_sum - alpha));
    }
  }
  return max_cv;
}

void ConnectivityMatrix::set_obj_fn_cv_args(const double delta,
					    const double pi) {
  this->delta = delta;
  this->pi = pi;
}
