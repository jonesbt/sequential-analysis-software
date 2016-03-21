#ifndef __INCLUDE_SEQUENTIAL_ANALYSIS_HPP__
#define __INCLUDE_SEQUENTIAL_ANALYSIS_HPP__

#include <string>
#include <sys/types.h>
#include <vector>

class ConnectivityMatrix {
public:
  /** Saves the names of the origins and destinations and initializes the 
   * connectivity matrix. All of the counts are set to 0, and the objective 
   * function values are set to their defaults (delta = 0.005, pi=0.05).
   *
   * \param origins The names of the origin sites.
   * \param destinations The names of the destination sites.
   */
  ConnectivityMatrix(
    const std::vector<std::string> origins =
      std::vector<std::string>(0, std::string("")),
    const std::vector<std::string> destinations =
      std::vector<std::string>(0, std::string("")));

  /** Allocates a batch of n particles to origin sites in order to minimize the
   * expected posterior value of the objective function.
   *
   * \param n The number of particles in this batch to allocate.
   * \return A vector giving the number of particles to release from each origin
   * site.
   */
  std::vector<uint> allocate(const uint n);

  /** Returns the matrix of counts. Each element count[i][j] is equal to the 
   * total number of particles that were released from origin i and arrived at
   * destination j. */
  inline std::vector< std::vector<uint> > get_counts() const {
    return this->counts;
  }
  
  /** Computes the coefficient of variance based objective function as given
   * below. 
   *
   * max_{i,j}(CV_{ij} : Prob(P_{ij} > delta) > pi)
   */
  double obj_fn_cv();
  
  /** Sets the arguments of the coefficient of variance based objective 
   * function. 
   *
   * \param delta The parameter delta for ConnectivityMatrix::obj_fn_cv().
   * \param pi The parameter pi for ConnectivityMatrix::obj_fn_cv().
   * \see ConnectivityMatrix::obj_fn_cv()
   */
  void set_obj_fn_cv_args(const double delta, const double pi);
  
  /** Updates the connectivity matrix to include the newly observed data by 
   * adding the counts to the existing ones. 
   *
   * \param The newly observed counts, where counts[i][j] gives the number of 
   * particles released from origin i and arriving at destination j.
   */
  void update(const std::vector< std::vector<uint> > counts);
private:
  /** Creates a copy of this connectivity matrix. */
  ConnectivityMatrix(const ConnectivityMatrix &conn_mat);
  /** Uniformly allocates the particles across the origin sites. 
   *
   * \param n The number of particles to allocate.
   */
  std::vector<uint> allocate_uniform(const uint n);
  /** Allocates the particles across the origin sites, attempting to minimize 
   * the posterior expected value of the objective function.
   *
   * \param n The number of particles to allocate.
   */
  std::vector<uint> allocate_optimized(const uint n);
  /** Computes the expected value of the objective function if n additional
   * particles were to be allocated from origin i. 
   * 
   * \param i The release site for the new particles.
   * \param n The number of new particles to release.
   */
  double expected_obj_fn_cv(const uint i, const uint n);
  /** Computes the coefficient of variation based objective function for origin
   * i. See ConnectivityMatrix::obj_fn_cv() for details. 
   *
   * \param i The origin for which the objective function should be computed.
   */
  double obj_fn_cv(const uint i);
  /* Returns the origin sites. */
  inline std::vector<std::string> get_origins() const { return origins; }
  /* Returns the destination sites. */
  inline std::vector<std::string> get_destinations() const {
    return destinations;
  }
  /** The origin sites. */
  std::vector<std::string> origins;
  /** The destination sites. */
  std::vector<std::string> destinations;
  /** The counts. counts[i][j] is the number of particles released from origin
   * i that arrived at destination j. */
  std::vector< std::vector<uint> > counts;
  /** The parameter delta for ConnectivityMatrix::obj_fn_cv(). */
  double delta;
  /** The parameter pi for ConnectivityMatrix::obj_fn_cv(). */
  double pi;
};

#endif
