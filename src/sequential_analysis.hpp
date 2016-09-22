#ifndef __INCLUDE_SEQUENTIAL_ANALYSIS_HPP__
#define __INCLUDE_SEQUENTIAL_ANALYSIS_HPP__

#include <string>
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
   * \param block_size The number of particles to allocate as a block.
   * \return A vector giving the number of particles to release from each origin
   * site.
   */
  std::vector<int> allocate(const int n, const int block_size=1);

  /** Returns the matrix of counts. Each element count[i][j] is equal to the 
   * total number of particles that were released from origin i and arrived at
   * destination j. */
  inline std::vector< std::vector<int> > get_counts() const {
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
  void update(const std::vector< std::vector<int> > counts);
  /** Updates the connectivity matrix to include the newly observed data by 
   * adding the counts to the existing ones. 
   *
   * \param The index of the origin for which the counts are being provided.
   * \param The newly observed counts, where counts[j] gives the number of 
   * particles released from origin i and arriving at destination j.
   */
  void update(const int i, const std::vector<int> counts);
private:
  /** Creates a copy of this connectivity matrix. */
  ConnectivityMatrix(const ConnectivityMatrix &conn_mat);
  /** Uniformly allocates the particles across the origin sites. 
   *
   * \param n The number of particles to allocate.
   */
  std::vector<int> allocate_uniform(const int n);
  /** Allocates the particles across the origin sites, attempting to minimize 
   * the posterior expected value of the objective function.
   *
   * \param n The number of particles to allocate.
   * \param block_size The number of particles to allocate as a block.
   */
  std::vector<int> allocate_optimized(const int n, const int block_size);
  /** Computes the expected value of the objective function if n additional
   * particles were to be allocated from origin i. 
   * 
   * \param i The release site for the new particles.
   * \param n The number of new particles to release.
   */
  double expected_obj_fn_cv(const int i, const int n);
  /** Computes the coefficient of variation based objective function for origin
   * i. See ConnectivityMatrix::obj_fn_cv() for details. 
   *
   * \param i The origin for which the objective function should be computed.
   */
  double obj_fn_cv(const int i);
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
  std::vector< std::vector<int> > counts;
  /** The parameter delta for ConnectivityMatrix::obj_fn_cv(). */
  double delta;
  /** The parameter pi for ConnectivityMatrix::obj_fn_cv(). */
  double pi;
};

#endif
