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
  std::vector<uint> allocate(const uint n) {
    return std::vector<uint>(origins.size(), 0);
  }
  
  inline std::vector< std::vector<int> > get_counts() const {
    return this->counts;
  }
  
  /** Computes the coefficient of variance based objective function as given
   * below. 
   *
   * max_{i,j}(CV_{ij} : Prob(P_{ij} > delta) > pi)
   */
  double obj_fn_cv() { return 0.; }
  
  /** Sets the arguments of the coefficient of variance based objective 
   * function. 
   *
   * \param delta The parameter delta for ConnectivityMatrix::obj_fn_cv().
   * \param pi The parameter pi for ConnectivityMatrix::obj_fn_cv().
   * \see ConnectivityMatrix::obj_fn_cv()
   */
  void set_obj_fn_cv_args(const double delta, const double pi) {}
  
  /** Updates the connectivity matrix to include the newly observed data by 
   * adding the counts to the existing ones. 
   *
   * \param The newly observed counts, where counts[i][j] gives the number of 
   * particles released from origin i and arriving at destination j.
   */
  void update(const std::vector< std::vector<int> > counts);
private:
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
