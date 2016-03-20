#define BOOST_TEST_MODULE SequentialAnalysis

#include <boost/test/unit_test.hpp>
#include "sequential_analysis.hpp"

class ConnectivityMatrixFixture {
public:
  ConnectivityMatrixFixture();
  std::vector<std::string> origins;
  std::vector<std::string> destinations;
  std::vector< std::vector<int> > counts;
  ConnectivityMatrix conn_mat;
};

ConnectivityMatrixFixture::ConnectivityMatrixFixture() {
  /* Set the origins. */
  origins.resize(3);
  origins[0] = std::string("1");
  origins[1] = std::string("2");
  origins[2] = std::string("3");
  /* Set the destinations. */
  destinations.resize(5);
  destinations[0] = std::string("Mid-Maine");
  destinations[1] = std::string("Three States");
  destinations[2] = std::string("Mass Bay");
  destinations[3] = std::string("Nantucket");
  destinations[4] = std::string("Other");
  /* Create the connectivity matrix. */
  conn_mat = ConnectivityMatrix(origins, destinations);
  /* Create a matrix with the counts from the first step of the Gulf of Maine
   * simulation. */
  int counts_mat[15] = {
    5, 3, 15, 2, 142,
    0, 1, 30, 1, 135,
    0, 0, 57, 2, 107
  };
  counts.resize(origins.size());
  for(uint i = 0; i < origins.size(); ++i)
    counts[i].assign(counts_mat + i * destinations.size(),
		     counts_mat + (i + 1) * destinations.size());
}

BOOST_AUTO_TEST_SUITE(connectivity_matrix)
/** Tests that the connectivity matrix has the correct dimensions. */
BOOST_FIXTURE_TEST_CASE(conn_mat_dim, ConnectivityMatrixFixture) {
  const std::vector< std::vector<int> > counts = conn_mat.get_counts();
  BOOST_CHECK_EQUAL(counts.size(), origins.size());
  for(uint i = 0; i < counts.size(); ++i)
    BOOST_CHECK_EQUAL(counts[i].size(), destinations.size());
}

/** Tests that the connectivity matrix is initialized with the correct values.*/
BOOST_FIXTURE_TEST_CASE(conn_mat_counts, ConnectivityMatrixFixture) {
  const std::vector< std::vector<int> > counts = conn_mat.get_counts();
  for(uint i = 0; i < counts.size(); ++i) {
    std::vector<int> zeros(counts[i].size(), 0);
    BOOST_CHECK_EQUAL_COLLECTIONS(counts[i].begin(), counts[i].end(),
				  zeros.begin(), zeros.end());
  }
}

/** Tests that update sets the values for the connectivity matrix correctly. */
BOOST_FIXTURE_TEST_CASE(conn_mat_update, ConnectivityMatrixFixture) {
  /* Call update, then load the new counts. */
  conn_mat.update(counts);
  std::vector< std::vector<int> > new_counts = conn_mat.get_counts();
  /* Check that values. */
  for(uint i = 0; i < counts.size(); ++i) {
    BOOST_CHECK_EQUAL_COLLECTIONS(new_counts[i].begin(), new_counts[i].end(),
				  counts[i].begin(), counts[i].end());
  }
  
}
/** Tests that the objective function returns the correct value. */
BOOST_FIXTURE_TEST_CASE(obj_fn_val, ConnectivityMatrixFixture) {}
/** Tests that set_obj_fn_cv_args sets the correct arguments for the objective 
 * function. */
BOOST_FIXTURE_TEST_CASE(obj_fn_args, ConnectivityMatrixFixture) {}
/** Tests that ConnectivityMatrix::allocate() assigns the correct number of 
 * particles to each origin. */
BOOST_FIXTURE_TEST_CASE(allocate, ConnectivityMatrixFixture) {}
BOOST_AUTO_TEST_SUITE_END()

