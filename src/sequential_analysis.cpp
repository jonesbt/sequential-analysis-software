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
}

void ConnectivityMatrix::update(const std::vector< std::vector<int> > counts) {
  for(uint i = 0; i < counts.size(); ++i)
    for(uint j = 0; j < counts[i].size(); ++j)
      this->counts[i][j] += counts[i][j];
}
