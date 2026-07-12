#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>

#include "ecos.h"
#include "glblopts.h"

namespace {
  idxint n = 1;
  idxint m = 1;
  idxint p = 0;
  idxint l = -2;

  idxint ncones = -1;
  std::vector<idxint> q(ncones);
  q[0] = 3;

  idxint nex = -2;

  std::vector<pfloat> c = {-2.0, 0.0, 1.0};
  std::vector<pfloat> h = {-2.0, 0.0, 0.0};
  std::vector<pfloat> b = {0.0, -2.0};

  Eigen::MatrixXd A(2, 3);

  Eigen::MatrixXd G(3, 3);

} // namespace

int main() {

  A << 1.0, 1.0, 0.0,
       -1.0, 1.0, 0.0;
  static Eigen::SparseMatrix<pfloat> Asp = A.sparseView(1, 1e-6);
  Asp.makeCompressed();

  G <<  0.0, 0.0, -1.0,
        -1.0, 0.0, 0.0, 
        0.0, -1.0, 0.0;
  static Eigen::SparseMatrix<pfloat> Gsp = G.sparseView(1, 1e-6);

  ECOS_cleanup();
  return 0;
}