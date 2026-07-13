#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include <iostream>

#include "ecos.h"
#include "glblopts.h"

namespace {
  idxint n = 3;
  idxint m = 3;
  idxint p = 2;
  idxint l = 0;

  idxint ncones = 1;
  std::vector<idxint> q(ncones);

  idxint nex = 0;

  std::vector<pfloat> c = {0.0, 0.0, 1.0};
  std::vector<pfloat> h = {0.0, 0.0, 0.0};
  std::vector<pfloat> b = {2.0, 0.0};

  Eigen::MatrixXd A(2, 3);

  Eigen::MatrixXd G(3, 3);

} // namespace

int main() {

  q[0] = 3;

  A << 1.0, 1.0, 0.0,
       -1.0, 1.0, 0.0;
  static Eigen::SparseMatrix<pfloat, Eigen::ColMajor, idxint> Asp = A.sparseView(1, 1e-6);
  Asp.makeCompressed();

  G <<  0.0, 0.0, -1.0,
        -1.0, 0.0, 0.0, 
        0.0, -1.0, 0.0;
  static Eigen::SparseMatrix<pfloat, Eigen::ColMajor, idxint> Gsp = G.sparseView(1, 1e-6);
  Gsp.makeCompressed();

  pwork* workspace = ECOS_setup(n, m, p, 
        l, ncones, q.data(), nex, 
        Gsp.valuePtr(), Gsp.outerIndexPtr(), Gsp.innerIndexPtr(), 
        Asp.valuePtr(), Asp.outerIndexPtr(), Asp.innerIndexPtr(),
        c.data(), h.data(), b.data());

  if (workspace == nullptr) {
        std::cerr << "Error: Could not set up ECOS workspace." << std::endl;
        return -1;
  }

  idxint exitflag = ECOS_solve(workspace); /*[cite: 1, 2] */
  
  if (exitflag == ECOS_OPTIMAL) { /*[cite: 1] */
      std::cout << "Optimal solution found!" << std::endl;
      std::cout << "x1 = " << workspace->x[0] << std::endl; /*[cite: 1] */
      std::cout << "x2 = " << workspace->x[1] << std::endl; /*[cite: 1] */
      std::cout << "t  = " << workspace->x[2] << std::endl; /*[cite: 1] */
  } else {
      std::cout << "ECOS exited with flag " << exitflag << std::endl;
  }

  ECOS_cleanup(workspace, 5U);

  return 0;
}