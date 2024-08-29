#include <iostream>
#include <Eigen/Dense>
#include <Eigen/SparseCore>
using namespace std;

typedef Eigen::Matrix<long, Eigen::Dynamic, Eigen::Dynamic> MatrixXi;
typedef Eigen::SparseMatrix<double> SpMat; 


//  mss_polynomial_span
//
//  Given:
//
//  coeffs  Sparse double  p-by-m matrix of coefficients.
//                                        no zero columns.
//  degree  Dense  integer p-by-n matrix, each row a vector degree.
//                                        unique rows.
//
//  Describing a set of polynomials p_i(x) = sum_{j=1}^m coeffs(i,j) x^degree(j,:)
//
//  Returns:
//
//  bcoeffs Sparse  double p-by-r matrix of coefficients
//                                       no zero columns.
//
//
//  q_i(x) = sum_{j=1}^m' coeffs(i,j) x^degree(j,:)
//
int linear_span(const SpMatd &coeffs, const MatrixXi &degree)
{
  //  TODO:  Speed up.
  //  Treat polys as graph with an edge for shared monomials.
  //  Run algorithm below indep. on the connected components.

  int m;
  SpMatd *bcoeffs;

  m = 0;
  bcoeffs = new SpMatd(coeffs.rows(),coeffs.cols());

  for(int i = 0; i < coeffs.rows(); i++){
    // Ax = b

    // support == non-zero columns of bcoeffs up to row m.
    // b  should be i-th row of coeffs, restricted to supp.
    // A  should be first m rows of bcoeffs, restricted to supp.

    // Identify the monomials in bcoeffs so var.



  }
}



int main()
{
  cout << "Hello World\n";
  return 0;
}
