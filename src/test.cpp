#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
void testFunc(NumericMatrix mtx){
  int row = mtx.nrow();
  int col = mtx.ncol();
  double *d = new double[row*col];
  std::cout<<"Row: "<<row<<" Col: "<<col<<std::endl;
  for(int c = 0; c < col; c++)
    for(int r = 0; r < row; r++)
      d[c*row + r] = mtx.at(r,c);
  for(int i = 0; i < row*col; i++){
    std::cout<<d[i]<<" ";
  }
}
