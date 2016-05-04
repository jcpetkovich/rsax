#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//


void timesTwo(NumericVector *x) {
  vector<float> test;
  ofstream outputFile("output.csv");
  outputFile << "s1" << endl;
  for(int i = 0; i < x -> size();i++){
      outputFile << x -> at(i)*2 <<endl;
  }
  outputFile.close();
}
void inputData(NumericVector x){
  timesTwo(&x);
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
timesTwo(42)
*/
