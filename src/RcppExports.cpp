// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// runSAX
RObject runSAX(std::vector<float> orgData, int segmentSize, int alphabetSize, bool iSAX);
RcppExport SEXP rsax_runSAX(SEXP orgDataSEXP, SEXP segmentSizeSEXP, SEXP alphabetSizeSEXP, SEXP iSAXSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< std::vector<float> >::type orgData(orgDataSEXP);
    Rcpp::traits::input_parameter< int >::type segmentSize(segmentSizeSEXP);
    Rcpp::traits::input_parameter< int >::type alphabetSize(alphabetSizeSEXP);
    Rcpp::traits::input_parameter< bool >::type iSAX(iSAXSEXP);
    __result = Rcpp::wrap(runSAX(orgData, segmentSize, alphabetSize, iSAX));
    return __result;
END_RCPP
}