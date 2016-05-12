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
// eucDis
float eucDis(std::vector<float> rawData1, std::vector<float> rawData2);
RcppExport SEXP rsax_eucDis(SEXP rawData1SEXP, SEXP rawData2SEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< std::vector<float> >::type rawData1(rawData1SEXP);
    Rcpp::traits::input_parameter< std::vector<float> >::type rawData2(rawData2SEXP);
    __result = Rcpp::wrap(eucDis(rawData1, rawData2));
    return __result;
END_RCPP
}
// minDis
float minDis(std:: vector<int> SAXData1, std::vector<int> card1, std::vector<int> SAXData2, std::vector<int> card2, int rawDataSize, bool suppressWarnings);
RcppExport SEXP rsax_minDis(SEXP SAXData1SEXP, SEXP card1SEXP, SEXP SAXData2SEXP, SEXP card2SEXP, SEXP rawDataSizeSEXP, SEXP suppressWarningsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< std:: vector<int> >::type SAXData1(SAXData1SEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type card1(card1SEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type SAXData2(SAXData2SEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type card2(card2SEXP);
    Rcpp::traits::input_parameter< int >::type rawDataSize(rawDataSizeSEXP);
    Rcpp::traits::input_parameter< bool >::type suppressWarnings(suppressWarningsSEXP);
    __result = Rcpp::wrap(minDis(SAXData1, card1, SAXData2, card2, rawDataSize, suppressWarnings));
    return __result;
END_RCPP
}
// runNormData
RObject runNormData(std::vector<float> data);
RcppExport SEXP rsax_runNormData(SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< std::vector<float> >::type data(dataSEXP);
    __result = Rcpp::wrap(runNormData(data));
    return __result;
END_RCPP
}
// runToPAA
RObject runToPAA(std::vector<float> data, int segSize);
RcppExport SEXP rsax_runToPAA(SEXP dataSEXP, SEXP segSizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< std::vector<float> >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type segSize(segSizeSEXP);
    __result = Rcpp::wrap(runToPAA(data, segSize));
    return __result;
END_RCPP
}
// runToSAX
RObject runToSAX(std::vector<float> data, int brkPtNum);
RcppExport SEXP rsax_runToSAX(SEXP dataSEXP, SEXP brkPtNumSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< std::vector<float> >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type brkPtNum(brkPtNumSEXP);
    __result = Rcpp::wrap(runToSAX(data, brkPtNum));
    return __result;
END_RCPP
}
// testFunc
void testFunc(NumericMatrix mtx);
RcppExport SEXP rsax_testFunc(SEXP mtxSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type mtx(mtxSEXP);
    testFunc(mtx);
    return R_NilValue;
END_RCPP
}
