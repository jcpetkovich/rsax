#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <Rcpp.h>
using namespace Rcpp;
#include "zlib.h"

//size of SAX series data
uLongf cData;

struct window {
  int head;
  int size;
};

//calculates the Compression-based Dissimilarity Measure (CDM) between two strings dataW and dataD
double CDM(Bytef *dataW, Bytef *dataDNW, int windowSize, int dataSize){
  uLongf cWindow = compressBound(windowSize);
  Bytef *destW = new Bytef[cWindow];
  compress(destW, &cWindow,dataW, windowSize);

  uLongf cDNW = compressBound(windowSize+dataSize);
  Bytef *destDNW = new Bytef[cDNW];
  compress(destDNW, &cDNW,dataDNW, windowSize+dataSize);

  return (double)cDNW/((double)cData + (double)cWindow);
}

//concatenates of dataD and dataD
Bytef* cat(Bytef* dataD, Bytef * dataW, int dataSize, int windowSize){
  Bytef *dataDNW = new Bytef[windowSize+dataSize];
  for(int i = 0; i < dataSize; i++){
    if(i < windowSize)
      dataDNW[i+dataSize] = dataW[i];
    dataDNW[i] = dataD[i];
  }
  return dataDNW;
}

//' Compute indeces of anomalies in SAX series based on Window Comparison Anomaly Detection Algorithm
//' @title Window Comparison Anomaly Detection (WCAD)
//' @param data A vector containing the SAX series to be tested for anomalies
//' @param windowNum The inital number of windows to start off with (The number of pieces we split the original SAX series into)
//' @param cdmThreshold The minimum a CDM value can be to be considered "interesting" or containing an anomaly
//' @return An unordered NumericVector containing the indices of the detected anomalies
// [[Rcpp::export]]
NumericVector WCAD(std::vector<int> data, int windowNum, float cdmThreshold){
  //Store windows marked as "interesting"
  std::vector<window> winChosen;
  //store anomaly indices
  std::vector<int> anomalies;

  int dataSize = data.size();
  int windowSize = dataSize/windowNum;
  //Saving SAX series data into Bytef form, to make compatiable with zLib
  Bytef *dataD = new Bytef[dataSize];
  std::transform(data.begin(), data.end(), dataD, [](int i){return '0' + i;});
  //Compute compression of SAX series data
  cData = compressBound(dataSize);
  Bytef *destD = new Bytef[cData];
  compress(destD, &cData, dataD, dataSize);

  //iterate through all windows and compute the CDM between each window and the SAX series data
  for(int w = 0; w < windowNum; w++){
    Bytef *dataW = new Bytef[windowSize];
    std::transform(data.begin()+windowSize*w, data.begin()+windowSize*(w+1), dataW, [](int i){return '0' + i;});

    Bytef *dataDNW = cat(dataD,dataW,dataSize,windowSize);
    //Any window with a CDM larger than the threshold will be stored and considered as "interesting"
    if(CDM(dataW,dataDNW,windowSize,dataSize)>=cdmThreshold){
      //window temp;
      //temp.head = windowSize*w;
      //temp.size = windowSize;
      //winChosen.push_back(temp);
      anomalies.push_back(w);
    }
    //Rcout<<CDM(dataW,dataDNW,windowSize,dataSize)<<" ";
    if((w+1)%10 == 0)
      Rcout<<std::endl;
  }
  return wrap(anomalies);
}
