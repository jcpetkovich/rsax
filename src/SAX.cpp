#include <iostream>
#include <vector>
#include <fstream>
#include <math.h>
#include <omp.h>
#include <Rcpp.h>
#include <Rmath.h>

using namespace Rcpp;

float* getTable(int n); // hardcoded table

//normalize data to have mean of 0 and standard deviation of 1 for SAX conversion
void normData(std::vector<float> *data, std::vector<float> *norm){
  float mean = 0; //mean of the data
  float StDev = 0;  //standard deviation of the data
  //calculating mean
#pragma omp parallel for reduction(+:mean)
  for(int i = 0; i < data -> size(); i ++){
    mean += data -> at(i);
  }
  mean /= data -> size();
  //std:: cout<<"Mean: "<<mean<<std::endl;

  //calculating standard deviation
#pragma omp parallel for reduction(+:StDev)
  for(int i = 0; i < data -> size(); i++){
    StDev += std::pow((data -> at(i)-mean),2);
  }
  StDev /= data -> size();
  StDev = std::sqrt(StDev);

  //std::cout<<"Standard Deviation: "<<StDev<<std::endl;

  //normalizing data
  int dataSize = data -> size();
  int tile = std::sqrt(dataSize);
  norm -> reserve(dataSize);
  norm -> resize(dataSize);

#pragma omp parallel for
  for(int ii = 0; ii < dataSize; ii += tile){
    for(int i = ii; i < ii + tile && i < dataSize; i++){
      norm -> at(i) = ((data -> at(i)) - mean)/StDev;
    }
  }
  //
  //    for(int i = 0; i < data -> size(); i++){
  //        //data.at(i) = (data.at(i)-mean)/StDev;
  //        norm -> push_back((data -> at(i)-mean)/StDev);
  //    }
}

//dimensionality reduction of data into segment sizes of the user's choice
void toPAA(std::vector<float> *data, std::vector<float> *PAA, int segSize){
  //iterating for each segment of data in Time Series data

  int dataSize = data -> size();
  int perThread = 10; //tbd
  int tile =  segSize * perThread;

  if(dataSize % segSize == 0){
    PAA -> reserve(dataSize/segSize);
    PAA -> resize(dataSize/segSize);
  }
  else{
    PAA -> reserve((dataSize/segSize)+1);
    PAA -> resize((dataSize/segSize)+1);
  }

# pragma omp parallel for
  for(int ii = 0; ii < dataSize;ii += tile){
    for(int i = ii; i < ii+tile && i < dataSize; i+=segSize){
      float sum = 0;
      //averaging for each segment
      for(int m = i; m < (i + segSize); m++){
        sum += data -> at(m);

      }
      //PAA -> push_back(sum/segSize);
      PAA -> at(i/segSize) = sum/segSize;
    }
  }
}

//Converting PAA to SAX

/*
 Note: Interval Numbering
 Interval number will be the value of the index greater than the value

 Example:
 Data: 0.2
 intervals: interval [2] = 0.1, interval [3] = 0.4
 Interval Number: 3
 */
void toSAX(std::vector<float> *data, std::vector<int> *SAXWord, std::vector<int> *card, int numBkPts){
  float* bkPts = getTable(numBkPts); //get table values for number of breakpoints choosen

  int dataSize = data -> size();

  SAXWord -> reserve (dataSize);
  SAXWord -> resize (dataSize);
  card -> reserve (dataSize);
  card -> resize (dataSize);
  int tile = std::sqrt(dataSize);
  //iterating through PAA data
#pragma omp parallel for
  for(int ii = 0; ii < dataSize;ii += tile){
    for(int i = ii; i < ii+tile && i < dataSize; i++){
      int counter = 0;
      for(int m = 0; m < numBkPts-1; m++){
        if(data -> at(i) > bkPts[m])
          counter++;
      }
      //SAXWord -> push_back(counter);
      //card -> push_back(numBkPts);
      SAXWord -> at(i) = counter;
      card -> at(i) = numBkPts;
    }
  }
}

//function that returns table based on normal gaussian distribution
float* getTable(int alphabetSize) {
  float* table = new float[alphabetSize - 1];

  double interval = 1/(double)alphabetSize;
  int tile = 10;
#pragma omp parallel for
  for(int ii = 1; ii < alphabetSize; ii += tile)
    for (int i = ii; i< ii+tile && i < alphabetSize; i++)
      table[i-1] = ((float)R::qnorm(i*interval, 0.0, 1.0, 1, 0));

  return table;
}

//Call this from R to convert input data to a SAX Word
// [[Rcpp::export]]
RObject runSAX(std::vector<float> orgData, int segmentSize, int alphabetSize, bool iSAX = false){
  //normalize data
  std::vector<float> nrmData;
  normData(&orgData,&nrmData);
  //calculating PAAs
  std:: vector<float> PAA;

  toPAA(&nrmData,&PAA,segmentSize);
  //Comparing and assigning SAX values
  std::vector<int> card;
  std::vector<int> SAXWord;

  toSAX(&PAA,&SAXWord,&card,alphabetSize);

  if(iSAX)
    return List::create(_["SAXWord"] = SAXWord,_["Cardinality"] = card);
  else
    return wrap(SAXWord);
}

//This function is used to test the outputs of the normData function
// [[Rcpp::export]]
RObject runNormData(std::vector<float> data){
  std::vector<float> norm;
  normData(&data,&norm);
  return wrap(norm);
}

//This function is used to test the outputs of the to PAA function
// [[Rcpp::export]]
RObject runToPAA(std::vector<float> data, int segSize){
  std::vector<float> PAA;
  toPAA(&data,&PAA,segSize);
  return wrap(PAA);
}

//This function is used to test the outputs of the to SAX function
// [[Rcpp::export]]
RObject runToSAX(std::vector<float> data,int brkPtNum){
  std::vector<int> SAXWord;
  std::vector<int> card;
  toSAX(&data,&SAXWord,&card,brkPtNum);
  return wrap(SAXWord);
}
