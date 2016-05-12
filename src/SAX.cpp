#include <iostream>
#include <vector>
#include <fstream>
#include <map>
#include <math.h>
#include <omp.h>
#include <Rcpp.h>
#include <Rmath.h>

using namespace Rcpp;

std::map<int,float*> tableMap;
static bool cardError = false;

float* getTable(int alphabetSize); //Function returns table

// normalize data to have mean of 0 and standard deviation of 1 for SAX
// conversion
void normData(std::vector<float> *data, std::vector<float> *norm) {
  float mean = 0;  // mean of the data
  float StDev = 0; // standard deviation of the data
  int tile = 4096 * std::pow(2, 2);

  float sumx = 0.0f, sumx2 = 0.0f;
#pragma vector aligned
#pragma omp parallel for default(shared) reduction(+ : sumx,                   \
  sumx2) schedule(guided)
    for (decltype(data->size()) i = 0; i < data->size(); i++) {
      sumx += data->at(i);
      sumx2 += data->at(i) * data->at(i);
    }
    mean = sumx / data->size();
  StDev = sqrtf(sumx2 / data->size() - mean * mean);

  // normalizing data
  int dataSize = data->size();
  norm->reserve(dataSize);
  norm->resize(dataSize);

  float invStDev = 1 / StDev;
#pragma vector aligned
#pragma omp parallel for
  for (int ii = 0; ii < dataSize; ii += tile) {
    for (int i = ii; i < ii + tile && i < dataSize; i++) {
      norm->at(i) = ((data->at(i)) - mean) * invStDev;
    }
  }
  //
  //    for(int i = 0; i < data -> size(); i++){
  //        //data.at(i) = (data.at(i)-mean)/StDev;
  //        norm -> push_back((data -> at(i)-mean)/StDev);
  //    }
}

// dimensionality reduction of data into segment sizes of the user's choice
void toPAA(std::vector<float> *data, std::vector<float> *PAA, int segSize) {
  // iterating for each segment of data in Time Series data

  int dataSize = data->size();
  int perThread = 4096; // tbd
  int tile = segSize * perThread;

  PAA->reserve(ceil((float)dataSize / (float)segSize));
  PAA->resize(ceil((float)dataSize / (float)segSize));

#pragma omp parallel for
  for (int ii = 0; ii < dataSize; ii += tile) {
    for (int i = ii; i < ii + tile && i < dataSize; i += segSize) {
      float sum = 0;
      // averaging for each segment
      int m;
      for (m = i; m < (i + segSize) && m < dataSize; m++) {
        sum += data->at(m);
      }
      // PAA -> push_back(sum/segSize);

      PAA->at(i / segSize) = sum / (m-i);
    }
  }
}

// Converting PAA to SAX

/*
 Note: Interval Numbering
 Interval number will be the value of the index greater than the value

 Example:
 Data: 0.2
 intervals: interval [2] = 0.1, interval [3] = 0.4
 Interval Number: 3
 */
void toSAX(std::vector<float> *data, std::vector<int> *SAXWord,
           std::vector<int> *card, int numBkPts) {
  float *bkPts =
    getTable(numBkPts); // get table values for number of breakpoints chosen

  int dataSize = data->size();

  SAXWord->reserve(dataSize);
  SAXWord->resize(dataSize);
  card->reserve(dataSize);
  card->resize(dataSize);
  int tile = 4096;
  // iterating through PAA data
#pragma omp parallel for
  for (int ii = 0; ii < dataSize; ii += tile) {
    for (int i = ii; i < ii + tile && i < dataSize; i++) {
      int counter = 0;
      for (int m = 0; m < numBkPts - 1; m++) {
        if (data->at(i) > bkPts[m])
          counter++;
      }
      // SAXWord -> push_back(counter);
      // card -> push_back(numBkPts);
      SAXWord->at(i) = counter;
      card->at(i) = numBkPts;
    }
  }
}

// function that returns table based on normal gaussian distribution
float *getTable(int alphabetSize) {

  if(tableMap.find(alphabetSize) == tableMap.end()){
    float *tempTable = new float[alphabetSize - 1];
    double interval = 1 / (double)alphabetSize;
    int tile = 10; //tbd

#pragma omp parallel for
    for (int ii = 1; ii < alphabetSize; ii += tile)
      for (int i = ii; i < ii + tile && i < alphabetSize; i++)
        tempTable[i - 1] = ((float)R::qnorm(i * interval, 0.0, 1.0, 1, 0));
    tableMap.insert (std::pair<int,float*> (alphabetSize,tempTable));
    return tempTable;
  }
  else
    return tableMap.find(alphabetSize) -> second;
}

//function that promotes the data of lower cardinality to the data of higher cardinality
int promoteCard(int cardLow, int dataLow, int cardHigh, int dataHigh, bool suppressWarnings){
  int shiftNum = ceil(log2(cardHigh)) - ceil(log2(cardLow));

  if(suppressWarnings)
    cardError = true;
  if((((int)log2(cardHigh) != ceil(log2(cardHigh)))||((int)log2(cardLow) != ceil(log2(cardLow))))&&(!cardError)){
    cardError = true;
    Rf_warning("One or more of your Cardinalities are not of base 2. This may cause some error in the result.");
  }

  dataLow = dataLow << shiftNum;
  if(dataLow < dataHigh){
    int diff = dataHigh - dataLow;
    for(int i = shiftNum - 1; i >= 0; i --){
      int bitAdd = i*2;
      if(bitAdd <= diff){
        diff -= bitAdd;
        dataLow += bitAdd;
      }
      else
        continue;
      if(diff == 0)
        break;
    }
  }
  return dataLow;
}

// Call this from R to convert input data to a SAX Word
// [[Rcpp::export]]
RObject runSAX(std::vector<float> orgData, int segmentSize, int alphabetSize,
               bool iSAX = false) {
  //normalize data
  std::vector<float> nrmData;
  normData(&orgData, &nrmData);

  //calculate PAA for data
  std::vector<float> PAA;

  toPAA(&nrmData, &PAA, segmentSize);

  //assign SAX values based on PAA value
  std::vector<int> card;
  std::vector<int> SAXWordFinal;

  toSAX(&PAA, &SAXWordFinal, &card, alphabetSize);
  //return values
  if (iSAX)
    return List::create(_["SAXWord"] = SAXWordFinal, _["Cardinality"] = card);
  else
    return wrap(SAXWordFinal);
}

//Euclidean distance calculated using raw values
// [[Rcpp::export]]
float eucDis(std::vector<float> rawData1, std::vector<float> rawData2){
  float distance = 0;
  if(rawData1.size()==rawData2.size()){
    int dataSize = rawData1.size();

#pragma omp parallel for reduction(+:distance)
    for(int i = 0; i < dataSize; i++){
      distance += pow((rawData1.at(i)-rawData2.at(i)), 2);
    }
    distance = sqrt(distance);
  }
  else
    std::cerr<<"Error the length of the SAX words are not the same"<<std::endl;

  return distance;
}

//Minimum distance calculated using SAX Words
// [[Rcpp::export]]
float minDis(std:: vector<int> SAXData1, std::vector<int> card1, std::vector<int> SAXData2, std::vector<int> card2, int rawDataSize,bool suppressWarnings = false){
  float distance = 0;
  if((SAXData1.size()==SAXData2.size())&&(card1.size() == card2.size()) && (SAXData1.size() == card1.size())){
    float* bkPts;
    int dataSize = SAXData1.size();
    for(int i = 0; i < dataSize; i++){
      if(card1.at(i) > card2.at(i)){
        SAXData2.at(i) = promoteCard(card2.at(i),SAXData2.at(i), card1.at(i), SAXData1.at(i),suppressWarnings);
        card2.at(i) = card1.at(i);
      }
      else if(card1.at(i) < card2.at(i)){
        SAXData1.at(i) = promoteCard(card1.at(i), SAXData1.at(i), card2.at(i), SAXData2.at(i),suppressWarnings);
        card1.at(i) = card2.at(i);
      }
      bkPts = getTable(card1.at(i));
      if(std::abs(SAXData1.at(i) - SAXData2.at(i)) > 1){
        distance += pow(bkPts[(int)fmax(SAXData1.at(i), SAXData2.at(i))-1] - bkPts[(int)fmin(SAXData1.at(i), SAXData2.at(i))],2);
      }
    }
    distance = std::sqrt(distance);
    distance *= std::sqrt(rawDataSize/SAXData1.size());
  }
  return distance;
}

// This function is used to test the outputs of the normData function
// [[Rcpp::export]]
RObject runNormData(std::vector<float> data) {
  std::vector<float> norm;
  normData(&data, &norm);
  return wrap(norm);
}

// This function is used to test the outputs of the to PAA function
// [[Rcpp::export]]
RObject runToPAA(std::vector<float> data, int segSize) {
  std::vector<float> PAA;
  toPAA(&data, &PAA, segSize);
  return wrap(PAA);
}

// This function is used to test the outputs of the to SAX function
// [[Rcpp::export]]
RObject runToSAX(std::vector<float> data, int brkPtNum) {
  std::vector<int> SAXWord;
  std::vector<int> card;
  toSAX(&data, &SAXWord, &card, brkPtNum);
  return wrap(SAXWord);
}
