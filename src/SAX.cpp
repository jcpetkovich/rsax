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

//Custom Container for omp_lock_t
//Written by: Joel Yliluoma
//http://bisqwit.iki.fi/story/howto/openmp/
#ifdef _OPENMP
struct MutexType{
  MutexType(){omp_init_lock(&lock);}
  ~MutexType(){omp_destroy_lock(&lock);}
  void Lock(){omp_set_lock(&lock);}
  void Unlock(){omp_unset_lock(&lock);}

  MutexType(const MutexType&){omp_init_lock(&lock);}
  MutexType& operator = (const MutexType&){return *this;}
public:
  omp_lock_t lock;
};
#else
/* A dummy mutex that doesn't actually exclude anything,
 * but as there is no parallelism either, no worries. */
struct MutexType{
  void Lock(){}
  void Unlock(){}
};
#endif
//end of custom container

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
void toSAX(std::vector<float> *data, std::vector<int> *SAXWord, std::vector<int> *card, int numBkPts) {
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
RObject runSAX(std::vector<float> orgData, int segmentSize, int alphabetSize, bool iSAX = false) {
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
    int tile = 4096; //tbd
#pragma omp parallel for
    for(int ii = 0; ii < dataSize; ii+=tile){
      for(int i = ii; i < (ii+tile) && i < dataSize; i++){
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
    }
    distance = std::sqrt(distance);
    distance *= std::sqrt(rawDataSize/SAXData1.size());
  }
  return distance;
}

struct element{
  std::vector<int> seq;
  std::vector<int> card;
  float dist = 0;
};

std::vector<element> kMedian(std::vector<int> seq, std::vector<int> card, int step, int windowSize,int clusterNum, int rawDataSize, std::vector<element> centroidIn, bool prevCentroid = false){
  int index = 0;
  if(seq.size() != card.size()){
    Rf_error("Error: The sequence does not match the cardinality size");
    return centroidIn;
  }
  //Assignment of centroids
  if(!prevCentroid){
    centroidIn.clear();
    if(clusterNum > seq.size()-windowSize){
      Rf_error("There are too many clusters and too little groups in the sequence");
      return centroidIn;
    }
    centroidIn.reserve(clusterNum);
    centroidIn.resize(clusterNum);
    for(int i = 0; i < clusterNum; i++){
      element tempCent;
      tempCent.seq.assign(seq.begin()+index,seq.begin()+index+windowSize);
      tempCent.card.assign(card.begin()+index,card.begin()+index+windowSize);
      centroidIn.at(i) = (tempCent);
      index += step;
    }
    if(clusterNum == seq.size()-windowSize)
      return centroidIn;
  }
  else if((centroidIn.at(0).seq.size() != windowSize)||(centroidIn.at(0).card.size() != windowSize)){
    Rf_error("Error: The median sequence or the median cardinality is not the same size as the window size");
    return centroidIn;
  }
  //End of assignment of centroids

  //locks for open mp when inserting into clusters to prevent race conditions/bad sorting
  MutexType *lock = new MutexType[clusterNum];

  //array of vectors to hold each cluster
  std::vector<element> *cluster = new std::vector<element>[clusterNum];
  //iterate through each possible group of size 'windowSize' and interval 'step' in the sequence
  int tile = 4096*step; //tbd
#pragma omp parallel for
  for(int ii = index; ii <= (seq.size()-windowSize); ii+=tile){
    for(int i = ii; i < (ii+tile) && i <= (seq.size()-windowSize);i+=step){
      //temp element to hold the current element
      element temp;
      temp.seq.assign(seq.begin()+i,seq.begin()+i+windowSize);
      temp.card.assign(card.begin()+i,card.begin()+i+windowSize);

      //vector of distance calculations to find the min between each centroid of clusters
      std::vector<float> dist;

      //iterate through centroids and save distance between the current vector and the centroid
      for(int c = 0; c < clusterNum; c++)
        dist.push_back(minDis(centroidIn.at(c).seq,centroidIn.at(c).card,temp.seq,temp.card,rawDataSize,true));

      //get the iterator for the minimum of the vector
      std::vector<float>::iterator it= std::min_element(dist.begin(), dist.end());
      //set the value of the minimum distance of the vector
      temp.dist = *it;
      //determine which cluster the element belongs in
      int clusterIndex = std::distance(dist.begin(),it);
      bool done = false;

      //iterate though elements in each cluster and insert the new element to have the cluster orded from smallest to largest
      for(int v = 0; v < cluster[clusterIndex].size(); v++){
        if(cluster[clusterIndex].at(v).dist > temp.dist){
          lock[clusterIndex].Lock();
          cluster[clusterIndex].insert(cluster[clusterIndex].begin()+v,temp);//omp lock this
          lock[clusterIndex].Unlock();
          done = true;
          break;
        }
      }

      //if every element in the cluster is smaller set the new element to the last index in the cluster
      if(!done){
        lock[clusterIndex].Lock();
        cluster[clusterIndex].push_back(temp);
        lock[clusterIndex].Unlock();
      }
    }
  }

  for(int i = 0; i < clusterNum; i++){
    if(cluster[i].size() != 0)
      centroidIn.at(i) = cluster[i].at((int)((cluster[i].size()-1)/2));
    else
      return kMedian(seq,card,step,windowSize,clusterNum,rawDataSize,centroidIn,false);
  }
  return centroidIn;
}

// [[Rcpp::export]]
List runKMedian(NumericMatrix seq, NumericMatrix card, int step, int windowSize, int clusterNum, int rawDataSize, NumericMatrix centroidSeq, NumericMatrix centroidCard, bool prevCentroid = false){
  std::vector<element> centroid;

  //if the user wants to use inputted centroids
  if(prevCentroid){
    //test if the centroid dimensions are valid
    if((centroidSeq.nrow() != centroidCard.nrow())||(centroidSeq.ncol() != centroidCard.ncol())||(centroidSeq.ncol() != windowSize)){
      Rf_error("Error: Your input Centroid sequence or cardinality don't have the correct dimensions");
      return List::create();
    }
    //put inputted centroids in centroid vector
    for(int i = 0; i < centroidSeq.nrow(); i++){
      NumericVector tempCentSeq = centroidSeq.row(i);
      NumericVector tempCentCard = centroidCard.row(i);
      element temp;
      temp.seq = as<std::vector<int>>(tempCentSeq);
      temp.card = as<std::vector<int>>(tempCentCard);
      centroid.push_back(temp);
    }
  }
  //iterate through all sequences given and feed into kMedian function
  for(int i = 0; i < seq.nrow(); i++){
    NumericVector tempSeq = seq.row(i);
    NumericVector tempCard = card.row(i);
    centroid = kMedian(as<std::vector<int>>(tempSeq),as<std::vector<int>>(tempCard),step,windowSize,clusterNum, rawDataSize, centroid, prevCentroid);
    prevCentroid = true;
  }

  List resultCentroids(clusterNum);
  for(int i = 0; i < clusterNum; i++)
    resultCentroids.at(i) = List::create(Named("Sequence") = wrap(centroid.at(i).seq), Named("Cardinality") = wrap(centroid.at(i).card));

  return resultCentroids;
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
