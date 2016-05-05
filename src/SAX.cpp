#include <iostream>
#include <vector>
#include <fstream>
#include <math.h>
#include <Rcpp.h>

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

#pragma omp parallel for
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
    float* bkPts = getTable(numBkPts-3); //get table values for number of breakpoints choosen

    //std::cout<<"Table: ";
    //for(int i = 0; i < numBkPts-1;i++){
    //  std:: cout<<bkPts[i]<<" ";
    //}
    //std::cout<<std::endl;

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


void outData(std::vector<int> *outData){
    std::ofstream outputFile("output.csv");
    outputFile << "s1" << std::endl;
    for(int i = 0; i < outData -> size();i++){
        outputFile << outData -> at(i) <<std::endl;
    }
    outputFile.close();
}

//Call this from R to convert input data to a SAX Word
// [[Rcpp::export]]
RObject runSAX(std::vector<float> orgData, int segmentSize, int alphabetSize, bool iSAX = false){
  //for(int i = 0; i < orgData.size();i++){
   // Rcout<<orgData.at(i)<<std::endl;
  //}
    std::vector<float> nrmData;
    normData(&orgData,&nrmData);
    std:: vector<float> PAA;

    toPAA(&nrmData,&PAA,segmentSize);
    std::vector<int> card;
    std::vector<int> SAXWordFinal;

    toSAX(&PAA,&SAXWordFinal,&card,alphabetSize);
    //outData(&SAXWordFinal);
    if(iSAX)
        return List::create(_["SAXWord"] = SAXWordFinal,_["Cardinality"] = card);
    else
        return wrap(SAXWordFinal);
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

    std::vector<float> norm;
    normData(&data,&norm);

    toPAA(&norm,&PAA,segSize);
    return wrap(PAA);
}

//hardcoded lookup table from 3 to 10 for now.
float* getTable(int n){
    float** bPts = new float*[8];

    bPts[0] = new float[2];
    bPts[0][0] = -0.43;
    bPts[0][1] = 0.43;

    bPts[1] = new float[3];
    bPts[1][0] = -0.67;
    bPts[1][1] = 0;
    bPts[1][2] = 0.67;

    bPts[2] = new float[4];
    bPts[2][0] = -0.84;
    bPts[2][1] = -0.25;
    bPts[2][2] = 0.25;
    bPts[2][3] = 0.84;

    bPts[3] = new float[5];
    bPts[3][0] = -0.97;
    bPts[3][1] = -0.43;
    bPts[3][2] = 0;
    bPts[3][3] = 0.43;
    bPts[3][4] = 0.97;

    bPts[4] = new float[6];
    bPts[4][0] = -1.07;
    bPts[4][1] = -0.57;
    bPts[4][2] = -0.18;
    bPts[4][3] = 0.18;
    bPts[4][4] = 0.57;
    bPts[4][5] = 1.07;

    bPts[5] = new float[7];
    bPts[5][0] = -1.15;
    bPts[5][1] = -0.67;
    bPts[5][2] = -0.32;
    bPts[5][3] = 0;
    bPts[5][4] = 0.32;
    bPts[5][5] = 0.67;
    bPts[5][6] = 1.15;


    bPts[6] = new float[8];
    bPts[6][0] = -1.22;
    bPts[6][1] = -0.76;
    bPts[6][2] = -0.43;
    bPts[6][3] = -0.14;
    bPts[6][4] = 0.14;
    bPts[6][5] = 0.43;
    bPts[6][6] = 0.76;
    bPts[6][7] = 1.22;

    bPts[7] = new float[9];
    bPts[7][0] = -1.28;
    bPts[7][1] = -0.84;
    bPts[7][2] = -0.52;
    bPts[7][3] = -0.25;
    bPts[7][4] = 0;
    bPts[7][5] = 0.25;
    bPts[7][6] = 0.52;
    bPts[7][7] = 0.84;
    bPts[7][8] = 1.28;

    return bPts[n];
}
