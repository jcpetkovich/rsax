#include <iostream>
#include <vector>
#include <math.h>

using namespace std;

float* getTable(int n); // hardcoded table

//normalize data to have mean of 0 and standard deviation of 1 for SAX conversion
void normData(vector<float> *data, vector<float> *norm){
    float mean = 0; //mean of the data
    float StDev = 0;  //standard deviation of the data
    
    //calculating mean
    for(int i = 0; i < data -> size(); i ++){
        mean += data -> at(i);
    }
    mean /= data -> size();
    
    //calculating standard deviation
    for(int i = 0; i < data -> size(); i++){
        StDev += pow((data -> at(i)-mean),2);
    }
    StDev /= data -> size();
    StDev = sqrt(StDev);
    
    //normalizing data
    for(int i = 0; i < data -> size(); i++){
        //data.at(i) = (data.at(i)-mean)/StDev;
        norm -> push_back((data -> at(i)-mean)/StDev);
    }
}

//dimensionality reduction of data into segment sizes of the user's choice
void toPAA(vector<float> *data, vector<float> *PAA, int segSize){
    
    //testing if size of data is divisible by the choosen segment size
    if((data -> size()/segSize)%1 == 0){
        //iterating for each segment of data in Time Series data
        for(int i = 0; i < data -> size()/segSize; i++){
            float sum = 0;
            //averaging for each segment
            for(int m = i*segSize; m < (i+1)*segSize; m++){
                sum += data -> at(m);
            }
            PAA -> push_back(sum/segSize);
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
void toSAX(vector<float> *data, vector<int> *SAXWord, vector<int> *card, int numBkPts){
    float* bkPts = getTable(numBkPts-3); //get table values for number of breakpoints choosen
    
    cout<<"Table: ";
    for(int i = 0; i < numBkPts-1;i++){
        cout<<bkPts[i]<<" ";
    }
    cout<<endl;
    
    //iterating through PAA data
    for(int i = 0; i < data -> size(); i++){
        //if the segment's PAA value is less than or equal to 0, compare to negative breakpoints
        if(data -> at(i) <= 0){
            //iterating through negative breakpoint values
            for(int m = (numBkPts/2)-0.5; m >= 0; m--){
                //when the interval is found set it's SAX value to the interval
                if(data -> at(i)>bkPts[m]){
                    SAXWord -> push_back(m+1);
                    break;
                }
                //if the value is less than the lowest breakpoint
                else if((m==0) && (data -> at(i)< bkPts[0])){
                    SAXWord -> push_back(0);
                }
            }
        }
        //if the segment's PAA value is more than 0, compare to positive breakpoints
        else if(data -> at(i) > 0){
            //iterating through positive breakpoint values
            for(int m = (numBkPts/2)-1; m < numBkPts - 1; m++){
                //when the interval is found set it's SAX value to the interval
                if(data -> at(i)<bkPts[m]){
                    SAXWord -> push_back(m);
                    break;
                }
                //if the value is higher than the highest breakpoint
                else if((m==numBkPts - 2) && (data -> at(i)>bkPts[m])){
                    SAXWord -> push_back(numBkPts - 1);
                }
            }
        }
    }
}



//used for testing
int main(){
    vector<float> orgData;
    //------------------------------- Input data from terminal
    float c;
    for(int i = 0; i < 3431; i ++){
        cin>>c;
        orgData.push_back(c);
    }
    cout<<"Data: ";
    for(int i = 0; i < orgData.size(); i++){
        cout << orgData.at(i)<<" ";
    }
    cout<<endl;
    //------------------------------- Normalize Data
    vector<float> nrmData;
    
    normData(&orgData,&nrmData);
    
    cout<<"Normalized Data: ";
    for(int i = 0; i < nrmData.size();i++)
        cout<<nrmData.at(i)<<" ";
    cout<<endl;
    //------------------------------- Calculate PAA
    vector<float> PAA;
    
    toPAA(&nrmData,&PAA,8);
    
    cout<<"PAA: ";
    for(int i = 0; i< PAA.size();i++)
        cout<<PAA.at(i)<<" ";
    cout<<endl;
    //------------------------------- Calculate SAX Word and set corresponding Cardinality
    vector<int> card;
    vector<int> SAXWordFinal;
    
    toSAX(&PAA,&SAXWordFinal,&card,4);
    for(int i = 0; i < SAXWordFinal.size();i++)
        cout<<SAXWordFinal.at(i) << endl;
    cout<<endl;
    //-------------------------------
    return 0;
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