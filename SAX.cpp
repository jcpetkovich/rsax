#include <iostream>
#include <vector>
#include <math.h>

using namespace std;

//normalize data to have mean of 0 and standard deviation of 1 for SAX conversion
vector<float> normData(vector<float> data){
    float mean = 0; //mean of the data
    float StDev = 0;  //standard deviation of the data
    
    //calculating mean
    for(int i = 0; i < data.size(); i ++){
        mean += data.at(i);
    }
    mean /= data.size();
    
    //calculating standard deviation
    for(int i = 0; i < data.size(); i++){
        StDev += pow((data.at(i)-mean),2);
    }
    StDev /= data.size();
    StDev = sqrt(StDev);
    
    //normalizing data
    for(int i = 0; i < data.size(); i++){
        data.at(i) = (data.at(i)-mean)/StDev;
    }
    return data;
}

//dimensionality reduction of data into segment sizes of the user's choosing
vector<float> calcPAA(vector<float> data, int segSize){
    data = normData(data); //normalize data
    vector<float> PAA; //resulting data
    
    //testing if size of data is divisible by the choosen segment size
    if((data.size()/segSize)%1 == 0){
        //iterating for each segment of data in Time Series data
        for(int i = 0; i < data.size()/segSize; i++){
            float sum = 0;
            //averaging for each segment
            for(int m = i*segSize; m < (i+1)*segSize; m++){
                sum += data.at(m);
            }
            PAA.push_back(sum/segSize);
        }
    }
    return PAA;
}


//used for testing
int main(){
    vector<float> data1;
    float c;
    for(int i = 0; i < 8; i ++){
        cin>>c;
        data1.push_back(c);
    }
    cout<<"Data: ";
    for(int i = 0; i < data1.size(); i++){
        cout << data1.at(i)<<" ";
    }
    cout<<endl;
    vector<float> data = normData(data1);
    vector<float> PAA = calcPAA(data,2);
    cout<<"PAA: ";
    for(int i = 0; i< PAA.size();i++){
        cout<<PAA.at(i)<<" ";
    }
    cout<<endl;
    return 0;
}