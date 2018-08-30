#include<bits//stdc++.h>
using namespace std;

void TfcStart(float polyRefractiveIndex, float polyThickness, int spectroResolution);

//#define MAX

int main()
{
    //input section
    float n_poly = 1.7; // This must be inputted by the user in WaferCal % get from thin film data
    float t=20000; // Thickness of polyimide [nanometers] %physical not optical % get from thin film data
    int M=512; // spectrometer resolution %for VITE, its 512

    TfcStart(n_poly,t,M);

    return 0;
}

void TfcStart(float polyRefractiveIndex, float polyThickness, int spectroResolution)
{
    //Apodization matrix
    float *triangular;
    triangular=new float[spectroResolution];
    memset(triangular,0,spectroResolution);

    float value=0, increment=1/((float)spectroResolution/(float)2-1);

    for(int i=0;i<spectroResolution/2;i++)
    {
        triangular[i]=value;
        value+=increment;
    }

    value=1.0;
    for(int i=spectroResolution/2;i<spectroResolution;i++)
    {
        //cout<<i<<" "<<value<<endl;
        triangular[i]=value;
        value-=increment;   //511/  last index contain approximately zero but not actually 0
    }


    float armLength=500000;//ig guess  it from ini file
    //[ X,Y ] =WarpSignalGraph( t, armLength,M,n_poly );

    return;
}
