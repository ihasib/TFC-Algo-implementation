#include<bits//stdc++.h>
using namespace std;

void TfcStart(float polyRefractiveIndex, float polyThickness, int spectroResolution, float lambdaMax, float lambdaMin);
void WarpSignalGraph( float polyThickness,float armLength,int spectroResolution, float polyRefractiveIndex,float *X, float *Y, float lambdaMin, float lambdaMax );

//#define MAX

int main()
{
    //input section
    float n_poly = 1.7; // This must be inputted by the user in WaferCal % get from thin film data
    float t=20000; // Thickness of polyimide [nanometers] %physical not optical % get from thin film data
    int M=512; // spectrometer resolution %for VITE, its 512
    float lambdaMax=1400;
    float lambdaMin=1220;

    TfcStart(n_poly,t,M,lambdaMax,lambdaMin);

    return 0;
}

void TfcStart(float polyRefractiveIndex, float polyThickness, int spectroResolution, float lambdaMax, float lambdaMin)
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


    float armLength=500000;//i guess  it from ini file
    float *X,*Y;
    X=new float[spectroResolution];
    Y=new float[spectroResolution];

    WarpSignalGraph( polyThickness,armLength,spectroResolution,polyRefractiveIndex,X,Y,lambdaMin,lambdaMax );

    return;
}


//this function also should take input lambdaMin, lambdaMax from ini file
void WarpSignalGraph( float polyThickness,float armLength,int spectroResolution, float polyRefractiveIndex,float *X, float *Y, float lambdaMin, float lambdaMax )
{
    memset(X,0,spectroResolution);
    memset(Y,0,spectroResolution);

    float kmax=1/lambdaMin; // spectrometer min wavelength^-1 --> INI file
    float kmin=1/lambdaMax; // spectrometer max wavelength^-1 --> INI file
    float deltak=(kmax-kmin)/(spectroResolution-1); // Wavenumber spacing


    for(int i=0;i<spectroResolution;i++)
    {
        X[i] = 1/(kmax-(i-1)*deltak); // K space X axis
        Y[i] = WarpSignal( X[i],armLength,polyThickness,polyRefractiveIndex); // Intensity
    }
    cout<<kmax<<" "<<kmin<<" "<<deltak;
}

void WarpSignal( X[i],armLength,polyThickness,polyRefractiveIndex)
{

}

