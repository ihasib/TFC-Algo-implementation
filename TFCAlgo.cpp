#include<bits//stdc++.h>
#include<complex.h>
using namespace std;

void TfcStart(float polyRefractiveIndex, float polyThickness, int spectroResolution, float lambdaMax, float lambdaMin);
void WarpSignalGraph( float polyThickness,float armLength,int spectroResolution, float polyRefractiveIndex,float *X, float *Y, float lambdaMin, float lambdaMax );
float WarpSignal( float lambda,float armLength,float polyThickness,float polyRefractiveIndex);
float _Complex** TransferFunctionPolyCu( float lambda,float polyThickness,float polyRefractiveIndex);
float _Complex** matrixMultiplication(float _Complex** mat1, float _Complex** mat2, int noOfRowsMat1, int noOfColsMat1, int noOfRowsMat2, int noOfColsMat2);
float _Complex RefractiveIndexCu(float lambda);

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
        X[i] = 1/(kmax-(i+1-1)*deltak); // K space X axis //since our's is 0 indexed
        Y[i] = WarpSignal( X[i],armLength,polyThickness,polyRefractiveIndex); // Intensity
        cout<<"X["<<i<<"] = "<<X[i]<<"\tY["<<i<<"] = "<<Y[i]<<"\n";
    }
    //cout<<kmax<<" "<<kmin<<" "<<deltak;
}



float WarpSignal( float lambda,float armLength,float polyThickness,float polyRefractiveIndex)
{
    float _Complex **T,Rtemp;
    float pi=3.1416;
    T=TransferFunctionPolyCu( lambda,polyThickness,polyRefractiveIndex );//lambda=X[i]
    Rtemp = T[1][0]/T[0][0]; // Reflected R = T_out/T_In

    float TwoPhi=2*2*pi*armLength/lambda;//phase
    float _Complex TempSignal = 0.5*(cexp(1i*TwoPhi)+Rtemp); // electrical field

    //float _Complex Signal = TempSignal*conj(TempSignal); // Intensity abs(E)^2;
    float  Signal = creal( TempSignal*conj(TempSignal) ); // Intensity abs(E)^2;
    /*cout<<"\n\nRtemp = "<<creal(Rtemp)<<" "<<cimag(Rtemp)<<"\n";
    cout<<TwoPhi<<"\n";
    cout<<"\n\nT signal = "<<creal(TempSignal)<<" "<<cimag(TempSignal)<<"\n";
    //cout<<"\n\nsignal = "<<creal(Signal)<<" "<<cimag(Signal)<<"\n";
    cout<<"\n\nsignal = "<<Signal<<"\n";
    */
    return Signal;
}



// T is 2x2 matrix
float _Complex** TransferFunctionPolyCu( float lambda,float polyThickness,float polyRefractiveIndex)
{
    int rows=2,colms=2;

    float _Complex nc = RefractiveIndexCu(lambda);

    //--------------------------------------------------
    float pi=3.1416;
    float Phi = 2*pi*polyRefractiveIndex*polyThickness/lambda;
    float r12 = (1.0-polyRefractiveIndex)/(1.0+polyRefractiveIndex);
    float t12 = 1.0 + r12;

    /*cout<<"Lambda="<<lambda<<"\n";
    cout<<"phi="<<Phi<<"\n";
    cout<<"r12 = "<<r12<<"\n";*/
    //TODO
    float _Complex r23 = (polyRefractiveIndex-nc)/(polyRefractiveIndex+nc);//complex add
    float _Complex t23 = 1 + r23;//complex add

    //array declarartion  for matrix multiplication
    float _Complex **T2, **T23, **T,**T12;

    // array size definition
    T12=new float _Complex*[rows];
    T2=new float _Complex*[rows];
    T23=new float _Complex*[rows];

    for(int i=0;i<rows;i++)
    {
        T12[i]=new float _Complex[colms];
        T2[i]=new float _Complex[colms];
        T23[i]=new float _Complex[colms];
    }

    int noOfRows=2,noOfCols=2;

    for(int i=0;i<noOfRows;i++)
    {
        for(int j=0;j<noOfCols;j++)
        {
            T12[i][j]=1/t12;
            T2[i][j]=0;
            T23[i][j]=1/t23;
        }
    }

    T12[0][1] = r12*T12[0][1];//real
    T12[1][0] = r12*T12[1][0];//real

    T2[0][0] = cexp(I*Phi);//complex
    T2[1][1] = cexp(-I*Phi);//complex

    T23[0][1] = r23*T23[0][1];//complex r23
    T23[1][0] = r23*T23[1][0];//complex r23

    int noOfRowsMat1=2,noOfColsMat1=2,noOfRowsMat2=2,noOfColsMat2=2;
    T=matrixMultiplication(T12,T2,noOfRowsMat1,noOfColsMat1,noOfRowsMat2,noOfColsMat2);
    T=matrixMultiplication(T,T23,noOfRowsMat1,noOfColsMat1,noOfRowsMat2,noOfColsMat2);


    /*cout<<"\n\n T\n";
    for(int i=0;i<noOfRows;i++)
    {
        for(int j=0;j<noOfCols;j++)
        {
            cout<<creal(T[i][j])<<" + "<<cimag(T[i][j])<<"i"<<"\t";
        }
        cout<<"\n";
    }*/
    return T;
}

//input: takes matrix 1 and 2 and their number of rows and columns
//output: multiplied matrix
float _Complex** matrixMultiplication(float _Complex** mat1, float _Complex** mat2, int noOfRowsMat1, int noOfColsMat1, int noOfRowsMat2, int noOfColsMat2)
{

    float _Complex** multiMatrix;

    multiMatrix=new float _Complex *[noOfRowsMat1];

    for(int i=0;i<noOfRowsMat1;i++)
    {
        multiMatrix[i]=new float _Complex [noOfColsMat2];
    }

    //reset multiMatrix

    for(int i=0;i<noOfRowsMat1;i++)
    {
        for(int j=0;j<noOfColsMat2;j++)
        {
            multiMatrix[i][j]=0+0i;
        }
    }

    //matrix multiplication in operation
    for(int row=0;row<noOfRowsMat1;row++)
    {
        for(int col=0;col<noOfColsMat2;col++)
        {
            for(int i=0;i<noOfColsMat1;i++)
            {
                multiMatrix[row][col] += mat1[row][i] * mat2[i][col];
            }
        }
    }

    return multiMatrix;
}



float _Complex RefractiveIndexCu(float lambda)
{
    float r = 0.4822 + (lambda-1305.26)*(0.5701-0.4822)/(1458.82-1305.26); //only this one
    float im = 7.9111 + (lambda-1305.26)*(7.6694-7.9111)/(1458.82-1305.26); //only this one

    return r+(im)*I;
}
