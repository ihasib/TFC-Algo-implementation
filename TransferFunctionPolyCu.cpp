#include<bits//stdc++.h>
#include<complex.h>
using namespace std;


float _Complex MakeComplex(float a, float b)
{
	return  a+b*I;
}


double _Complex RefractiveIndexCu(float lambda)
{
    float r = 0.4822 + (lambda-1305.26)*(0.5701-0.4822)/(1458.82-1305.26); //only this one
    float im = 7.9111 + (lambda-1305.26)*(7.6694-7.9111)/(1458.82-1305.26); //only this one

    return MakeComplex(r,im);
}

int main()
{
    float lambda,polyThickness,polyRefractiveIndex,**T;
    lambda=1220;
    polyThickness=20000;
    polyRefractiveIndex=1.7;
    int rows=2,colms=2;

    T=new float *[rows];

    for(int i=0;i<rows;i++)
    {
        T[i]=new float [colms];
    }

    //T[0][0]=12.2;
    //cout<<T[0][0];

    float _Complex nc = RefractiveIndexCu(lambda);
    cout<<creal(nc)<<" "<<cimag(nc)<<"\n";
    //--------------------------------------------------
    float pi=3.1416;
    float Phi = 2*pi*polyRefractiveIndex*polyThickness/lambda;
    cout<<"phi="<<Phi<<"\n";
    float r12 = (1.0-polyRefractiveIndex)/(1.0+polyRefractiveIndex);
    float t12 = 1.0 + r12;
    cout<<"r12 = "<<r12<<"\n";
    //TODO
    float _Complex r23 = (polyRefractiveIndex-nc)/(polyRefractiveIndex+nc);//complex add
    cout<<"r23 = "<<creal(r23)<<" "<<cimag(r23)<<"\n";
    float _Complex t23 = 1 + r23;//complex add
    cout<<"t23 = "<<creal(t23)<<" "<<cimag(t23)<<"\n";

    float T12[2][2];
    float _Complex T2[2][2], T23[2][2];

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

    for(int i=0;i<noOfRows;i++)
    {
        for(int j=0;j<noOfCols;j++)
        {
            cout<<"i="<<i<<"j="<<j<<"T = "<<T12[i][j]<<"\n";
        }
    }


    T2[0][0] = cexp(I*Phi);//complex
    T2[1][1] = cexp(-I*Phi);//complex


    for(int i=0;i<noOfRows;i++)
    {
        for(int j=0;j<noOfCols;j++)
        {
            cout<<"i="<<i<<"j="<<j<<"T = "<<creal(T2[i][j])<<" + "<<cimag(T2[i][j])<<"i"<<"\n";
        }
    }

    T23[0][1] = r23*T23[0][1];//complex r23
    T23[1][0] = r23*T23[1][0];//complex r23

    for(int i=0;i<noOfRows;i++)
    {
        for(int j=0;j<noOfCols;j++)
        {
            cout<<"i="<<i<<"j="<<j<<"T = "<<creal(T23[i][j])<<" + "<<cimag(T23[i][j])<<"i"<<"\n";
        }
    }
    //T = T12*T2*T23;//real*complex*complex

    return 0;
}
