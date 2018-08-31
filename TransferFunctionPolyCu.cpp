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
    int rows=2,colms=2;

    T=new float *[rows];

    for(int i=0;i<rows;i++)
    {
        T[i]=new float [colms];
    }

    //T[0][0]=12.2;
    //cout<<T[0][0];

    float _Complex nc = RefractiveIndexCu(lambda);
    cout<<creal(nc)<<" "<<cimag(nc);
    //--------------------------------------------------
    float pi=3.1416;
    float Phi = 2*pi*polyRefractiveIndex*polyThickness/lambda;

    float r12 = (1-polyRefractiveIndex)/(1+polyRefractiveIndex);
    float t12 = 1.0 + r12;

    //TODO
    float _Complex r23 = (np-nc)/(np+nc);//complex add
    float _Complex t23 = 1 + r23;//complex add


    return 0;
}
