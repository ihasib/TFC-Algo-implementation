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

    return r+(im)*I;
}

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

int main()
{
    float lambda,polyThickness,polyRefractiveIndex;

    lambda=1220;
    polyThickness=20000;
    polyRefractiveIndex=1.7;
    int rows=2,colms=2;

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


    //array declarartion  for matrix multiplication
    float _Complex **T2, **T23, **T,**T12;

    // array size definition
   // T=new float _Complex*[rows];
    T12=new float _Complex*[rows];
    T2=new float _Complex*[rows];
    T23=new float _Complex*[rows];

    for(int i=0;i<rows;i++)
    {
        //T[i]=new float _Complex[colms];
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


    cout<<"\n\n T12\n";
    for(int i=0;i<noOfRows;i++)
    {
        for(int j=0;j<noOfCols;j++)
        {
            cout<<creal(T12[i][j])<<" + "<<cimag(T12[i][j])<<"i"<<"\t";
        }
        cout<<"\n";
    }


    T2[0][0] = cexp(I*Phi);//complex
    T2[1][1] = cexp(-I*Phi);//complex


    cout<<"\n\n T2\n";
    for(int i=0;i<noOfRows;i++)
    {
        for(int j=0;j<noOfCols;j++)
        {
            cout<<creal(T2[i][j])<<" + "<<cimag(T2[i][j])<<"i"<<"\t";
        }
        cout<<"\n";
    }

    T23[0][1] = r23*T23[0][1];//complex r23
    T23[1][0] = r23*T23[1][0];//complex r23


    cout<<"\n\n T23\n";
    for(int i=0;i<noOfRows;i++)
    {
        for(int j=0;j<noOfCols;j++)
        {
            cout<<creal(T23[i][j])<<" + "<<cimag(T23[i][j])<<"i"<<"\t";
        }
        cout<<"\n";
    }

    int noOfRowsMat1=2,noOfColsMat1=2,noOfRowsMat2=2,noOfColsMat2=2;
    T=matrixMultiplication(T12,T2,noOfRowsMat1,noOfColsMat1,noOfRowsMat2,noOfColsMat2);


    cout<<"\n\n T\n";
    for(int i=0;i<noOfRows;i++)
    {
        for(int j=0;j<noOfCols;j++)
        {
            cout<<creal(T[i][j])<<" + "<<cimag(T[i][j])<<"i"<<"\t";
        }
        cout<<"\n";
    }

    T=matrixMultiplication(T,T23,noOfRowsMat1,noOfColsMat1,noOfRowsMat2,noOfColsMat2);


    cout<<"\n\n T\n";
    for(int i=0;i<noOfRows;i++)
    {
        for(int j=0;j<noOfCols;j++)
        {
            cout<<creal(T[i][j])<<" + "<<cimag(T[i][j])<<"i"<<"\t";
        }
        cout<<"\n";
    }
    return 0;
}
