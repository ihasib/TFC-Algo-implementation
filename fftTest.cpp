#include<bits//stdc++.h>
#include <vector>
#include <complex.h>

using namespace std;

typedef complex<double> cd;
typedef vector<cd> vcd;

vcd fft(const vcd &as) {
  int n = as.size();
  if (n == 1) return vcd(1, as[0]);

  vcd w(n);
  for (int i = 0; i < n; i++) {
    double alpha = 2 * M_PI * i / n;
    w[i] = cd(cos(alpha), sin(alpha));
  }


  vcd A(n / 2), B(n / 2);
  for (int i = 0; i < n / 2; i++) {
    A[i] = as[i * 2];
    B[i] = as[i * 2 + 1];
  }
  vcd Av = fft(A);
  vcd Bv = fft(B);
  vcd res(n);
  for (int i = 0; i < n; i++)
    res[i] =   Av[i % (n / 2)] +
        w[i] * Bv[i % (n / 2)];
  return res;
}

int main()
{
    float ft,a;
    vcd ar;
    int i=0;
    FILE *f=fopen("fftinp.txt","r");
    while(fscanf(f,"%f",&ft)!=EOF)
    {
        //fscanf(f,"%lf",&ft);
        //printf("%f\n",ft);
        //ar[i++]=ft+0i;
        ar.push_back(ft+0i);
    }
    fclose(f);

    vcd vec=fft(ar);

    FILE *fl=fopen("fftout.txt","w");
    for(i=0;i<512;i++)
        fprintf(fl,"%f %fi\n",vec[i].real(),vec[i].imag());

    complex<float> fval;
    fval=3-9i;
    cout<<fval.real()<<" "<<fval.imag();
        //fprintf(fl,"%f %fi\n",creal(vec[i]),cimag(vec[i]));
    fclose(fl);

    return 0;
}
