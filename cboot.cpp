#include <math.h>
#include <fstream>
#include <string.h>
#include <random> 

using namespace std;

struct fatique{
 double pv;
 double alpha;
 double an;
 int k;
 double *cx;
 double *y;
 double *w;
 double Q;
};
fatique fat;

 void fatiq(int k,double *w,double *cx,double *y);
 double fatiq_fun(double pvx);
 double bisectionMethod(double (*fct)(double),double a, double b, double tolerance);
 double normal_inv(double p, double mu, double sigma);


int main () {

 int i,k,k1,k2,km,j,n1,n2,p;
 int *m1,*m2;
 long int mm,jj;
 string st,ff;
 double *cx1,*cx2,*lgN1,*lgN2,*a1,*a2,*sd1,*sd2,*w1,*w2,*logb,*cb1,*cb2;
 double s1,s2,cb7,c001,c005,c01,pcrit,alpha; 
 double *nplus,*nminus,*nnull;

 ff="cboot";
 ifstream fopen(ff+".inp");
 ofstream fout(ff+".out");

////////////////Read Data///////////////////////////////////////

 fopen>>st;
 fopen>>mm;
 fopen>>st;
 fopen>>alpha;
 //c001=1.2879;c005=0.98;c01=0.8224;
 //pcrit=((mm-1.)/2.-c005*sqrt(mm+1.));
 pcrit=normal_inv(alpha,0.5*mm,0.5*sqrt(mm));

 fat.cx=new double[20]; 
 fat.w=new double[20]; 
 fat.y=new double[20]; 

 fopen>>st;
 fopen>>k;
 logb = new double[k];
 fopen>>st;
 for(i=0;i<k;i++) fopen>>logb[i];
 cb1=new double[k]; 
 cb2=new double[k]; 
 nplus=new double [k];
 nnull=new double [k];
 nminus=new double [k];
 for (i=0;i<k;i++) {nplus[i]=0;nnull[i]=0;nminus[i]=0;}
 
 fopen>>st;
 fopen>>k1;
 cx1=new double[k1];
 a1=new double[k1];
 sd1=new double[k1];
 w1=new double[k1];

 fopen>>st;
 for(i=0;i<k1;i++) fopen>>cx1[i];
 m1=new int[k1];
 fopen>>st;
 for(i=0;i<k1;i++) fopen>>m1[i]; 

 n1=0;km=0;
 for(i=0;i<k1;i++) n1+=m1[i]; 
 lgN1=new double[n1]; 
 for(i=0;i<k1;i++) {
   fopen>>st; 
   for(j=0;j<m1[i];j++) {
      fopen>>lgN1[j+km];
   }
   km+=m1[i];
 } 

//////////////////////Out Data/////////////////////////////////////////

fout<<"Alloy_1"<<endl;
for(i=0;i<k1;i++) fout<<cx1[i]<<"    ";
fout<<endl;
 km=0;
 for(i=0;i<k1;i++) {
   fout<<"Sample_"<<i+1<<endl;
   for(j=0;j<m1[i];j++) {
     fout<<lgN1[j+km]<<"  ";
   }
   fout<<endl;
   km+=m1[i];
 }

 fopen>>st;
 fopen>>k2;
 cx2=new double[k2];
 a2=new double[k2];
 sd2=new double[k2];
 w2=new double[k2];

 fopen>>st;
 for(i=0;i<k2;i++) fopen>>cx2[i];
 m2=new int[k2];
 fopen>>st;
 for(i=0;i<k2;i++) fopen>>m2[i]; 
 n2=0;km=0;
 for(i=0;i<k2;i++) n2+=m2[i]; 
 lgN2=new double[n2]; 
 for(i=0;i<k2;i++) {
   fopen>>st; 
   for(j=0;j<m2[i];j++)  fopen>>lgN2[j+km];
   km+=m2[i];
 } 
 
 fout<<"Alloy_2"<<endl;
 km=0;
 for(i=0;i<k2;i++) fout<<cx2[i]<<"    ";
 fout<<endl;
 for(i=0;i<k2;i++) {
   fout<<"Sample_"<<i+1<<endl;
   for(j=0;j<m2[i];j++)   fout<<lgN2[j+km]<<"  ";
   fout<<endl;
   km+=m2[i];
 }

 fout<<endl; fout<<endl;       //End Input,Output Data


/////////////////Simulation///////////////////////////////////////////

    random_device rd;
    mt19937 g(rd());

for(jj=0;jj<mm;jj++) {       // start cycle for simulation
 
////////////Alloy 1////////////////////////////

 km=0;
 for(i=0;i<k1;i++) {
    uniform_int_distribution<>dist(0,m1[i]-1);
    s1=0;s2=0;
    for(j=0;j<m1[i];j++) {
       p=dist(g);
       s1+=lgN1[p+km];s2+=pow(lgN1[p+km],2);      //bootstrep
    }
    a1[i]=s1/m1[i];
    sd1[i]=(s2-(a1[i]*a1[i])*m1[i])/(m1[i]-1.);
    w1[i]=m1[i]*a1[i]*a1[i]/sd1[i];
    km+=m1[i];
 }

////////////////////////Fatique Curve _1///////////////////////////////////////
 
   fatiq(k1,w1,cx1,a1); 
   cb7=fat.pv+fat.an*pow(7.,-fat.alpha);
   for(i=0;i<k;i++) cb1[i]=(fat.pv+fat.an*pow(logb[i],-fat.alpha))/cb7;
  
////////////Alloy 2//////////////////////////////////////////
  
 km=0;
 for(i=0;i<k2;i++) {
    uniform_int_distribution<>dist(0,m2[i]-1);
    s1=0;s2=0;
    for(j=0;j<m2[i];j++) {
       p=dist(g);
       s1+=lgN2[p+km]; s2+=pow(lgN2[p+km],2);
    }
    a2[i]=s1/m2[i];
    sd2[i]=(s2-(a2[i]*a2[i])*m2[i])/(m2[i]-1.);
    w2[i]=m2[i]*a2[i]*a2[i]/sd2[i];
    km+=m2[i];
 }

////////////////////////Fatique Curve _2///////////////////////////////////////

   fatiq(k2,w2,cx2,a2); 
   cb7=fat.pv+fat.an*pow(7,-fat.alpha);
   for(i=0;i<k;i++) cb2[i]=(fat.pv+fat.an*pow(logb[i],-fat.alpha))/cb7;
 
/////////////////////Signe Criterion////////////////////////////

    for(i=0;i<k;i++) {
      if(cb1[i]>cb2[i]) nplus[i]=nplus[i]+1;
      if(cb1[i]==cb2[i]) nnull[i]=nnull[i]+1;
      if(cb1[i]<cb2[i])  nminus[i]=nminus[i]+1;
    }

} //end jj

 fout<<"pcrit="<<pcrit<<endl;  
 fout<<"nplus=";    
 for(i=0;i<k;i++) fout<<nplus[i]<<"    ";
 fout<<endl;
 fout<<"nnull=";
 for(i=0;i<k;i++) fout<<nnull[i]<<"    ";
 fout<<endl;
 fout<<"nminus=";
 for(i=0;i<k;i++) fout<<nminus[i]<<"    ";
 fout<<endl;

 for(i=0;i<k;i++) {
    if (fmin(nplus[i],nminus[i])>=pcrit) fout<<"H0+"<<"   ";
    if (fmin(nplus[i],nminus[i])<pcrit) fout<<"H0-"<<"   ";
 }

/////////////////////////////////////////////////////////////////////////////////////
 
 fopen.close();fout.close();
 delete [] cx1,cx2,lgN1,lgN2,m1,m2,a1,a2,sd1,sd2,w1,w2,logb,cb1,cb2,nplus,nminus,nnull;
 return 0;

}

///////////////////////////////////////////////////////////

void fatiq(int k,double *w,double *cx,double *y) {
    int i;
    double eps,s1,s2,s3,s5,xcp,a,b,pvx;
 
    eps=1e-10;

    fat.k=k;
    for (i=0;i<k;i++) {
       fat.w[i]=w[i];
       fat.y[i]=y[i];
       fat.cx[i]=cx[i];
   }

      pvx=bisectionMethod(fatiq_fun,0,cx[k-1],eps);

     s1=0.;s2=0.;s3=0.;
    for (i=0;i<k;i++) {
      s1+=w[i]*log(cx[i]-pvx);
      s2+=w[i];
      s3+=w[i]*log(y[i]);
    }
  xcp=s1/s2;a=s3/s2;s1=0.;s3=0.;
  for (i=0;i<k;i++) {
   s1=s1+w[i]*pow(log(cx[i]-pvx)-xcp,2);
   s3=s3+w[i]*(log(cx[i]-pvx)-xcp)*log(y[i]);
  }
  b=s3/s1;
  fat.pv=pvx;
  fat.alpha=-1./b;
  fat.an=exp(xcp-a/b);
  s5=0;
  for (i=0;i<k;i++) {
   s5=s5+w[i]*pow(log(y[i])-a-b*(log(fat.cx[i]-pvx)-xcp),2);
  }
  fat.Q=s5/s2;
}

///////////////////////////////////////////////////////////////

 double fatiq_fun(double x) {
    int i;
    double s1,s2,s3,s4,s5,xcp,a,b;

   if(x<0 || x>=fat.cx[fat.k-1]) return 1e10;

   s1=0.;s2=0.;s3=0.;
   for (i=0;i<fat.k;i++) {
     s1+=fat.w[i]*log(fat.cx[i]-x);
     s2+=fat.w[i];
     s3+=fat.w[i]*log(fat.y[i]);
   }
  xcp=s1/s2;
  a=s3/s2;
  s1=0.;s4=0.;
  for (i=0;i<fat.k;i++) {
   s1=s1+fat.w[i]*pow(log(fat.cx[i]-x)-xcp,2);
   s4=s4+fat.w[i]*(log(fat.cx[i]-x)-xcp)*log(fat.y[i]);
  }
  b=s4/s1;
  s5=0;
 for (i=0;i<fat.k;i++) {
   s5=s5+fat.w[i]*pow(log(fat.y[i])-a-b*(log(fat.cx[i]-x)-xcp),2);
  }
  return s5/s2;
 }


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double bisectionMethod(double (*fct)(double),double a, double b, double tolerance) {
  // h - шаг для численного дифференцирования
  double h = 0.00001;
 
  while (abs(b - a) > tolerance) {
    double mid = (a + b) / 2;
    double derivative = (fct(mid + h)-fct(mid - h)) / (2 * h);
    if (derivative > 0) {
      b = mid; 
    } else if (derivative < 0) {
      a = mid;
    } else {
      return mid;
    }
  }
  return (a + b) / 2;
}

///////////////////////////////////////////////

double normal_inv(double p, double mu, double sigma) {

    double q,r,num,x,den;
  
    q = p - 0.5;
    if(fabs(q) <= 0.425) {
        r = 0.180625 - q * q;
        num = (((((((2.5090809287301226727e+3 * r +3.3430575583588128105e+4) * r + 6.7265770927008700853e+4) * r + 4.5921953931549871457e+4) * r +
                     1.3731693765509461125e+4) * r +1.9715909503065514427e+3) * r +1.3314166789178437745e+2) * r +3.3871328727963666080e+0) * q;
        den = (((((((5.2264952788528545610e+3 * r +2.8729085735721942674e+4) * r + 3.9307895800092710610e+4) * r + 2.1213794301586595867e+4) * r +
                     5.3941960214247511077e+3) * r + 6.8718700749205790830e+2) * r + 4.2313330701600911252e+1) * r +1.0);
        x = num / den;
        return(mu + (x * sigma));
     }

    if(q <= 0.0)  {
      r=p;
   }
    else {
      r=1.0-p;
    }
    r =sqrt(-log(r));
    if(r <= 5.0) {
        r = r - 1.6;
        num = (((((((7.74545014278341407640e-4 * r + 2.27238449892691845833e-2) * r +2.41780725177450611770e-1) * r +1.27045825245236838258e+0) * r +
                     3.64784832476320460504e+0) * r + 5.76949722146069140550e+0) * r +4.63033784615654529590e+0) * r +1.42343711074968357734e+0);
        den = (((((((1.05075007164441684324e-9 * r +5.47593808499534494600e-4) * r +1.51986665636164571966e-2) * r + 1.48103976427480074590e-1) * r +
                     6.89767334985100004550e-1) * r +1.67638483018380384940e+0) * r +2.05319162663775882187e+0) * r +
                     1.0);
    } else {
        r = r - 5.0;
        num = (((((((2.01033439929228813265e-7 * r + 2.71155556874348757815e-5) * r +1.24266094738807843860e-3) * r + 2.65321895265761230930e-2) * r +
                     2.96560571828504891230e-1) * r +1.78482653991729133580e+0) * r + 5.46378491116411436990e+0) * r + 6.65790464350110377720e+0);
        den = (((((((2.04426310338993978564e-15 * r +1.42151175831644588870e-7) * r + 1.84631831751005468180e-5) * r +7.86869131145613259100e-4) * r +
                     1.48753612908506148525e-2) * r +1.36929880922735805310e-1) * r + 5.99832206555887937690e-1) * r +1.0);
    }
    x = num / den;
    if(q < 0.0)  x = -x;
    return(mu + (x * sigma));

}
