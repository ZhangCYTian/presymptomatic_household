#include <RcppArmadilloExtensions/sample.h>
#include <Rcpp.h>
#include <RcppParallel.h>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppParallel)]]

using namespace Rcpp ;
using namespace RcppParallel;

//########################################################################################################################################
//function for computing the serial interval density
// [[Rcpp::export]]
NumericVector serial_density(double p1, double p2, int shift){
int idx;
NumericVector a(shift+14);
for(idx=0;idx<shift+14;++idx){
  a[idx]=idx+1;
}

NumericVector density=pgamma(a+0.5,p1,p2)-pgamma(a-0.5,p1,p2);

double sum=0;
for(idx=0;idx<shift+14;++idx){
  sum+=density[idx];
}
for(idx=0;idx<shift+14;++idx){
  density[idx]=density[idx]/sum;
}

return density;
}

//########################################################################################################################################
//function for computing the incubation period density
// [[Rcpp::export]]
NumericVector incubation(double p1, double p2){
  NumericVector a=NumericVector::create(1.0,2.0,3.0,4.0,5.0);
  NumericVector density=plnorm(a+0.5,p1,p2)-plnorm(a-0.5,p1,p2);
  int i;
  double sum;
  
  sum=density[0]+density[1]+density[2]+density[3]+density[4];
  
  for(i=0;i<=4;++i){
    density[i]=density[i]/sum;
  }
  
  return density;
}


//########################################################################################################################################
// function to generate normal random variable
// [[Rcpp::export]]
double rnorm(double a, double b) { // a: mean, b: s.d.
	double c=a+b*sum(rnorm(1));
    return c;
}

//########################################################################################################################################
// function to generate binomial random number
// [[Rcpp::export]]
int gen_binom(double p){
double cut=(double)rand()/(RAND_MAX);
int out=0;
if (cut<p){
out=1;
}
return out;
}

//########################################################################################################################################
// function to general multiviariable normal given sigma
// [[Rcpp::export]]
NumericVector rmnorm(arma::mat sigma) {
int ncols=sigma.n_cols;
arma::rowvec c=arma::randn(1,ncols);
arma::rowvec a=c*arma::chol(sigma);   
NumericVector b=NumericVector(a.begin(),a.end());   
return b;
}

//########################################################################################################################################
// function to calculate log-sum
// [[Rcpp::export]]
double log_sum(double a, double b) {
  double exp_a = exp(a);
  double exp_b = exp(b);
  double LL = log(exp_a+exp_b);
  return LL;
}

//########################################################################################################################################
//function to compute the prior likelihood 
// [[Rcpp::export]]
double prior_loglik(NumericVector para, int shift){
// check if the para are within their possible range
NumericVector out(para.length());
int b1;

//out(0)=R::dgamma(pow(1.0/para(0),2.0),1.5,1/0.0001,1); 
//for (b1=1;b1>=0;--b1){
//  out(b1)=R::dunif(para(b1),0.000000000000000001,15,1);
//}

for (b1=0;b1>=0;--b1){
  out(b1)=R::dunif(para(b1),-5,5,1);
}

for (b1=1;b1>=1;--b1){
  out(b1)=R::dunif(para(b1),0.000000000000000001,5,1);
}

for (b1=2;b1>=2;--b1){
  out(b1)=R::dunif(para(b1),shift-5,15,1);
}

for (b1=3;b1>=3;--b1){
  out(b1)=R::dunif(para(b1),0,15,1);
}

for (b1=5;b1>=4;--b1){
  out(b1)=R::dunif(para(b1),0.000000000000000001,10,1);
}

for (b1=30;b1>=6;--b1){
  out(b1)=R::dunif(para(b1),-5.00, 5.00,1);
}

double output=sum(out);
// if the prior is outside the parameter space
if (output< -9999999){
output=-9999999;
}
return output;
}




//########################################################################################################################################
//########################################################################################################################################
// the main body of the parallel for doing simulation
struct SimData:public Worker{
  // source vector
  RMatrix<double> data11;
  RMatrix<double> data1;
  RVector<double> para;
  RVector<double> SI;
  RVector<double> incu;
  RVector<double> inci;
  int shift;
  int sep1;
  int sep2;
  RMatrix<int> record;
  // destination vector
  // initialize with source and destination
  SimData(NumericMatrix data11,
          NumericMatrix data1,
          NumericVector para,
          NumericVector SI,
          NumericVector incu,
          NumericVector inci,
          int shift,
          int sep1,
          int sep2,
          IntegerMatrix record) 
    :data11(data11),data1(data1),para(para),SI(SI),incu(incu),inci(inci),shift(shift),sep1(sep1),sep2(sep2),record(record){}
  void operator()(std::size_t begin, std::size_t end) {
    
    // section to write the parallel version
    // functor (pass input and output matrixes)
    
    for (unsigned int b1=begin;b1<end;++b1){	
      int b2;
      int b3;
      int b4;
      double hazard;
      double sus[int(data11(b1,1))];
      double SI_i[14+shift];
      int count=0;
      int candi[8]={-1,-1,-1,-1,-1,-1,-1,-1};;
      
      // first set every household contacts have 0 infection status
      for (b2=data11(b1,1)-1;b2>=0;--b2){
        
        
        // initial point to set all people is uninfected at the beginning	
        if(data11(b1,b2*sep2+sep1+9)<1){
          data11(b1,b2*sep2+sep1)=0;
          data11(b1,b2*sep2+sep1+1)=-1;
          data11(b1,b2*sep2+sep1+2)=-1;
          data11(b1,b2*sep2+sep1+9)=0;
        }
        
        if(data11(b1,b2*sep2+sep1+9)>0){
          //Rcout << "Member: " << b2 << std::endl;
          candi[count]=b2;
          count+=1;
        }
        
        // data
        //## 1. hhID, 2.size, 3 influenza type ,4-5: follow start and end date 
        //## 5 + 1. infectious status, 2. infection time, 3,4.start and end of follow-up 5.id 
        //## 6.age, 7.sex, 8.vaccination, 9.oseltamivir treatment, 10.community type
        sus[b2]=0;
        sus[b2]+=para[6]*(data11(b1,b2*sep2+sep1+4)<=5);
        sus[b2]+=para[7]*(data11(b1,b2*sep2+sep1+4)>5)*(data11(b1,b2*sep2+sep1+4)<18);
        sus[b2]+=para[8]*(data11(b1,b2*sep2+sep1+4)>50);
        //sus[b2]+=para[9]*(data11(b1,2)==2)*(data11(b1,b2*sep2+sep1+4)>50);
        //sus[b2]+=para[10]*(data11(b1,2)==3)*(data11(b1,b2*sep2+sep1+4)<18);
        //sus[b2]+=para[11]*(data11(b1,2)==3)*(data11(b1,b2*sep2+sep1+4)>50);
        sus[b2]+=para[12]*(data11(b1,b2*sep2+sep1+5)==1);
        //sus[b2]+=para[13]*(data11(b1,2)==2)*(data11(b1,b2*sep2+sep1+5)==1);
        //sus[b2]+=para[14]*(data11(b1,2)==3)*(data11(b1,b2*sep2+sep1+5)==1);
        //sus[b2]+=para[7]*(data11(b1,b2*sep2+sep1+5)==1);
        //sus[b2]+=para[4]*(data11(b1,b2*sep2+sep1+3)==3);
        // age
        //sus[b2]+=para[3]*(data11(b1,1)==2);
      }
      
      //Rcout << "Count: " << count << std::endl;
      int index=(rand()%count);
      b2=candi[index];
      
      double cut=(double)rand()/(RAND_MAX);
      double prob_den=incu[0];
      int incu_time=1;
      while(cut>prob_den&&incu_time<5){
        prob_den=prob_den+incu[incu_time];
        incu_time=incu_time+1;
      }
      data11(b1,b2*sep2+sep1+2)=data11(b1,b2*sep2+sep1+1)-incu_time;
      data11(b1,3)=data11(b1,b2*sep2+sep1+2);
      data11(b1,b2*sep2+sep1+9)=1;
      
      if(count>1){
        for(b3=0;b3<=7;++b3){
          if(b3!=index&&candi[b3]>=0){
            b2=candi[b3];
            data11(b1,b2*sep2+sep1)=0;
            data11(b1,b2*sep2+sep1+1)=-1;
            data11(b1,b2*sep2+sep1+2)=-1;
            data11(b1,b2*sep2+sep1+9)=0;
          }
        }
      }
      
      
      // first need to have a time index
      for (b2=data11(b1,3)+1;b2<=data11(b1,4);++b2){
        // participant index
        // index case do not need update
        for (b3=data11(b1,1)-1;b3>=0;--b3){
          if (data11(b1,b3*sep2+sep1)==0){
            // compute the community hazard
            hazard=para[4]*inci[b2-1];
            // compute the household hazard
            for (b4=data11(b1,1)-1;b4>=0;--b4){
              if ((b4!=b3)&&(data11(b1,b4*sep2+sep1)==1)){
                int idx;
                for(idx=0;idx<=shift+13;++idx){
                  SI_i[idx]=SI[idx];
                }
                if(data11(b1,b4*sep2+sep1+1)-data11(b1,b4*sep2+sep1+2)<shift){
                  int ic;
                  for(ic=shift-int(data11(b1,b4*sep2+sep1+1)-data11(b1,b4*sep2+sep1+2))-1;ic>=0;--ic){
                    SI_i[ic]=0;
                  }
                }
                double sum=0;
                for(idx=0;idx<=shift+13;++idx){
                  sum=sum+SI_i[idx];
                }
                for(idx=0;idx<=shift+13;++idx){
                  SI_i[idx]=SI_i[idx]/sum;
                }
                // need to with the range of serial interval to have contribution
                if ((b2-data11(b1,b4*sep2+sep1+2)>0)&&(b2-data11(b1,b4*sep2+sep1+1)>-shift)&&(b2-data11(b1,b4*sep2+sep1+1)<=14)){
                  double hrisk=para[5];	
                  // here need to add factor affecting infectivity
                  hrisk*=exp(para[15]*(data11(b1,b4*sep2+sep1+7)==1));
                  //hrisk*=exp(para[16]*(data11(b1,2)==2)*(data11(b1,b4*sep2+sep1+7)==1));
                  //hrisk*=exp(para[17]*(data11(b1,2)==3)*(data11(b1,b4*sep2+sep1+7)==1));
                  hrisk*=exp(para[18]*(data11(b1,2)==5));
                  hrisk*=exp(para[19]*(data11(b1,2)==6));
                  hrisk*=exp(para[20]*(data11(b1,1)>=4)*(data11(b1,1)<=5));
                  hrisk*=exp(para[21]*(data11(b1,1)>=6));
                  hrisk*=exp(para[22]*(data11(b1,b4*sep2+sep1+4)<=5));
                  hrisk*=exp(para[23]*(data11(b1,b4*sep2+sep1+4)>5)*(data11(b1,b4*sep2+sep1+4)<18));
                  //hrisk*=exp(para[24]*(data11(b1,2)==3)*(data11(b1,b4*sep2+sep1+4)<18));
                  //hrisk*=exp(para[25]*data11(b1,b4*sep2+sep1+6));
                  //hrisk*=exp(para[26]*(data11(b1,2)==2)*data11(b1,b4*sep2+sep1+6));
                  //hrisk*=exp(para[27]*(data11(b1,2)==3)*data11(b1,b4*sep2+sep1+6));
                  //hrisk*=exp(para[28]*(data11(b1,b4*sep2+sep1+8)==1));
                  //hrisk*=exp(para[29]*(data11(b1,2)==2)*(data11(b1,b4*sep2+sep1+8)==1));
                  //hrisk*=exp(para[30]*(data11(b1,2)==3)*(data11(b1,b4*sep2+sep1+8)==1));
                  hazard+=hrisk*SI_i[b2-int(data11(b1,b4*sep2+sep1+1))+shift-1]/pow(data11(b1,1)-1.0,0);;
                  // need to -1 beacuse the index is from 0 to 9
                }
              }
            }
            // then need to add factor affecting transmission
            hazard*=exp(sus[b3]);
            // generate infeciton
            // infection
            if (gen_binom(1-exp(-hazard))){	
              data11(b1,b3*sep2+sep1)=1;
              data11(b1,b3*sep2+sep1+2)=b2;
              double cut2=(double)rand()/(RAND_MAX);
              double prob_den2=incu[0];
              int incu_time2=1;
              while(cut2>prob_den2&&incu_time2<5){
                prob_den2=prob_den2+incu[incu_time2];
                incu_time2=incu_time2+1;
              }
              data11(b1,b3*sep2+sep1+1)=data11(b1,b3*sep2+sep1+2)+incu_time2;
              if(data11(b1,b3*sep2+sep1+6)<0){
                data11(b1,b3*sep2+sep1+6)=R::rnorm(6.168,0.526)*(data11(b1,b3*sep2+sep1+4)<=17)+R::rnorm(5.442,0.696)*(data11(b1,b3*sep2+sep1+4)>17);
                data11(b1,b3*sep2+sep1+6)=(data11(b1,b3*sep2+sep1+6)-6.678)/0.672;
              }
            }
          }
        }
      }
      
      
      
      
    }
  }
};

//########################################################################################################################################
//function to do simulation
// [[Rcpp::export]]
List sim_data(NumericMatrix data1,
              NumericVector inci,
              NumericVector para,
              int shift,
              int sep1,      // sep1=5
              int sep2){     // sep2=10
  int b1;
  int b2;
  int b3;
  // clone the data first
  NumericMatrix data11(clone(data1));
  
  //NumericVector SI=serial_density(para[2],para[3]);
  NumericVector SI=serial_density(pow(para[2],2)/para[3],para[3]/para[2],shift);
  //NumericVector SI=serial_density(para[2],1,shift);
  
  NumericVector incu=incubation(para[0],para[1]);
  
  IntegerMatrix record(data11.nrow(),100);
  
  
  //1. hhID, 2.size, 3 influenza type ,4-5: follow start and end date 
  //5 + 1. infectious status, 2. infection time, 3,4.start and end of follow-up 5.id 
  //6.age, 7.sex, 8.vaccination, 9.oseltamivir treatment, 10.community type
  
  // compute the serial_density to use
  //NumericVector SI=serial_density(para[0],para[1]);
  
  // call parallel program
  SimData simdata1(data11,data1,para,SI,incu,inci,shift,sep1,sep2,record);
  
  // call parallelFor to do the work
  parallelFor(0,data1.nrow(),simdata1);
  
  
  return List::create(_[""]=data11,
                      _[""]=record);
} 



//########################################################################################################################################
//########################################################################################################################################
// the main body of the parallel for Digraphlikelihood
struct LogLik:public Worker{
// source vector
RMatrix<double> out;
RMatrix<double> out2;
RMatrix<double> data;
RVector<double> para;
RVector<double> SI;
RVector<double> incu;
RVector<double> inci;
int shift;
int sep1;
int sep2;
RMatrix<double> record;
// destination vector
// initialize with source and destination
LogLik(NumericMatrix out,
NumericMatrix out2,	
NumericMatrix data,
NumericVector para,
NumericVector SI,
NumericVector incu,
NumericVector inci,
int shift,
int sep1,
int sep2,
NumericMatrix record) 
  :out(out),out2(out2),data(data),para(para),SI(SI),incu(incu),inci(inci),shift(shift),sep1(sep1),sep2(sep2),record(record){}
void operator()(std::size_t begin, std::size_t end) {

// section to write the parallel version
// functor (pass input and output matrixes)
for (unsigned int b1=begin;b1<end;++b1){
int b2;
int b3;
int b4;
int b5;
double SI_i[shift+14];

// b3 is the participant index
for (b3=data(b1,1)-1;b3>=0;--b3){
if(data(b1,b3*sep2+sep1+9)==0){
//if (!((data(b1,b3*sep2+sep1)==1)&(data(b1,b3*sep2+sep1+1)==data(b1,sep1+1)))){
// need to add factor addecting susceptibility
double sus=0;
  sus+=para[6]*(data(b1,b3*sep2+sep1+4)<=5);
  sus+=para[7]*(data(b1,b3*sep2+sep1+4)>5)*(data(b1,b3*sep2+sep1+4)<18);
  sus+=para[8]*(data(b1,b3*sep2+sep1+4)>50);
  //sus+=para[9]*(data(b1,2)==2)*(data(b1,b3*sep2+sep1+4)>50);
  //sus+=para[10]*(data(b1,2)==3)*(data(b1,b3*sep2+sep1+4)<18);
  //sus+=para[11]*(data(b1,2)==3)*(data(b1,b3*sep2+sep1+4)>50);
  sus+=para[12]*(data(b1,b3*sep2+sep1+5)==1);
  //sus+=para[13]*(data(b1,2)==2)*(data(b1,b3*sep2+sep1+5)==1);
  //sus+=para[14]*(data(b1,2)==3)*(data(b1,b3*sep2+sep1+5)==1);

// the final date for contribution from non-infection
int finaltime=data(b1,4)+1;
if (data(b1,b3*sep2+sep1)==1){
  finaltime=data(b1,b3*sep2+sep1+2);	
}

if(finaltime>data(b1,3)){
  double h[finaltime-int(data(b1,3))];
  // fill the community risk
  for (b2=finaltime-data(b1,3)-1;b2>=0;--b2){
    h[b2]=para[4]*inci[data(b1,3)+b2-1];	
  }
  
  
  for (b4=data(b1,1)-1;b4>=0;--b4){
    if ((b4!=b3)&&(data(b1,b4*sep2+sep1)==1)){
      
      int idx;
      for(idx=0;idx<=shift+13;++idx){
        SI_i[idx]=SI[idx];
      }
      if(data(b1,b4*sep2+sep1+1)-data(b1,b4*sep2+sep1+2)<shift){
        int ic;
        for(ic=shift-int(data(b1,b4*sep2+sep1+1)-data(b1,b4*sep2+sep1+2))-1;ic>=0;--ic){
          SI_i[ic]=0;
        }
      }
      double sum=0;
      for(idx=0;idx<=shift+13;++idx){
        sum=sum+SI_i[idx];
      }
      for(idx=0;idx<=shift+13;++idx){
        SI_i[idx]=SI_i[idx]/sum;
      }
      
      for (b5=14;b5>=-shift+1;--b5){
        //if (data(b1,b4*sep2+sep1+1)+b5<=data(b1,4)){
        // if infection date of individual b3 is smaller than the final time	
        if ((data(b1,b4*sep2+sep1+1)+b5<=finaltime)&&(data(b1,b4*sep2+sep1+1)+b5>data(b1,b4*sep2+sep1+2))){ 
          double hrisk=para[5];
          // here need to add factor affecting infectivity
          hrisk*=exp(para[15]*(data(b1,b4*sep2+sep1+7)==1));
          //hrisk*=exp(para[16]*(data(b1,2)==2)*(data(b1,b4*sep2+sep1+7)==1));
          //hrisk*=exp(para[17]*(data(b1,2)==3)*(data(b1,b4*sep2+sep1+7)==1));
          hrisk*=exp(para[18]*(data(b1,2)==5));
          hrisk*=exp(para[19]*(data(b1,2)==6));
          hrisk*=exp(para[20]*(data(b1,1)>=4)*(data(b1,1)<=5));
          hrisk*=exp(para[21]*(data(b1,1)>=6));
          hrisk*=exp(para[22]*(data(b1,b4*sep2+sep1+4)<=5));
          hrisk*=exp(para[23]*(data(b1,b4*sep2+sep1+4)>5)*(data(b1,b4*sep2+sep1+4)<18));
          //hrisk*=exp(para[24]*(data(b1,2)==3)*(data(b1,b4*sep2+sep1+4)<18));
          //hrisk*=exp(para[25]*data(b1,b4*sep2+sep1+6));
          //hrisk*=exp(para[26]*(data(b1,2)==2)*data(b1,b4*sep2+sep1+6));
          //hrisk*=exp(para[27]*(data(b1,2)==3)*data(b1,b4*sep2+sep1+6));
          //hrisk*=exp(para[28]*(data(b1,b4*sep2+sep1+8)==1));
          //hrisk*=exp(para[29]*(data(b1,2)==2)*(data(b1,b4*sep2+sep1+8)==1));
          //hrisk*=exp(para[30]*(data(b1,2)==3)*(data(b1,b4*sep2+sep1+8)==1));
          h[int(data(b1,b4*sep2+sep1+1)-data(b1,3))+b5-1]+=hrisk*SI_i[b5+shift-1]/pow(data(b1,1)-1.0,0);;
        }
      }
    }
  }
  //}
  
  
  for (b2=finaltime-data(b1,3)-2;b2>=0;--b2){
    out(b1,b3)-=h[b2]*exp(sus);	
  }
  if (data(b1,b3*sep2+sep1)==1){
    out(b1,b3)+=log(1-exp(-h[finaltime-int(data(b1,3))-1]*exp(sus)));
  }
  
  if (b3==0){
    for (b2=finaltime-data(b1,3)-1;b2>=0;--b2){
      record(b1,b2)=h[b2];
    }
  }
}
if(finaltime<=data(b1,3)){
  out(b1,b3)=-9999999;
}


}

if(data(b1,b3*sep2+sep1+9)==2){
  double sus=0;
  sus+=para[6]*(data(b1,b3*sep2+sep1+4)<=5);
  sus+=para[7]*(data(b1,b3*sep2+sep1+4)>5)*(data(b1,b3*sep2+sep1+4)<18);
  sus+=para[8]*(data(b1,b3*sep2+sep1+4)>50);
  //sus+=para[9]*(data(b1,2)==2)*(data(b1,b3*sep2+sep1+4)>50);
  //sus+=para[10]*(data(b1,2)==3)*(data(b1,b3*sep2+sep1+4)<18);
  //sus+=para[11]*(data(b1,2)==3)*(data(b1,b3*sep2+sep1+4)>50);
  sus+=para[12]*(data(b1,b3*sep2+sep1+5)==1);
  //sus+=para[13]*(data(b1,2)==2)*(data(b1,b3*sep2+sep1+5)==1);
  //sus+=para[14]*(data(b1,2)==3)*(data(b1,b3*sep2+sep1+5)==1);
  
  out(b1,b3)+=log(1-exp(-para[4]*inci[data(b1,b3*sep2+sep1+2)-1]*exp(sus)));
}

}



for (b3=data(b1,1)-1;b3>=0;--b3){
// input the rm for infectivity
if (data(b1,b3*sep2+sep1)==1){
//out2(b1,b3)=R::dnorm(data(b1,b3*sep2+sep1+2),0.0,para[0],1);
  out2(b1,b3)=0;
}
}


}
}
};

//########################################################################################################################################
//function to likelihood
// [[Rcpp::export]]
List loglik(NumericMatrix data,
NumericVector inci,
NumericVector para,
int shift,
int sep1,
int sep2){
// check if the para are within their possible range
NumericMatrix out(data.nrow(),30);
NumericMatrix out2(data.nrow(),30);
NumericMatrix record(data.nrow(),100);
NumericMatrix total_ll(data.nrow(),2);
NumericVector SI=serial_density(pow(para[2],2)/para[3],para[3]/para[2],shift);
//NumericVector SI=serial_density(para[2],1,shift);
NumericVector incu=incubation(para[0],para[1]);
// call parallel program
LogLik loglik(out,out2,data,para,SI,incu,inci,shift,sep1,sep2,record);
// call parallelFor to do the work
parallelFor(0,data.nrow(),loglik);

total_ll(_,0)=data(_,0);
total_ll(_,1)=data(_,5);

int bi;
for(bi=data.nrow()-1;bi>=0;--bi){
  total_ll(bi,1)+=sum(out(bi,_));
}

NumericMatrix out_ll(data.nrow(),2);
out_ll(0,0)=total_ll(0,0);
out_ll(0,1)=total_ll(0,1);
int count=0;
for(bi=1;bi<=data.nrow()-1;++bi){
  if(total_ll(bi,0)==out_ll(count,0)){
    out_ll(count,1)=log_sum(out_ll(count,1),total_ll(bi,1));
  }
  if(total_ll(bi,0)!=out_ll(count,0)){
    count+=1;
    out_ll(count,0)=total_ll(bi,0);
    out_ll(count,1)=total_ll(bi,1);
  }
}

double sum_ll = sum(out_ll(_,1));

return List::create(_[""]=sum_ll,
  _[""]=out,
  _[""]=total_ll,
	_[""]=out_ll);

}




//##############################################################################################################################################
//##############################################################################################################################################
// function for mcmc
// [[Rcpp::export]]
List mcmc(NumericMatrix data1,
NumericVector inci,
int shift,
int mcmc_n,             // length of mcmc stain
NumericVector int_para, // initial parameter
NumericVector move,     // which one should move in the model
NumericVector sigma){            

// create the vector for use
int b0;
int b1;
int b2;
int b3;
int b4;
int sep1=6;
int sep2=10;
int moveindex;

//####################################################################################################################################
// backup data 
NumericMatrix data11(clone(data1));



// matrix to record LL
// need to set number of parameter here
NumericMatrix p_para(mcmc_n,int_para.length());
NumericMatrix p_para_r(mcmc_n,sum(move));
p_para(0,_)=int_para;
moveindex=sum(move)-1;
for (b1=int_para.length()-1;b1>=0;--b1){
if (move(b1)){
p_para_r(0,moveindex)=p_para(0,b1);
--moveindex;
}	
}

// first row is the overall matrix, other three row is the indiviudal likelihood
NumericVector acceptrate(int_para.length());
NumericMatrix LL1(mcmc_n,3);
NumericMatrix LL2(mcmc_n,3);

//####################################################################################################################################
// initial step

//####################################################################################################################################
// compute likelihood
List loglikall=loglik(data11,inci,p_para(0,_),shift,sep1,sep2);
List loglikallpro;
//NumericMatrix loglik1=loglikall(0);
//NumericMatrix loglik2=loglikall(1);
//NumericMatrix loglik1pro;
//NumericMatrix loglik2pro;

LL1(0,1)=loglikall(0);
LL1(0,2)=0;
LL1(0,0)=LL1(0,1)+LL1(0,2);

NumericVector temploglik(3);
NumericVector newloglik(3);
temploglik(0)=LL1(0,0)+prior_loglik(p_para(0,_),shift);
temploglik(1)=LL1(0,1);
temploglik(2)=LL1(0,2);


double loglikeratio;
double accept_pro;
NumericVector pro_para(int_para.length());



//####################################################################################################################################
// main mcmc step
NumericMatrix updateacceptrate(mcmc_n,20);

//####################################################################################################################################
for (b0=1;b0<mcmc_n;++b0){

// after 500 step, then set the sigma to be the empirical sigma
if (b0>500){
for (b1=int_para.length()-1;b1>=0;--b1){
if (move(b1)){	
NumericVector temp1(b0-1);
for (b2=b0-2;b2>=0;--b2){
temp1(b2)=p_para(b2,b1);	
}
sigma(b1)=sd(temp1);
// tuning
if (acceptrate(b1)<0.1){
sigma(b1)*=0.5;
}	
if ((acceptrate(b1)<0.15)&(acceptrate(b1)>0.1)){
sigma(b1)*=0.8;
}
if ((acceptrate(b1)<0.2)&(acceptrate(b1)>0.15)){
sigma(b1)*=0.95;
}
if ((acceptrate(b1)<0.4)&(acceptrate(b1)>0.3)){
sigma(b1)*=1.05;
}
if ((acceptrate(b1)<0.9)&(acceptrate(b1)>0.4)){
sigma(b1)*=1.2;
}
if (acceptrate(b1)>0.9){
sigma(b1)*=2;
}
}
}
}

// metorpolis-hasing update on parameter
for (b1=0;b1<int_para.length();++b1){
if (move(b1)){
pro_para=p_para(b0-1,_);
for (b2=b1-1;b2>=0;--b2){
pro_para(b2)=p_para(b0,b2);	
}
pro_para(b1)+=rnorm(0.0,sigma(b1));
newloglik(0)=prior_loglik(pro_para,shift);
if (newloglik(0)> -9999999){
loglikallpro=loglik(data11,inci,pro_para,shift,sep1,sep2);
//NumericMatrix tempoutput=loglikallpro(0);
//NumericMatrix tempoutput2=loglikallpro(1);
//loglik1pro=clone(tempoutput);
//loglik2pro=clone(tempoutput2);
newloglik(1)=loglikallpro(0);
newloglik(2)=0;
newloglik(0)+=newloglik(1)+newloglik(2);
loglikeratio=newloglik(0)-temploglik(0);
accept_pro=pow(exp(1),loglikeratio);
}
else{
accept_pro=0;	
}
if(gen_binom(accept_pro)){
p_para(b0,b1)=pro_para(b1);
//loglik1=clone(loglik1pro);
//loglik2=clone(loglik2pro);		
temploglik(2)=newloglik(2);
temploglik(1)=newloglik(1);
temploglik(0)=newloglik(0);
acceptrate(b1)*=(b0-1);
acceptrate(b1)+=1;
acceptrate(b1)/=b0;
}
else{
p_para(b0,b1)=p_para(b0-1,b1);
acceptrate(b1)*=(b0-1);
acceptrate(b1)/=b0;
}
}
else {
p_para(b0,b1)=p_para(b0-1,b1);
}
}

LL1(b0,0)=temploglik(0)-prior_loglik(p_para(b0,_),shift);
LL1(b0,1)=temploglik(1);
LL1(b0,2)=temploglik(2);

// move the matirx to another matrix to store the parameter and compute the correlation matrix
moveindex=sum(move)-1;
for (b1=int_para.length()-1;b1>=0;--b1){
if (move(b1)){
p_para_r(b0,moveindex)=p_para(b0,b1);
--moveindex;
}	
}


if (b0%100==0){
Rcout << b0 << std::endl;
}

}


return List::create(_[""]=p_para,
_[""]=LL1);
} 




//##############################################################################################################################################
//##############################################################################################################################################
// function for compute DIC
// [[Rcpp::export]]
NumericMatrix DIC(NumericMatrix data1,
                  NumericMatrix para,
                  NumericMatrix para2, // it is the rm
                  NumericVector inci,
                  int simtime,
                  int shift,
                  int sep1,
                  int sep2){
  
  NumericMatrix LL(para.nrow(), data1.nrow());
  
  double match;
  int b0;
  int b1;
  int b2;
  int b3;
  
  for (b3=para.nrow()-1;b3>=0;--b3){
    Rcout << "b3: " << b3 << std::endl;	
    //NumericVector mat1(data1.nrow());
    
    for (b0=simtime-1;b0>=0;--b0){
      //List datasim=cond_sim_data(data1,SI,para2(b3,_),para(b3,_),sep1,sep2);	
      List datasim=sim_data(data1,inci,para(b3,_),shift,sep1,sep2);	
      NumericMatrix data11=datasim(0);
      for (b1=data1.nrow()-1;b1>=0;--b1){
        match=1;
        for (b2=data1(b1,1)-1;b2>=0;--b2){
          if (data1(b1,b2*sep2+sep1)!=-1){	
            if (data11(b1,b2*sep2+sep1)!=data1(b1,b2*sep2+sep1)){
              match*=0.001;
            }
            else{
              match*=0.999;	
            }	
          }
        }
        LL(b3,b1)+=match;
        //mat1(b1)+=match;
      }
    }
    //LL(b3)=sum(log(mat1/simtime));	
    
  }
  
  return LL;
} 


