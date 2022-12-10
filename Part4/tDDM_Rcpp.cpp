#include <Rcpp.h>
#include <RcppParallel.h>
#include <math.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <time.h>

using namespace Rcpp;
using namespace RcppParallel;

typedef boost::mt19937                     ENG;    // Mersenne Twister
typedef boost::normal_distribution<double> DIST;   // Normal Distribution
typedef boost::variate_generator<ENG,DIST> GEN;    // Variate generator



// [[Rcpp::depends(RcppParallel)]]
struct ddm2w : public Worker {
  
  // input vectors/matrices to read from
  const double drate_m; // the drift rate (weight) for monetary reward
  const double drate_c; // the drift rate (weight) for consistency

  const double lat_m; // latency for monetary reward: the time money comes into the decision process
  const double lat_c; // latency for consistency

  const double moneydiff;
  const double consdiff;
  
  
  GEN gen;
  
  // output vector to write to
  RVector<double> vecOut;
  
  // initialize from Rcpp input and output matrices/vectors (the RMatrix/RVector class
  // can be automatically converted to from the Rcpp matrix/vector type)
  ddm2w(const double drate_m, const double drate_c,  const double lat_m, const double lat_c, const double moneydiff, const double consdiff, NumericVector vecOut , GEN gen)
    : drate_m(drate_m), drate_c(drate_c), lat_m(lat_m),lat_c(lat_c),  moneydiff(moneydiff), consdiff(consdiff), gen(gen), vecOut(vecOut) {}
  
  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {
    
    double T = 4, dt , lt=1024;
    //double T = 50, dt = 0.001, lt;
    dt = (T/lt);
    
    
    std::vector<double> vec_tCons(lt,0);
    std::vector<double> vec_tMoney(lt,0);
    
    
    
    for (int t=lat_c/dt; t<lt; t++) {
        vec_tCons[t] = 1;
    }
    for (int t=lat_m/dt; t<lt; t++) {
    vec_tMoney[t] = 1;
    }
   
    for (std::size_t i = begin; i < end; i++) {
      vecOut[i] = T;
      double X = 0;
      int flag = 0;
      double cont = 0;
      double noise = 0;
      double thres=2.5;
      
      while (flag==0 && cont<lt) {
        
        noise=gen()*sqrt(dt);
        X = X + (drate_m * moneydiff * vec_tMoney[cont] + drate_c * vec_tCons[cont])*dt + noise;
        
        if (X > thres) {
          flag=1;
          //vecOut[i] = nDT + cont*dt;
          vecOut[i] = cont*dt;
        }
        else if (X < -thres) {
          flag=1;
          //vecOut[i] = -nDT -cont*dt;
          vecOut[i] = -cont*dt;
        }
        cont++;
        
      }
    }
  }
};


// [[Rcpp::export]]
NumericVector tddm_parallel(double drate_m, double drate_c, double lat_m, double lat_c, double moneydiff, double consdiff, unsigned int N) {
  
  const double sd_n = 1.4;
  struct timespec time;
  clock_gettime(CLOCK_REALTIME, &time);
  ENG  eng;
  eng.seed(time.tv_nsec);
  DIST dist(0,sd_n);
  GEN  gen(eng,dist);
  
  //output vector
  NumericVector vecOut(N);
  
  // create the worker
  ddm2w ddm2w(drate_m, drate_c, lat_m, lat_c, moneydiff, consdiff, vecOut, gen);
  
  // call the worker
  parallelFor(0, N, ddm2w);
  
  return vecOut;
}
