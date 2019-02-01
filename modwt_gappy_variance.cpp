

#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
double SigmaHat (NumericVector w1, NumericVector w2,
		 int m1, int m2,
		 int l1, int l2) {

  int n = w1.size();
  int t0, t1, t, lag;
  double cterm, out;
  
  out = 0.0;
  
  for (lag = -(m1-1); lag <= (m2-1); lag++) {
    
    if ((l2) >= (l1-lag))
      t0 = l2-1;
    else
      t0 = l1-lag-1;
    
    if ((n-lag) >= n)
      t1 = n-1;
    else
      t1 = n-lag-1;
    
    cterm = 0.0;
    for (t = t0; t <= t1; t++)
      cterm += w1[t+(lag)] * w2[t];
    
    out += (cterm * cterm);
  }

  return (out / (double(m1) * double(m2)));
}





// [[Rcpp::export]]
NumericVector CalcYjtHat (NumericVector x, IntegerVector eta,
			  NumericVector hj) {
  
  int l1, l2, s, u, pi_inv, t1, t2;
  int N  = x.size();
  int Lj = hj.size();
  int Mj = N-Lj+1;
  double diffX, alpha;
  NumericVector y(Mj);
  double K = -0.5 * double(Mj);

  for (s=0; s<Mj; s++) {
    y[s] = 0.0;
  }
      
  for (l1=0; l1<Lj; l1++) {
    
    for (l2=0; l2<Lj; l2++) {

      pi_inv = 0;
      for (u=(Lj-1); u<N; u++) {
	pi_inv += eta[u-l1] * eta[u-l2];
      }
      
      alpha = K * hj[l1] * hj[l2] / double(pi_inv);

      for (s=0, t1=Lj-1-l1, t2=Lj-1-l2; s<Mj; s++, t1++, t2++) {
	
	diffX = x[t1] - x[t2];
	y[s] += alpha * double(eta[t1] * eta[t2]) * diffX * diffX;
      }
    }
  }

  return y;
}
