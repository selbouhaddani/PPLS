#include <iostream>
#include <RcppEigen.h>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <Eigen/Cholesky>
//#include <boost/timer/timer.hpp>
//#include <cmath>
// [[Rcpp::depends(RcppEigen)]]
//// [[Rcpp::depends(BH)]]

using namespace Eigen;
using namespace std;
using namespace Rcpp;


// [[Rcpp::export]]
Eigen::MatrixXd invUpdate(double sig, Eigen::VectorXd W)
{
  int p = W.size();   // number of X variables, can be 1e1 to 1e4
  VectorXd SinvW = W / sig;
  MatrixXd numer = SinvW * SinvW.transpose();
  double denom = 1 + (W.transpose() * SinvW)(0,0);
  
  MatrixXd Sinv = - numer / denom;
  for(int i=0;i<p;i++){Sinv(i,i) += sig;} // add double to diagonal
  
  return Sinv;
}

// ---------------- This is one of the two CORE functions
// ---------------- It evaluates the likelihood at each step
// [[Rcpp::export]]
double loglC(Eigen::MatrixXd W,Eigen::MatrixXd C,Eigen::MatrixXd P_Yosc,Eigen::MatrixXd P_Xosc,Eigen::MatrixXd B_T, Eigen::MatrixXd Dat,
             double sig2X,double sig2Y,Eigen::MatrixXd sig2H,Eigen::MatrixXd sig2T,Eigen::MatrixXd sig2To,Eigen::MatrixXd sig2Uo)
{
  int N = Dat.rows(); // number of samples, typically of order 1e2 - 1e3 
  int p = W.rows();   // number of X variables, can be 1e1 to 1e4
  int q = C.rows();   // number of Y variables, say can be 1e1 to 1e4
  
  // ----------- initialize, fill and invert theoretical cov of cbind(X,Y)
  MatrixXd SX; // LARGE will be p times p
  SX = W*sig2T*W.transpose() + P_Yosc*sig2To*P_Yosc.transpose(); // rank low
  for(int i=0;i<p;i++){SX(i,i) += sig2X;} // add double to diagonal
  
  MatrixXd SY; // LARGE will be q times q
  SY = C*sig2T*B_T*B_T*C.transpose() + C*sig2H*C.transpose() + P_Xosc*sig2Uo*P_Xosc.transpose(); //rank low
  for(int i=0;i<q;i++){SY(i,i) += sig2Y;} // add double to diagonal
  
  MatrixXd SXY; // LARGE will be p times q
  SXY = SXY = W*sig2T* B_T*C.transpose(); // only rank low
  
  // initialize NewSX with dimensions (p+q) times (p+q)
  MatrixXd NewSX(SX.rows()+SY.rows(),SX.cols() + SXY.cols()); 
  NewSX << SX,SXY,
           SXY.transpose(),SY; // equivalent to matlab [SX SXY ; SXY' SY]
  LLT<MatrixXd> lltSX(NewSX);   // calc cholesky decomposition (Eigen/Cholesky library)
  
  //Inversion step can be done maybe with "rank one update", in R I would use chol2inv()
  MatrixXd invS = lltSX.solve(MatrixXd::Identity(p+q,p+q)); // Invert using cholesky
  // -----------
  
  double som2 = (Dat.transpose()*Dat*invS).trace(); //trace of product of LARGE matrices
  // ----- I can calc this directly I think
  MatrixXd L = lltSX.matrixL(); 
  MatrixXd Logdiag;
  Logdiag = L.diagonal(); // The diagonal of L are the sqrt of the eigenvalues
  for(int j=0;j<(p+q);j++){
    Logdiag(j) = 2*log(Logdiag(j)); // 2*log since we need to square the sqrt(eigenval's)
  }
  // -----
  double Loglik = - 0.5*N*Logdiag.sum() - 0.5 * som2;
  return Loglik;
  
}

// [[Rcpp::export]]
List EMstepC(Eigen::MatrixXd W,Eigen::MatrixXd C, Eigen::MatrixXd P_Yosc, Eigen::MatrixXd P_Xosc,Eigen::MatrixXd B_T, Eigen::MatrixXd Dat,
            double sig2X,double sig2Y,Eigen::MatrixXd sig2H,Eigen::MatrixXd sig2T,Eigen::MatrixXd sig2To,Eigen::MatrixXd sig2Uo,bool onlyW = false)
{
  int N = Dat.rows();
  int p = W.rows();
  int q = C.rows();
  int a = W.cols();
  int nx = P_Yosc.cols();
  int ny = P_Xosc.cols();
  
  // copied from loglC, calc theoretical cov matrix, may be faster to reuse from loglC
  // ----------- initialize, fill and invert theoretical cov of cbind(X,Y)
  MatrixXd SX; // LARGE will be p times p
  SX = W*sig2T*W.transpose() + P_Yosc*sig2To*P_Yosc.transpose(); // rank low
  for(int i=0;i<p;i++){SX(i,i) += sig2X;} // add double to diagonal
  
  MatrixXd SY; // LARGE will be q times q
  SY = C*sig2T*B_T*B_T*C.transpose() + C*sig2H*C.transpose() + P_Xosc*sig2Uo*P_Xosc.transpose(); //rank low
  for(int i=0;i<q;i++){SY(i,i) += sig2Y;} // add double to diagonal
  
  MatrixXd SXY; // LARGE will be p times q
  SXY = SXY = W*sig2T* B_T*C.transpose(); // only rank low
  
  // initialize NewSX with dimensions (p+q) times (p+q)
  MatrixXd NewSX(SX.rows()+SY.rows(),SX.cols() + SXY.cols());
  NewSX << SX,SXY,
           SXY.transpose(),SY; // equivalent to matlab [SX SXY ; SXY' SY]
  LLT<MatrixXd> lltSX(NewSX);   // calc cholesky decomposition (Eigen/Cholesky library)
  
  //Inversion step can be done maybe with "rank one update", in R I would use chol2inv()
  MatrixXd invS = lltSX.solve(MatrixXd::Identity(p+q,p+q)); // Invert using cholesky
  // -----------
  
   // // // // //  Expectations------------------------------------------------------------------------------
   MatrixXd Cxxyy = Dat.transpose()*Dat / N;
   MatrixXd Cxx = Cxxyy.topLeftCorner(p, p);
   MatrixXd Cxy = Cxxyy.topRightCorner(p, q);
   MatrixXd Cyy = Cxxyy.bottomRightCorner(q, q);
   
   // // // For T ----------------------------
   // // Expected sample cov matrix for T
   MatrixXd covT(p+q,a);
   // equivalent to rbind(W,B_T*C) or c(W,B_T*C)
   covT << W*sig2T,
   C*sig2T*B_T;
   MatrixXd delT = invS * covT;
   MatrixXd DelT = sig2T - (covT.transpose() * invS * covT);
   MatrixXd Cxt = (MatrixXd(p,p+q) << Cxx,Cxy).finished() * delT;
   MatrixXd Ctt = (delT.transpose() * Cxxyy * delT) + DelT;
   
   // // Expected sample cov matrix for To
   MatrixXd covTo(p+q,nx);
   covTo << P_Yosc*sig2To,
            MatrixXd::Zero(q,nx);
   MatrixXd delTo = invS * covTo;
   MatrixXd DelTo = sig2To - (covTo.transpose() * invS * covTo);
   MatrixXd Cxto = (MatrixXd(p,p+q) << Cxx,Cxy).finished() * delTo;
   MatrixXd Ctoto = (delTo.transpose() * Cxxyy * delTo) + DelTo;
   
   MatrixXd Ctto =  - (covT.transpose() * invS * covTo) + (covT.transpose() * invS * Cxxyy * invS * covTo);
   
   // // For U ----------------------------
   // Expected sample cov matrix for U
   MatrixXd covU(p+q,a);
   // equivalent to rbind(W,B_T*C) or c(W,B_T*C)
   covU << W*sig2T*B_T,
   C*sig2T*B_T*B_T + C*sig2H;
   MatrixXd delU = invS * covU;
   MatrixXd DelU = sig2T*B_T*B_T+sig2H - (covU.transpose() * invS * covU);
   MatrixXd Cyu = (MatrixXd(q,p+q) << Cxy.transpose(),Cyy).finished() * delU;
   MatrixXd Cuu = (delU.transpose() * Cxxyy * delU) + DelU;
   
   // Expected sample cov matrix for Uo
   MatrixXd covUo(p+q,ny);
   covUo << MatrixXd::Zero(p,ny),
   sig2Uo * P_Xosc;
   MatrixXd delUo = invS * covUo;
   MatrixXd DelUo = sig2Uo - (covUo.transpose() * invS * covUo);
   MatrixXd Cyuo = (MatrixXd(q,p+q) << Cxy.transpose(),Cyy).finished() * delUo;
   MatrixXd Cuouo = (delUo.transpose() * Cxxyy * delUo) + DelUo;
   
   MatrixXd Cuuo =  - (covU.transpose() * invS * covUo) + (covU.transpose() * invS * Cxxyy * invS * covUo);
   
   // // For E, F and H -------------------
   MatrixXd SX0(p+q,p);MatrixXd SY0(p+q,q);MatrixXd SH0(p+q,a);
   SX0 << sig2X*MatrixXd::Identity(p,p),
          MatrixXd::Zero(q,p);
   SY0 << MatrixXd::Zero(p,q),
          sig2Y*MatrixXd::Identity(q,q);
   SH0 << MatrixXd::Zero(p,a),
          C*sig2H;
   
   // Can be done faster in some way...        
   MatrixXd Cee = sig2X*MatrixXd::Identity(p,p) - 
   SX0.transpose() * invS * SX0 + SX0.transpose() * invS * Cxxyy * invS * SX0;
   MatrixXd Cff = sig2Y*MatrixXd::Identity(p,p) - 
   SY0.transpose() * invS * SY0 + SY0.transpose() * invS * Cxxyy * invS * SY0;
   MatrixXd Chh = sig2H*MatrixXd::Identity(a,a) - 
   SH0.transpose() * invS * SH0 + SH0.transpose() * invS * Cxxyy * invS * SH0;
   
   MatrixXd Cut = sig2T*B_T - (covU.transpose()*invS*covT) + (delU.transpose()*Cxxyy*delT);
   
   // ------------- Maximization step is done in R ( with multiroot() )
   
   List ret;
   ret["Cxt"] = Cxt;
   ret["Cxto"] = onlyW*Cxto;
   ret["Ctt"] = Ctt;
   ret["Ctto"] = onlyW*Ctto;
   ret["Ctoto"] = onlyW*Ctoto;
   
   ret["Cyu"] = Cyu;
   ret["Cyuo"] = onlyW*Cyuo;
   ret["Cuu"] = Cuu;
   ret["Cuuo"] = onlyW*Cuuo;
   ret["Cuouo"] = onlyW*Cuouo;
   
   ret["Cut"] = Cut;
   
   ret["diag_Cee"] = Cee.diagonal();
   ret["diag_Cff"] = Cff.diagonal();
   ret["Chh"] = Chh;
   
   return ret;
}


// [[Rcpp::export]]
Eigen::MatrixXd simulC(int N, Eigen::MatrixXd W,Eigen::MatrixXd C, Eigen::MatrixXd P_Yosc, Eigen::MatrixXd P_Xosc,Eigen::MatrixXd B_T,
                       double sigX,double sigY,Eigen::MatrixXd sigH,Eigen::MatrixXd sigT,Eigen::MatrixXd sigTo,Eigen::MatrixXd sigUo)
{
  RNGScope scope;    // ensure RNG gets set/reset
  int p = W.rows();
  int q = C.rows();
  int a = W.cols();
  int nx = P_Yosc.cols();
  int ny = P_Xosc.cols();
  
  // Simulate LV's with sd's
  NumericMatrix t1(N,a);
  NumericMatrix to1(N,nx);
  NumericMatrix uo1(N,ny);
  NumericMatrix h1(N,a);
  
  NumericMatrix e1(N,p);
  NumericMatrix f1(N,q);
  
  for(int i=0;i<p;i++){ e1(_,i) = rnorm(N); }
  for(int i=0;i<q;i++){ f1(_,i) = rnorm(N); }
  for(int i=0;i<a;i++){ t1(_,i) = rnorm(N); }
  for(int i=0;i<nx;i++){ to1(_,i) = rnorm(N); }
  for(int i=0;i<ny;i++){ uo1(_,i) = rnorm(N); }
  for(int i=0;i<a;i++){ h1(_,i) = rnorm(N); }
  
  
  // convert to MatrixXd (Eigen)
  Map<MatrixXd> Tt = as<Map<MatrixXd> >(t1);
  Map<MatrixXd> TYo = as<Map<MatrixXd> >(to1);
  Map<MatrixXd> UXo = as<Map<MatrixXd> >(uo1);
  
  Map<MatrixXd> H = as<Map<MatrixXd> >(h1);
  Map<MatrixXd> E = as<Map<MatrixXd> >(e1);
  Map<MatrixXd> Ff = as<Map<MatrixXd> >(f1);
  
  Tt *= sigT;
  TYo *= sigTo;
  UXo *= sigUo;
  H *= sigH;
  E *= sigX;
  Ff *= sigY;
  
  VectorXd U = Tt * B_T + H;
  MatrixXd X = Tt*W.transpose() + TYo*P_Yosc.transpose() + E;
  MatrixXd Y = U*C.transpose() + UXo*P_Xosc.transpose() + Ff;
  MatrixXd Dat(N,p+q);
  Dat << X, Y;
  return Dat;
}
