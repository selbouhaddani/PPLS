#include <iostream>
#include <RcppEigen.h>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <Eigen/Cholesky>
// [[Rcpp::depends(RcppEigen)]]

using namespace Eigen;
using namespace std;
using namespace Rcpp;

// ---------------- This is one of the two CORE functions
// ---------------- It evaluates the likelihood at each step
// [[Rcpp::export]]
double loglC(Eigen::VectorXd W,Eigen::VectorXd C,Eigen::VectorXd P_Yosc,Eigen::VectorXd P_Xosc,double B_T, Eigen::MatrixXd Dat,
                  double sig2X,double sig2Y,double sig2H,double sig2T,double sig2To,double sig2Uo) 
{
  int N = Dat.rows(); // number of samples, typically of order 1e2 - 1e3 
  int p = W.size();   // number of X variables, can be 1e1 to 1e4
  int q = C.size();   // number of Y variables, say can be 1e1 to 1e4
  
  // ----------- initialize, fill and invert theoretical cov of cbind(X,Y)
  MatrixXd SX; // LARGE will be p times p
  SX = sig2T*W*W.transpose() + sig2To*P_Yosc*P_Yosc.transpose(); // rank two
  for(int i=0;i<p;i++){SX(i,i) += sig2X;} // add double to diagonal
  
  MatrixXd SY; // LARGE will be q times q
  SY = (sig2T*B_T*B_T + sig2H)*C*C.transpose() + sig2Uo*P_Xosc*P_Xosc.transpose(); //rank two
  for(int i=0;i<q;i++){SY(i,i) += sig2Y;} // add double to diagonal
  
  MatrixXd SXY; // LARGE will be p times q
  SXY = SXY = sig2T* B_T*W*C.transpose(); // only rank one
  
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
List EMstepC(Eigen::MatrixXd W,Eigen::MatrixXd C, Eigen::MatrixXd P_Yosc, Eigen::MatrixXd P_Xosc,double B_T, Eigen::MatrixXd Dat,
              double sig2X,double sig2Y,double sig2H,double sig2T, double sig2To, double sig2Uo)
{
  int N = Dat.rows();
  int p = W.size();
  int q = C.size();
  
  // copied from loglC, calc theoretical cov matrix, may be faster to reuse from loglC
  MatrixXd SX;
  SX = sig2T * W * W.transpose() + sig2To * P_Yosc * P_Yosc.transpose();
  for(int i=0;i<p;i++){SX(i,i) += sig2X;}
  MatrixXd SY;
  SY = (sig2T*B_T*B_T+sig2H) * C * C.transpose() + sig2Uo * P_Xosc * P_Xosc.transpose();
  for(int i=0;i<q;i++){SY(i,i) += sig2Y;}
  MatrixXd SXY;
  SXY = SXY = sig2T * B_T * W * C.transpose();
  MatrixXd NewSX(SX.rows()+SY.rows(),SX.cols() + SY.cols());
  NewSX << SX,SXY,
            SXY.transpose(),SY;
  LLT<MatrixXd> lltSX(NewSX);
  MatrixXd invS = lltSX.solve(MatrixXd::Identity(p+q,p+q));
  
  // // // // //  Expectations------------------------------------------------------------------------------
  MatrixXd Cxxyy = Dat.transpose()*Dat / N;
  MatrixXd Cxx = Cxxyy.topLeftCorner(p, p);
  MatrixXd Cxy = Cxxyy.topRightCorner(p, q);
  MatrixXd Cyy = Cxxyy.bottomRightCorner(q, q);
  // // // For T ----------------------------
  // // Expected sample cov matrix for T
  MatrixXd covT(p+q,1);
  // equivalent to rbind(W,B_T*C) or c(W,B_T*C)
  covT << sig2T*W,
        sig2T*B_T*C;
  MatrixXd delT = invS * covT;
  double DelT = sig2T - (covT.transpose() * invS * covT)(0,0);
  MatrixXd Cxt = (MatrixXd(p,p+q) << Cxx,Cxy).finished() * delT;
  double Ctt = (delT.transpose() * Cxxyy * delT)(0,0) + DelT;
  
  // // Expected sample cov matrix for To
  MatrixXd covTo(p+q,1);
  covTo << sig2To*P_Yosc,
        0*C;
  MatrixXd delTo = invS * covTo;
  double DelTo = sig2To - (covTo.transpose() * invS * covTo)(0,0);
  MatrixXd Cxto = (MatrixXd(p,p+q) << Cxx,Cxy).finished() * delTo;
  double Ctoto = (delTo.transpose() * Cxxyy * delTo)(0,0) + DelTo;
  
  double Ctto = 0 - (covT.transpose() * invS * covTo)(0,0) + (covT.transpose() * invS * Cxxyy * invS * covTo)(0,0);
  
  // // For U ----------------------------
  // Expected sample cov matrix for U
  MatrixXd covU(p+q,1);
  // equivalent to rbind(W,B_T*C) or c(W,B_T*C)
  covU << sig2T*B_T*W,
        (sig2T*B_T*B_T+sig2H)*C;
  MatrixXd delU = invS * covU;
  double DelU = sig2T*B_T*B_T+sig2H - (covU.transpose() * invS * covU)(0,0);
  MatrixXd Cyu = (MatrixXd(q,p+q) << Cxy.transpose(),Cyy).finished() * delU;
  double Cuu = (delU.transpose() * Cxxyy * delU)(0,0) + DelU;
  
  // Expected sample cov matrix for Uo
  MatrixXd covUo(p+q,1);
  covUo << 0*W,
        sig2Uo * P_Xosc;
  MatrixXd delUo = invS * covUo;
  double DelUo = sig2Uo - (covUo.transpose() * invS * covUo)(0,0);
  MatrixXd Cyuo = (MatrixXd(q,p+q) << Cxy.transpose(),Cyy).finished() * delUo;
  double Cuouo = (delUo.transpose() * Cxxyy * delUo)(0,0) + DelUo;
  
  double Cuuo = 0 - (covU.transpose() * invS * covUo)(0,0) + (covU.transpose() * invS * Cxxyy * invS * covUo)(0,0);
  
  // // For E, F and H -------------------
  MatrixXd SX0(p+q,p);MatrixXd SY0(p+q,q);MatrixXd SH0(p+q,1);
  SX0 << sig2X*MatrixXd::Identity(p,p),
          MatrixXd::Zero(q,p);
  SY0 << MatrixXd::Zero(p,q),
          sig2Y*MatrixXd::Identity(q,q);
  SH0 << MatrixXd::Zero(p,1),
          sig2H*C;
          
  // Can be done faster in some way...        
  MatrixXd Cee = sig2X*MatrixXd::Identity(p,p) - 
                  SX0.transpose() * invS * SX0 + SX0.transpose() * invS * Cxxyy * invS * SX0;
  MatrixXd Cff = sig2Y*MatrixXd::Identity(p,p) - 
                  SY0.transpose() * invS * SY0 + SY0.transpose() * invS * Cxxyy * invS * SY0;
  MatrixXd Chh = sig2H*MatrixXd::Identity(1,1) - 
                  SH0.transpose() * invS * SH0 + SH0.transpose() * invS * Cxxyy * invS * SH0;
  
  double Cut = sig2T*B_T - (covU.transpose()*invS*covT)(0,0) + (delU.transpose()*Cxxyy*delT)(0,0);
  
 // ------------- Maximization step is done in R with multiroot()
 
  List ret;
  ret["Cxt"] = Cxt;
  ret["Cxto"] = Cxto;
  ret["Ctt"] = Ctt;
  ret["Ctto"] = Ctto;
  ret["Ctoto"] = Ctoto;
  
  ret["Cyu"] = Cyu;
  ret["Cyuo"] = Cyuo;
  ret["Cuu"] = Cuu;
  ret["Cuuo"] = Cuuo;
  ret["Cuouo"] = Cuouo;
  
  ret["Cut"] = Cut;
  
  ret["diag_Cee"] = Cee.diagonal();
  ret["diag_Cff"] = Cff.diagonal();
  ret["Chh"] = Chh;
  
  return ret;
  
}


// [[Rcpp::export]]
Eigen::MatrixXd simulC(const int N, Eigen::VectorXd W,Eigen::VectorXd C,Eigen::VectorXd P_Yosc,Eigen::VectorXd P_Xosc,double B_T,
                  double sigX,double sigY,double sigH,double sigT,double sigTo,double sigUo)
{
  RNGScope scope;    // ensure RNG gets set/reset
  const int p = W.size();
  const int q = C.size();
  const int pq = max(p,q);
  
  // Simulate LV's with sd's
  NumericVector t1 = rnorm(N);
  NumericVector to1 = rnorm(N);
  NumericVector uo1 = rnorm(N);
  NumericVector h1 = rnorm(N);
  
  NumericMatrix e1(N,p);
  NumericMatrix f1(N,q);
  
  for(int i=0;i<pq;i++){
  if(i<p){
    e1(_,i) = rnorm(N);
  }
  if (i<q) {
    f1(_,i) = rnorm(N);
    }
  }
  
  // convert to VectorXd (Eigen)
  Map<VectorXd> Tt = as<Map<VectorXd> >(t1);
  Map<VectorXd> TYo = as<Map<VectorXd> >(to1);
  Map<VectorXd> UXo = as<Map<VectorXd> >(uo1);
  
  Map<VectorXd> H = as<Map<VectorXd> >(h1);
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
