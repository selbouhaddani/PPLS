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

/*
//' @export
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
*/

// ---------------- This is one of the two CORE functions
// ---------------- It evaluates the likelihood at each step

//' Title here
//'
//' @param X Numeric matrix.
//' @return ...
// [[Rcpp::export]]
List loglC(Eigen::MatrixXd W,Eigen::MatrixXd C,Eigen::MatrixXd P_Yosc,Eigen::MatrixXd P_Xosc,Eigen::MatrixXd B_T, Eigen::MatrixXd Dat,
             double sigX,double sigY,double sigH,Eigen::MatrixXd sigT,Eigen::MatrixXd sigTo,Eigen::MatrixXd sigUo,
             double c1 = 0, double c2 = 0, double c3 = 0)
{
  double sig2X = sigX*sigX;
  double sig2Y = sigY*sigY;
  MatrixXd sig2H = sigH*sigH*MatrixXd::Identity(1,1);
  MatrixXd sig2T = sigT*sigT;
  MatrixXd sig2To = sigTo*sigTo;
  MatrixXd sig2Uo = sigUo*sigUo;

  int N = Dat.rows(); // number of samples, typically of order 1e2 - 1e3
  int p = W.rows();   // number of X variables, can be 1e1 to 1e4
  int q = C.rows();   // number of Y variables, say can be 1e1 to 1e4
  MatrixXd invS(p+q,p+q);
  // ----------- initialize, fill and invert theoretical cov of cbind(X,Y)
  if(c1 == 0 && c2 == 0 && c3 == 0){
    MatrixXd SX; // LARGE will be p times p
    SX = W*sig2T*W.transpose() + P_Yosc*sig2To*P_Yosc.transpose(); // rank low
    for(int i=0;i<p;i++){SX(i,i) += sig2X;} // add double to diagonal

    MatrixXd SY; // LARGE will be q times q
    SY = C*sig2T*B_T*B_T*C.transpose() + C*sig2H*C.transpose() + P_Xosc*sig2Uo*P_Xosc.transpose(); //rank low
    for(int i=0;i<q;i++){SY(i,i) += sig2Y;} // add double to diagonal

    MatrixXd SXY; // LARGE will be p times q
    SXY = SXY = W*sig2T* B_T*C.transpose(); // only rank low

    // initialize NewSX with dimensions (p+q) times (p+q)
    MatrixXd NewSX(p+q,p+q);
    NewSX << SX,SXY,
             SXY.transpose(),SY; // equivalent to matlab [SX SXY ; SXY' SY]
    LLT<MatrixXd> lltSX(NewSX);   // calc cholesky decomposition (Eigen/Cholesky library)

    //Inversion step can be done maybe with "rank one update", in R I would use chol2inv()
    invS = lltSX.solve(MatrixXd::Identity(p+q,p+q)); // Invert using cholesky
    MatrixXd L = lltSX.matrixL();
    // ----- I can calc this directly I think
    MatrixXd Logdiag = L.diagonal(); // The diagonal of L are the sqrt of the eigenvalues
    for(int j=0;j<(p+q);j++){
      Logdiag(j) = 2*log(Logdiag(j)); // 2*log since we need to square the sqrt(eigenval's)
    }
    double som2 = 0;
    for(int k=0;k<N;k++){som2 += Dat.row(k) * invS * (Dat.row(k)).transpose();}
    double Loglik = - 0.5*N*Logdiag.sum() - 0.5 * som2;

    List ret;
    //cout << Loglik << endl;
    //cout << (invS).trace() << endl;
    ret["value"] = Loglik;
    ret["invS"] = invS;
    return ret;

  } else {
    MatrixXd C1inv = -c1*W*W.transpose();for(int i=0;i<p;i++){C1inv(i,i) += 1/sig2X;}
    MatrixXd C2inv = -c3*C*C.transpose();for(int j=0;j<q;j++){C2inv(j,j) += 1/sig2Y;}

    invS << C1inv , -c2*W*(C.transpose()),
            -c2*C*(W.transpose()) , C2inv;

    double Logdiag = log(sig2X+sig2T(0,0)) + (p-1)*log(sig2X) + log(sig2Y + (c3*sig2Y*sig2Y / (1 - c3 * sig2Y))) + (q-1)*log(sig2Y);

    double som2 = 0;
    for(int k=0;k<N;k++){som2 += Dat.row(k) * invS * (Dat.row(k)).transpose();}
    double Loglik = - 0.5*N*Logdiag - 0.5 * som2;

    List ret;
    //cout << Logdiag << endl;
    //cout << (invS).trace() << endl;

    ret["value"] = Loglik;
    ret["invS"] = invS;
    return ret;
  }
  // -----------

}

//' Title here
//'
//' @param X Numeric matrix.
//' @return ...
// [[Rcpp::export]]
List EMstepC(Eigen::MatrixXd W,Eigen::MatrixXd C, Eigen::MatrixXd P_Yosc, Eigen::MatrixXd P_Xosc,
             Eigen::MatrixXd B_T, Eigen::MatrixXd Dat,double sigX,double sigY,double sigH,
             Eigen::MatrixXd sigT,Eigen::MatrixXd sigTo,Eigen::MatrixXd sigUo,Eigen::MatrixXd invS = Eigen::MatrixXd::Zero(1,1))
{
  double sig2X = sigX*sigX;
  double sig2Y = sigY*sigY;
  MatrixXd sig2H = sigH*sigH*MatrixXd::Identity(1,1);
  MatrixXd sig2T = sigT*sigT;
  MatrixXd sig2To = sigTo*sigTo;
  MatrixXd sig2Uo = sigUo*sigUo;

  int N = Dat.rows();
  int p = W.rows();
  int q = C.rows();
  int a = W.cols();
  int nx = P_Yosc.cols();
  int ny = P_Xosc.cols();

  if( invS(0,0) == 0 ){
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
    invS = lltSX.solve(MatrixXd::Identity(p+q,p+q)); // Invert using cholesky
    // -----------
  }
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

//' Title here
//'
//' @param X Numeric matrix.
//' @return ...
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd simulC(int N, Eigen::MatrixXd W,Eigen::MatrixXd C, Eigen::MatrixXd P_Yosc, Eigen::MatrixXd P_Xosc,Eigen::MatrixXd B_T,
                       double sigX,double sigY,double sigH,Eigen::MatrixXd sigT,Eigen::MatrixXd sigTo,Eigen::MatrixXd sigUo)
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

  for(int i=0;i<p;i++){ e1(_,i) = sigX*rnorm(N); }
  for(int i=0;i<q;i++){ f1(_,i) = sigY*rnorm(N); }
  for(int i=0;i<a;i++){ t1(_,i) = rnorm(N); }
  for(int i=0;i<nx;i++){ to1(_,i) = rnorm(N); }
  for(int i=0;i<ny;i++){ uo1(_,i) = rnorm(N); }
  for(int i=0;i<a;i++){ h1(_,i) = sigH*rnorm(N); }


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

  VectorXd U = Tt * B_T + H;
  MatrixXd X = Tt*W.transpose() + TYo*P_Yosc.transpose() + E;
  MatrixXd Y = U*C.transpose() + UXo*P_Xosc.transpose() + Ff;
  MatrixXd Dat(N,p+q);
  Dat << X, Y;
  return Dat;
}

//' Title here
//'
//' @param X Numeric matrix.
//' @return ...
// [[Rcpp::export]]
double loglC_fast(Eigen::MatrixXd W,Eigen::MatrixXd C,Eigen::MatrixXd X, Eigen::MatrixXd Y,
           double sigX,double sigY,Eigen::VectorXd sig2T,
           Eigen::VectorXd c1, Eigen::VectorXd c2, Eigen::VectorXd c3 ,Eigen::VectorXd Kc)
{
  double sig2X = sigX*sigX;
  double sig2Y = sigY*sigY;
  //VectorXd sig2T = sigT.array()*sigT.array();
  int N = X.rows(); // number of samples, typically of order 1e2 - 1e3
  int p = W.rows();   // number of X variables, can be 1e1 to 1e4
  int q = C.rows();   // number of Y variables, say can be 1e1 to 1e4
  int a = W.cols();

  double Logdiag = (sig2X+sig2T.array()).log().sum() + (p-a)*log(sig2X) + (sig2Y + Kc.array()).log().sum() + (q-a)*log(sig2Y);
  MatrixXd XW = X*W;
  MatrixXd YC = Y*C;
  double traceL = 1/sig2X * X.array().square().sum() + 1/sig2Y * Y.array().square().sum();
  for(int i=0;i<a;i++) traceL += -c1(i)*(XW.col(i)).squaredNorm() - 2*c2(i)*(XW.col(i)).dot(YC.col(i)) - c3(i)*(YC.col(i)).squaredNorm();
  double Loglik = - 0.5*N*Logdiag - 0.5 * traceL;
  return Loglik;
}

//' Title here
//'
//' @param X Numeric matrix.
//' @return ...
// [[Rcpp::export]]
List EMstepC_fast(Eigen::VectorXd W,Eigen::VectorXd C,double B,
                  Eigen::MatrixXd X, Eigen::MatrixXd Y,double sigX,double sigY,double sigH,double sigT,double c1, double c2, double c3)
{
  double sig2X = sigX*sigX;
  double sig2Y = sigY*sigY;
  double sig2H = sigH*sigH;
  double sig2T = sigT*sigT;

  int N = X.rows();
  int p = W.size();
  int q = C.size();

  VectorXd Xw = X*W; VectorXd Yc = Y*C;
  // // // // //  Expectations------------------------------------------------------------------------------
  // // // For T ----------------------------
  VectorXd mu_T = Xw *sig2T*(-c1 + -c2*B + 1/sig2X) + Yc *sig2T*(-c2 + -c3*B + 1/sig2Y*B);
  VectorXd Cxt = X.transpose() * mu_T / N;
  double Ctt = sig2T - sig2T*sig2T*(-c1 - 2*B*c2 - B*B*(c3-1/sig2Y) + 1/sig2X) + mu_T.squaredNorm()/N;


  VectorXd mu_U = Xw *(-sig2T*B*c1 + -c2*(sig2T*B*B+sig2H) + 1/sig2X*B*sig2T) + Y * C *(-c2*B*sig2T + -c3*(sig2T*B*B+sig2H) + 1/sig2Y*(sig2T*B*B+sig2H));
  VectorXd Cyu = Y.transpose() * mu_U / N;
  //double Cuu = (sig2T*B*B+sig2H) - (-(c1-1/sig2X)*sig2T*sig2T*B*B - 2*sig2T*B*(sig2T*B*B+sig2H)*c2 - pow(sig2T*B*B+sig2H,2)*(c3-1/sig2Y)) + mu_U.squaredNorm() / N;

  double Cut = sig2T*B - (-sig2T*sig2T*B*(c1-1/sig2X) - sig2T*sig2T*B*B*c2 - sig2T*(sig2T*B*B+sig2H)*c2 - (sig2T*B*B+sig2H)*sig2T*B*(c3-1/sig2Y)) + mu_U.dot(mu_T)/N;

  double Ceetmp = c1*c1*sig2X*sig2X*(Xw).squaredNorm() + X.array().square().sum() + c2*c2*sig2X*sig2X*(Yc).squaredNorm() -
    2*c1*sig2X*(Xw).squaredNorm() + 2*c1*c2*sig2X*sig2X*(Xw).dot(Yc) - 2*c2*sig2X*(Xw).dot(Yc);
  double Cee = sig2X - (-sig2X*sig2X*(c1) + p*sig2X )/p + ( Ceetmp )/N/p;

  double Cfftmp = c3*c3*sig2Y*sig2Y*(Yc).squaredNorm() + Y.array().square().sum() + c2*c2*sig2Y*sig2Y*(Xw).squaredNorm() -
    2*c3*sig2Y*(Yc).squaredNorm() + 2*c3*c2*sig2Y*sig2Y*(Yc).dot(Xw) - 2*c2*sig2Y*(Yc).dot(Xw);
  double Cff = sig2Y - (-sig2Y*sig2Y*c3 + q*sig2Y )/q + ( Cfftmp )/N/q;

  double Chh = sig2H - (-sig2H*sig2H*(c3-1/sig2Y)) + (-c2*sig2H*Xw - (c3-1/sig2Y)*sig2H*Yc).squaredNorm()/N;

  Vector2d sighat;
  Vector2d siglathat;
  sighat(0) = sqrt(Cee); sighat(1) = sqrt(Cff);
  siglathat(0) = sqrt(Chh); siglathat(1) = sqrt(Ctt);
  List ret;

  ret["W"] = Cxt.normalized();
  ret["C"] = Cyu.normalized();
  ret["B"] = Cut / Ctt;
  ret["sighat"] = sighat;
  ret["siglathat"] = siglathat;
  return ret;

}
