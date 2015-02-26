#include <iostream>
#include <RcppEigen.h>
#include <Eigen/Dense>
#include <Eigen/Cholesky>
#include <boost/math/tools/roots.hpp>
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(BH)]]

using namespace Eigen;
using namespace std;
using namespace Rcpp;
using namespace boost;
//using Eigen::MatrixXd;                  // variable size matrix, double precision
//using Eigen::VectorXd;                  // variable size vector, double precision

// [[Rcpp::export]]
double loglC(VectorXd W,VectorXd C,VectorXd P_Yosc,VectorXd P_Xosc,double B_T, MatrixXd Dat,
                  double sig2X,double sig2Y,double sig2H,double sig2T,double sig2To,double sig2Uo) 
{
  unsigned int N = Dat.rows();
  unsigned int p = W.size();
  unsigned int q = C.size();
  
  MatrixXd SX;
  SX = sig2T*W*W.transpose() + sig2To*P_Yosc*P_Yosc.transpose();
  for(unsigned int i=0;i<p;i++){SX(i,i) += sig2X;}
  MatrixXd SY;
  SY = sig2T*B_T*B_T*C*C.transpose() + sig2H*C*C.transpose() + sig2Uo*P_Xosc*P_Xosc.transpose();
  for(unsigned int i=0;i<q;i++){SY(i,i) += sig2Y;}
  MatrixXd SXY;
  SXY = SXY = sig2T* B_T*W*C.transpose();
  MatrixXd NewSX(SX.rows()+SY.rows(),SX.cols() + SXY.cols());
  NewSX << SX,SXY,
            SXY.transpose(),SY;
  Eigen::LLT<MatrixXd> lltSX(NewSX);
  MatrixXd invS = lltSX.solve(MatrixXd::Identity(p+q,p+q));
  
  double som2 = (Dat.transpose()*Dat*invS).trace();
  MatrixXd L = lltSX.matrixL();
  MatrixXd Logdiag;
  Logdiag = L.diagonal();
  for(int j=0;j<(p+q);j++){
    Logdiag(j) = 2*log(Logdiag(j));
  }
  
  double Loglik = - 0.5*N*Logdiag.sum() - 0.5 * som2;
  return Loglik;
  
}
// loglC(W,C,PYo,PXo,c(B_T),cbind(X,Y),sigX^2,sigY^2,sigH^2,sigT^2,sigTo^2,sigUo^2)

// [[Rcpp::export]]
MatrixXd simulC(const int N, VectorXd W,VectorXd C,VectorXd P_Yosc,VectorXd P_Xosc,double B_T,
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
  Eigen::Map<Eigen::VectorXd> Tt = as<Eigen::Map<Eigen::VectorXd> >(t1);
  Eigen::Map<Eigen::VectorXd> TYo = as<Eigen::Map<Eigen::VectorXd> >(to1);
  Eigen::Map<Eigen::VectorXd> UXo = as<Eigen::Map<Eigen::VectorXd> >(uo1);
  
  Eigen::Map<Eigen::VectorXd> H = as<Eigen::Map<Eigen::VectorXd> >(h1);
  Eigen::Map<Eigen::MatrixXd> E = as<Eigen::Map<Eigen::MatrixXd> >(e1);
  Eigen::Map<Eigen::MatrixXd> Ff = as<Eigen::Map<Eigen::MatrixXd> >(f1);
  
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

// [[Rcpp::export]]
double findr1(double lambda2, VectorXd Cxt, VectorXd Cxto, double Ctt, double Ctto, double Ctoto)
  {
  
  double alp2 = Ctoto + lambda2;
  
  double Cxto2 = Cxto.squaredNorm();
  double Cxt2 = Cxt.squaredNorm();
  MatrixXd Cxtto2_Mat = Cxto.transpose()*Cxt;
  double Cxtto2 = Cxtto2_Mat(0,0);
  
  double a_rt = Cxto2 - alp2*alp2;
  double b_rt = 2*Ctt*Cxto2 - 2*(Cxtto2*Ctto + Ctt*alp2*alp2 - alp2*Ctto*Ctto);
  double c_rt = Cxto2*Ctt*Ctt - 2*Cxtto2*Ctt*Ctto + Cxt2*Ctto*Ctto - 
    Ctt*Ctt*alp2*alp2 + 2*Ctt*alp2*Ctto*Ctto - Ctto*Ctto*Ctto*Ctto;
  
  double D_rt = b_rt*b_rt - 4*a_rt*c_rt;
  double lambda1 = (-b_rt - sqrt(D_rt))/(2*a_rt);
  
  double alp1 = Ctt+lambda1;
  
  VectorXd Pl = (Cxto*alp1 - Cxt*Ctto) / (alp1*alp2 - Ctto*Ctto);
  VectorXd Wl = (Cxt - Pl*Ctto) / alp1;
  
  double tst1 = Wl.squaredNorm();
  //double tst2 = Pl.squaredNorm();
  
  return (tst1-1);
  
}
// [[Rcpp::export]]
MatrixXd findr2(double lambda2, VectorXd Cxt, VectorXd Cxto, double Ctt, double Ctto, double Ctoto)
  {
  
  double alp2 = Ctoto + lambda2;
  
  double Cxto2 = Cxto.squaredNorm();
  double Cxt2 = Cxt.squaredNorm();
  MatrixXd Cxtto2_Mat = Cxto.transpose()*Cxt;
  double Cxtto2 = Cxtto2_Mat(0,0);
  
  double a_rt = Cxto2 - alp2*alp2;
  double b_rt = 2*Ctt*Cxto2 - 2*(Cxtto2*Ctto + Ctt*alp2*alp2 - alp2*Ctto*Ctto);
  double c_rt = Cxto2*Ctt*Ctt - 2*Cxtto2*Ctt*Ctto + Cxt2*Ctto*Ctto - 
    Ctt*Ctt*alp2*alp2 + 2*Ctt*alp2*Ctto*Ctto - Ctto*Ctto*Ctto*Ctto;
  
  double D_rt = b_rt*b_rt - 4*a_rt*c_rt;
  double lambda1 = (-b_rt - sqrt(D_rt))/(2*a_rt);
  
  double alp1 = Ctt+lambda1;
  
  VectorXd Pl = (Cxto*alp1 - Cxt*Ctto) / (alp1*alp2 - Ctto*Ctto);
  VectorXd Wl = (Cxt - Pl*Ctto) / alp1;
  MatrixXd outp(Wl.rows(),2);
  outp << Wl, Pl;
  return outp;
}

// [[Rcpp::export]]
List EMstepC(VectorXd W,VectorXd C,VectorXd P_Yosc,VectorXd P_Xosc,double B_T, MatrixXd Dat,
              double sig2X,double sig2Y,double sig2H,double sig2T,double sig2To,double sig2Uo,
              double startroot_X, double startroot_Y)
{
  const int N = Dat.rows();
  const int p = W.size();
  const int q = C.size();
  
  VectorXd WC(p+q);
  VectorXd PY0(p+q);
  VectorXd PX0(p+q);
  // equivalent to rbind(W,B_T*C) or c(W,B_T*C) etc
  WC << W,
        B_T*C;
  PY0 << P_Yosc,
         0*P_Xosc;
  PX0 << 0*P_Yosc,
         P_Xosc;
  
  // copied from loglC, calc theoretical cov matrix
  MatrixXd SX;
  SX = sig2T*W*W.transpose() + sig2To*P_Yosc*P_Yosc.transpose();
  for(int i=0;i<p;i++){SX(i,i) += sig2X;}
  MatrixXd SY;
  SY = sig2T*B_T*B_T*C*C.transpose() + sig2H*C*C.transpose() + sig2Uo*P_Xosc*P_Xosc.transpose();
  for(int i=0;i<q;i++){SY(i,i) += sig2Y;}
  MatrixXd SXY;
  SXY = SXY = sig2T* B_T*W*C.transpose();
  MatrixXd NewSX(SX.rows()+SY.rows(),SX.cols() + SXY.cols());
  NewSX << SX,SXY,
            SXY.transpose(),SY;
  Eigen::LLT<MatrixXd> lltSX(NewSX);
  MatrixXd invS = lltSX.solve(MatrixXd::Identity(p+q,p+q));
  
  // Expected sample cov matrix for T and To
  MatrixXd Cxxyy = Dat.transpose()*Dat / N;
  MatrixXd Cxx = Cxxyy.topLeftCorner(p, p);
  MatrixXd Cxy = Cxxyy.topRightCorner(p, q);
  MatrixXd Cyy = Cxxyy.bottomRightCorner(q, q);
  
  MatrixXd Cxt = sig2T * Cxxyy.topLeftCorner(p,p+q) * invS * WC;
  MatrixXd Cxto = sig2To * Cxxyy.topLeftCorner(p,p+q) * invS * PY0;
  MatrixXd Ctto = sig2T * sig2To * (-1 * WC.transpose() * invS * PY0 + 
                  WC.transpose() * invS * Cxxyy * invS * PY0);
  MatrixXd Ctt = sig2T * (MatrixXd::Identity(1,1) - 
                  sig2T*WC.transpose()*invS*WC + sig2T*WC.transpose()*invS*Cxxyy*invS*WC);
  MatrixXd Ctoto = sig2To*(MatrixXd::Identity(1,1) - 
                  sig2To*PY0.transpose()*invS*PY0 + sig2To*PY0.transpose()*invS*Cxxyy*invS*PY0);
  MatrixXd Cyt = sig2T*(Cxxyy.topLeftCorner(p+q,p)*invS*WC);
  
  // Expected sample cov matrix for U and Uo
  MatrixXd Cyu = sig2T * Cxxyy.bottomLeftCorner(q,p+q) * invS * WC*B_T;
  MatrixXd Cyuo = sig2Uo * Cxxyy.bottomLeftCorner(q,p+q) * invS * PX0;
  MatrixXd Cuuo = sig2T * sig2Uo * (-1 * B_T*WC.transpose() * invS * PX0 + 
                  B_T*WC.transpose() * invS * Cxxyy * invS * PX0);
  MatrixXd Cuu = sig2T * (MatrixXd::Identity(1,1) - 
                  sig2T*B_T*WC.transpose()*invS*WC*B_T + sig2T*B_T*WC.transpose()*invS*Cxxyy*invS*WC*B_T);
  MatrixXd Cuouo = sig2Uo*(MatrixXd::Identity(1,1) - 
                  sig2Uo*PX0.transpose()*invS*PX0 + sig2Uo*PX0.transpose()*invS*Cxxyy*invS*PX0);
  
  // Expected sample cov matrix for E F H
  MatrixXd SX0(p+q,p);MatrixXd SY0(p+q,q);MatrixXd SH0(p+q,1);
  SX0 << sig2X*MatrixXd::Identity(p,p),
          MatrixXd::Zero(q,p);
  SY0 << MatrixXd::Zero(p,q),
          sig2Y*MatrixXd::Identity(q,q);
  SH0 << MatrixXd::Zero(p,1),
          sig2H*C;
          
  MatrixXd Cee = sig2X*MatrixXd::Identity(p,p) - 
                  SX0.transpose() * invS * SX0 + SX0.transpose() * invS * Cxxyy * invS * SX0;
  MatrixXd Cff = sig2Y*MatrixXd::Identity(p,p) - 
                  SY0.transpose() * invS * SY0 + SY0.transpose() * invS * Cxxyy * invS * SY0;
  MatrixXd Chh = sig2H*MatrixXd::Identity(1,1) - 
                  SH0.transpose() * invS * SH0 + SH0.transpose() * invS * Cxxyy * invS * SH0;
  
  VectorXd Wl; VectorXd PYo; VectorXd Cl; VectorXd PXo;
  
  Environment env = Environment::global_env();
  Function mroot = env["multiroot2"];
  Function findrt = env["findr1"];
  bool unormX = false; bool unormY = false;
  
  //startroot_X = 1/(P_Yosc.transpose()*P_Yosc)(0,0) * 
  //      (P_Yosc.transpose()*Cxto - P_Yosc.transpose()*W*Ctto - P_Yosc.transpose()*P_Yosc*Ctoto)(0,0);
  
  double lambda_2 = Rcpp::as<double>(mroot(findrt,startroot_X,Cxt,Cxto,Ctt,Ctto,Ctoto));
  if(lambda_2 != 12345){
  MatrixXd WPY = findr2(lambda_2,Cxt,Cxto,Ctt(0,0),Ctto(0,0),Ctoto(0,0));
  PYo = WPY.col(1);
  Wl = WPY.col(0);
  }else{
    PYo = (Cxto*Ctt(0,0) - Cxt*Ctto(0,0)) / (Ctt(0,0)*Ctoto(0,0) - Ctto(0,0)*Ctto(0,0));
    Wl = (Cxt - PYo*Ctto(0,0)) / Ctt(0,0);
    PYo /= PYo.norm();
    Wl /= Wl.norm();
    unormX = true;
  }
  
  //startroot_Y = 1/(P_Xosc.transpose()*P_Xosc)(0,0) * 
  //       (P_Xosc.transpose()*Cyuo - P_Xosc.transpose()*C*Cuuo - P_Xosc.transpose()*P_Xosc*Cuouo)(0,0);
  
  double lambda_4 = Rcpp::as<double>(mroot(findrt,startroot_Y,Cyu,Cyuo,Cuu,Cuuo,Cuouo));
  if(lambda_4 != 12345){
    MatrixXd CPX = findr2(lambda_4,Cyu,Cyuo,Cuu(0,0),Cuuo(0,0),Cuouo(0,0));
    PXo = CPX.col(1);
    Cl = CPX.col(0);
  }else{
    PXo = (Cyuo*Cuu(0,0) - Cyu*Cuuo(0,0)) / (Cuu(0,0)*Cuouo(0,0) - Cuuo(0,0)*Cuuo(0,0));
    Cl = (Cyu - PXo*Cuuo(0,0)) / Cuu(0,0);
    PXo /= PXo.norm();
    Cl /= Cl.norm();
    unormY = true;
  }
  
  VectorXd signew(6);
  signew(0) = (Cee.trace()/p); signew(1) = (Cff.trace()/q); signew(2) = (Chh.trace());
  signew(3) = (Ctt.trace()); signew(4) = (Ctoto.trace()); signew(5) = (Cuouo.trace());
  
  //cout << startroot_X << "   ";
  //cout << startroot_Y << endl;
  
  List ret;
  ret["W"] = Wl;
  ret["PYo"] = PYo;
  ret["C"] = Cl;
  ret["PXo"] = PXo;
  ret["B"] = 2;//Cyt.norm() / Cxt.norm();
  ret["sig"] = signew;
  ret["X_forcednorm"] = unormX;
  ret["Y_forcednorm"] = unormY;
  ret["lambda_2"] = lambda_2;
  ret["lambda_4"] = lambda_4;
  
  /*
  ret["Cxt"]=Cxt;
  ret["Cxto"]=Cxto;
  ret["Ctt"]=Ctt;
  ret["Ctto"]=Ctto;
  ret["Ctoto"]=Ctoto;
  
  ret["Cyu"]=Cyu;
  ret["Cyuo"]=Cyuo;
  ret["Cuu"]=Cuu;
  ret["Cuuo"]=Cuuo;
  ret["Cuouo"]=Cuouo;
  
  ret["Cee"]=Cee;
  ret["Cff"]=Cff;
  ret["Chh"]=Chh;
  */
  
  return ret;
  
}

// [[Rcpp::export]]
List EMC(int maxit, List ret, MatrixXd Dat,double startroot_X,double startroot_Y)
{
  VectorXd loglik(maxit+1);
  VectorXd signw = ret["sig"];
  
  loglik(0) = loglC(ret["W"],ret["C"],ret["PYo"],ret["PXo"],ret["B"],Dat,
                    (signw(0)),(signw(1)),(signw(2)),(signw(3)),(signw(4)),(signw(5)));
  for(int i=0;i<maxit;i++){
  signw = ret["sig"];
  ret = EMstepC(ret["W"],ret["C"],ret["PYo"],ret["PXo"],ret["B"],Dat,
                signw(0),signw(1),signw(2),signw(3),signw(4),signw(5),startroot_X,startroot_Y);
  loglik(i+1) = loglC(ret["W"],ret["C"],ret["PYo"],ret["PXo"],ret["B"],Dat,
                    (signw(0)),(signw(1)),(signw(2)),(signw(3)),(signw(4)),(signw(5)));
  }
  ret["loglik"] = loglik;
  return ret;
}
