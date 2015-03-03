
##### INSTALL.PACKAGES #####################################
#install.packages("ChemometricsWithR")
#install.packages("kohonen")
#install.packages('mixOmics')
#install.packages('VIM');
#install.packages('fBasics') #skewness
#install.packages('gplots')
#install.packages('Cairo')
#install.packages('bootstrap')
#install.packages('car')
#install.packages('magic')
#install.packages('numDeriv')
#install.packages('doSNOW')
#install.packages('rootSolve')
#install.packages('rgl')
#install.packages('mvtnorm')
#install.packages('Rcpp')
#install.packages('RcppEigen')
#install.packages("inline")
#install.packages('BH')
##### LIBRARY ####
#library('Cairo')
library('car')
library('fBasics')
library("MASS")
#library("VIM")
#library('ChemometricsWithR') #chemo..RData,MASS,pls
#library('mixOmics')
library('gplots')
#library('bootstrap')
#library('boot')
library('magic')
#library('numDeriv')
#library('doSNOW')
library('parallel')
library('rootSolve')
#library('rgl')
library('mvtnorm')
library("inline")
library('Rcpp')
library('RcppEigen')
library('BH')

source('functions.R')
Rcpp::sourceCpp('loglC.cpp')

#########################################################