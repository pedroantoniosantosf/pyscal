#include "atom.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <stdio.h>
#include "string.h"
#include <chrono>
#include <pybind11/stl.h>
#include <complex>

//functions for atoms
//-------------------------------------------------------------------------------------------------------------------------
Atom::Atom( vector<double> pos, int idd, int typ){

    posx = pos[0];
    posy = pos[1];
    posz = pos[2];
    id = idd;
    type = typ;

    //assign other values - the default ones
    belongsto = -1;
    issolid = 0;
    issurface = 0;
    loc = 0;
    isneighborset = 0;
    n_neighbors = 0;
    lcluster = 0;
    head = -1;


    for (int tn = 0; tn<MAXNUMBEROFNEIGHBORS; tn++){
        neighbors[tn] = -1;
        neighbordist[tn] = -1.0;
        neighborweight[tn] = -1.0;
        facevertices[tn] = -1;
        faceverticenumbers[tn] = -1;
        faceperimeters[tn] = -1.0;
        sij[tn] = -1.0;
        //edgelengths[tn] = -1.0;

    }

    for (int t =0; t < MAXRADFUNCS; t++){
      for (int tn = 0; tn<11; tn++){
          q[t][tn] = -1;
          aq[t][tn] = -1;

          for (int tnn =0; tnn<25; tnn++){
              realq[t][tn][tnn] = -1;
              imgq[t][tn][tnn] = -1;
              arealq[t][tn][tnn] = -1;
              aimgq[t][tn][tnn] = -1;
          }
      }
    }
}

Atom::~Atom(){ }


vector<int> Atom::gneighbors(){
    vector<int> nn;
    nn.reserve(n_neighbors);
    for(int i=0;i<n_neighbors;i++){
        nn.emplace_back(neighbors[i]);
    }
    return nn;
}

int Atom::gnneighbors(){
    return n_neighbors;
}
void Atom::snneighbors(int dd){

}

void Atom::sneighdist(vector<double> dd){
}

vector<double> Atom::gneighdist(){
  vector<double> neighdist;
  for(int i=0; i<n_neighbors; i++){
    neighdist.emplace_back(neighbordist[i]);
  }
  return neighdist;
}

void Atom::ssij(vector<double> dd){
}

vector<double> Atom::gsij(){
  vector<double> ss;
  for(int i=0; i<n_neighbors; i++){
    ss.emplace_back(sij[i]);
  }
  return ss;
}

int Atom::gid(){ return id; }
int Atom::gfrenkelnumber(){ return frenkelnumber; }
void Atom::sfrenkelnumber(int nn){ frenkelnumber=nn; }
void Atom::sid(int idd){ id=idd; }
int Atom::gloc(){ return loc; }
void Atom::sloc(int idd){ loc=idd; }
int Atom::gtype(){ return type; }
void Atom::stype(int idd){ type=idd; }
double Atom::gvolume(){ return volume; }
void Atom::svolume(double vv){ volume = vv; }
double Atom::gcutoff(){ return cutoff; }
void Atom::scutoff(double cc){ cutoff = cc; }
double Atom::gasij(){ return avq6q6; }
void Atom::sasij(double vv){ avq6q6 = vv; }
double Atom::gavgvolume(){ return avgvolume; }
void Atom::savgvolume(double vv){ avgvolume = vv; }
int Atom::gsolid(){ return issolid; }
void Atom::ssolid(int idd){
  if (!((idd == 0) || (idd == 1))){
    throw invalid_argument("surface should be 1 or 0");
  }
  issolid=idd; }

bool Atom::gmask(){ return mask; }
void Atom::smask(bool mm){
  if (!((mm == true) || (mm == false))){
    throw invalid_argument("mask should be true or false");
  }
  mask=mm;
}

int Atom::gstructure(){ return structure; }
void Atom::sstructure(int idd){ structure=idd; }
void Atom::scondition(int idd){ condition=idd; }
int Atom::gcondition(){ return condition; }

vector<vector<double>> Atom::gallq(){
    vector<vector<double>> allq;
    vector<double> dummy;
    for(int n=0; n<MAXRADFUNCS; n++){
      dummy.clear();
      for(int i=0; i<11; i++){
          dummy.emplace_back(q[n][i]);
      }
      allq.emplace_back(dummy);
    }
    return allq;
}

vector<vector<double>> Atom::gallaq(){
  vector<vector<double>> allq;
  vector<double> dummy;
  for(int n=0; n<MAXRADFUNCS; n++){
    dummy.clear();
    for(int i=0; i<11; i++){
        dummy.emplace_back(aq[n][i]);
    }
    allq.emplace_back(dummy);
  }
  return allq;
}

void Atom::sallq(vector<vector<double>> allq){
  for(int n=0; n<MAXRADFUNCS; n++){
    for(int i=0; i<11; i++){
        q[n][i] = allq[n][i];
    }
  }
}

void Atom::sallaq(vector<vector<double>> allaq){
  for(int n=0; n<MAXRADFUNCS; n++){
    for(int i=0; i<11; i++){
        aq[n][i] = allaq[n][i];
    }
  }
}

//aceesss funcs
vector<double> Atom::gx(){
    vector<double> pos;
    pos.emplace_back(posx);
    pos.emplace_back(posy);
    pos.emplace_back(posz);
    return pos;
}

void Atom::sx(vector<double> rls){
    posx = rls[0];
    posy = rls[1];
    posz = rls[2];
}


//structure properties
int Atom::gsurface() {return issurface; }
int Atom::gcluster() {return belongsto; }
void Atom::ssurface( int val) {
  if (!((val == 0) || (val == 1))){
    throw invalid_argument("surface should be 1 or 0");
  }
  issurface = val; }

void Atom::scluster( int val) {
  belongsto = val; }

int Atom::glcluster() {return lcluster; }
void Atom::slcluster( int val) {
  if (!((val == 0) || (val == 1))){
    throw invalid_argument("largest_cluster should be 1 or 0");
  }

  lcluster = val; }


py::dict Atom::gcustom(){

    return custom;
}

void Atom::scustom(py::dict cc){
    custom = cc;
}



//for q vals
double Atom::gq(int qq, int n){ return q[n][qq-2]; }
void Atom::sq(int qq, double qval, int n){ q[n][qq-2] = qval; }

double Atom::gaq(int qq, int n){ return aq[n][qq-2]; }
void Atom::saq(int qq, double qval, int n){ aq[n][qq-2] = qval; }

double Atom::gq_big(int qval, int n, bool averaged){

    if ((qval < 2) || (qval > 12)){
        throw invalid_argument("q value should be between 2-12");
    }
    if(averaged == true) { return gaq(qval, n);}
    else {return gq(qval, n);}
}

void Atom::sq_big(int qval, double val,  int n, bool averaged){

    if ((qval < 2) || (qval > 12)){
        throw invalid_argument("q value should be between 2-12");
    }
    if(averaged == true) { saq(qval, val, n);}
    else { sq(qval, val, n);}
}

//overloaded version which takes a vector
vector<double> Atom::gq_big(vector<int> qval, int n, bool averaged ){
    int d;
    if(averaged == true) {
        vector<double> retvals;
        for(int i=0; i<qval.size(); i++){
            if ((qval[i] < 2) || (qval[i] > 12)){
                throw invalid_argument("q value should be between 2-12");
            }
            retvals.push_back(gaq(qval[i], n));
        }
        return retvals;
    }
    else {
        vector<double> retvals;
        for(int i=0; i<qval.size(); i++){
            if ((qval[i] < 2) || (qval[i] > 12)){
                throw invalid_argument("q value should be between 2-12");
            }
            retvals.push_back(gq(qval[i], n));
        }
        return retvals;
    }
}

//overloaded version which takes a vector
void Atom::sq_big(vector<int> qval, vector<double> setvals, int n, bool averaged){

    if(averaged == true) {

        for(int i=0; i<qval.size(); i++){
            if ((qval[i] < 2) || (qval[i] > 12)){
                throw invalid_argument("q value should be between 2-12");
            }

            saq(qval[i], setvals[i], n);
        }
    }
    else {

        for(int i=0; i<qval.size(); i++){
            if ((qval[i] < 2) || (qval[i] > 12)){
                throw invalid_argument("q value should be between 2-12");
            }
            sq(qval[i], setvals[i], n);
        }
    }
}




//functions to set the neighbors for each atoms
void Atom::sneighbors(vector<int> nns){

    //first reset all neighbors
    for (int i = 0;i<MAXNUMBEROFNEIGHBORS;i++){
        neighbors[i] = NILVALUE;
        neighbordist[i] = -1.0;
    }

    //now assign the neighbors
    for(int i=0; i<nns.size(); i++){
        neighbors[i] = nns[i];
        //auto assign weight to 1
        neighborweight[i] = 1.00;
    }

    n_neighbors = nns.size();
    isneighborset = 1;

}

void Atom::sneighborweights(vector<double> nss){
    for(int i=0; i<nss.size(); i++){
        neighborweight[i] = nss[i];
    }
}

vector<double> Atom::gneighborweights(){
    vector <double> rqlms;
    for(int i=0; i<n_neighbors; i++){
        rqlms.emplace_back(neighborweight[i]);
    }
    return rqlms;
}

void Atom::sfacevertices(vector<int> nss){
    for(int i=0; i<nss.size(); i++){
        facevertices[i] = nss[i];
    }
}

vector<int> Atom::gfacevertices(){
    vector <int> rqlms;
    for(int i=0; i<n_neighbors; i++){
        rqlms.emplace_back(facevertices[i]);
    }
    return rqlms;
}

void Atom::svertex_numbers(vector<int> nss){
    vertex_numbers = nss;
}

vector<int> Atom::gvertex_numbers(){

    //loop over face vertices
    return vertex_numbers;
}

void Atom::svertex_vectors(vector<double> nss){
    vertex_vectors = nss;
}

vector<double> Atom::gvertex_vectors(){
    return vertex_vectors;
}

void Atom::sfaceperimeters(vector<double> nss){
    for(int i=0; i<nss.size(); i++){
        faceperimeters[i] = nss[i];
    }
}

vector<double> Atom::gfaceperimeters(){
    vector <double> rqlms;
    for(int i=0; i<n_neighbors; i++){
        rqlms.emplace_back(faceperimeters[i]);
    }
    return rqlms;
}

void Atom::sedgelengths(vector<vector<double>> nss){
    edgelengths.clear();
    edgelengths = nss;
}

vector<vector<double>> Atom::gedgelengths(){
    return edgelengths;
}

vector<int> Atom::gvorovector(){
    vector<int> voro;
    voro.emplace_back(n3);
    voro.emplace_back(n4);
    voro.emplace_back(n5);
    voro.emplace_back(n6);
    return voro;
}

void Atom::svorovector(vector<int> voro){

    n3 = voro[0];
    n4 = voro[1];
    n5 = voro[2];
    n6 = voro[3];
}

double Atom::gangular(){
    return angular;
}

void Atom::sangular(double dd){
    angular = dd;
}

double Atom::gavgangular(){
    return avg_angular;
}

void Atom::savgangular(double dd){
    avg_angular = dd;
}


vector<int> Atom::gchiparams(){
  return chiparams;
}

void Atom::schiparams(vector<int> nns){
  chiparams.clear();
  chiparams = nns;
}

double Atom::gdisorder(){
    return disorder;
}

void Atom::sdisorder(double dd){
    disorder = dd;
}

double Atom::gavgdisorder(){
    return avgdisorder;
}

void Atom::savgdisorder(double dd){
    avgdisorder = dd;
}

vector<double> Atom::gsro(){
    return sro;
}

void Atom::ssro(vector<double> dd){
    sro = dd;
}

vector<complex<double>> Atom::get_qcomps(int qq, int n, bool averaged){

  vector<complex<double>> qlms;
  qlms.reserve(2*qq+1);
  if (!averaged){
    for(int i=0;i<(2*qq+1);i++){
        complex<double> lmval(realq[n][qq-2][i], imgq[n][qq-2][i]);
        qlms.emplace_back(lmval);
    }
  }
  else{
    for(int i=0;i<(2*qq+1);i++){
        complex<double> lmval(arealq[n][qq-2][i], aimgq[n][qq-2][i]);
        qlms.emplace_back(lmval);
    }
  }
  return qlms;
}


//NEW ACE methods

void Atom::calculate_cheb(double x){
  //set the first two
  chebs[0] = 1.00;
  chebs[1] = x;

  for(int i=2; i<nk; i++){
    chebs[i] = 2*x*chebs[i-1] - chebs[i-2];
  }
}

void Atom::calculate_gks(double r){
  //lmb is lambda the exponential decay
  //r is the interatomic distance
  double factor = 1.0 + cos(PI*(r/cutoff));
  gks[0] = 1.00;
  gks[1] = factor;

  //scaled distances
  double num = (exp(-1*lmb*((r/cutoff) - 1)) - 1.)/(exp(lmb) - 1.);
  double x = 1.0 - 2.0*num;

  //now calculate chebs
  calculate_cheb(x);

  for(int i=2; i<nk; i++){
    gks[i] = 0.25*factor*(1.0 - chebs[i-1]);
  }
}

vector<double> Atom::gcheb(){
  vector<double> rs;
  for(int i=0; i<nk; i++){
    rs.emplace_back(chebs[i]);
  }
  return rs;
}

void Atom::scheb(vector<double> rs){
  for(int i=0; i<nk; i++){
    chebs[i] = rs[i];
  }
}

vector<double> Atom::ggks(){
  vector<double> rs;
  for(int i=0; i<nk; i++){
    rs.emplace_back(gks[i]);
  }
  return rs;
}

void Atom::sgks(vector<double> rs){
  for(int i=0; i<nk; i++){
    gks[i] = rs[i];
  }
}
