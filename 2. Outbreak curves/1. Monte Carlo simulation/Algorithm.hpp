#include"Data.hpp"
#include<random>
#include<cmath>
#include<utility>
#include<vector>
#include<map>
#include<algorithm>
#include<fstream>
using namespace std;

extern unsigned int year, fEnsembleNumber;
extern string correlation;
extern long double alpha, ilambda, dlambda, flambda, recoveryRate, outbreakThreshold;

class Country; class Product; class Hyperedge;
class Country {public: vector<Hyperedge*> exporterHyperedgeVector, importerHyperedgeVector, countryHyperedgeVector;};
class Product {public: vector<Hyperedge*> productHyperedgeVector;};
class Hyperedge {
public:
  Country *exporter, *importer; Product *product;
  long double weight, alphaWeight; unsigned int degree, state, newState;
};

class Hypergraph {
public:
  map<unsigned int, Country> countryMap; map<unsigned int, Product> productMap; vector<Hyperedge> hyperedgeVector;
  long double outbreakFraction, averageDegree=0, averageDegreeSquared=0, averageAlphaWeight=0;
  
  Hypergraph();
  void SIRDynamics(long double lambda);
};

Hypergraph::Hypergraph() {
  //Construct a trade hypergraph
  ifstream scanFile("UN"+to_string(year)+".txt");
  hyperedgeVector.reserve(1e6);
  unsigned int exporter, importer, product; long double weight; string logRatio;
  while (scanFile >> exporter >> importer >> product >> weight >> logRatio) {
    Hyperedge hyperedge=Hyperedge();
    hyperedge.exporter=&countryMap[exporter];
    hyperedge.importer=&countryMap[importer];
    hyperedge.product=&productMap[product];
    hyperedge.weight=weight;
    hyperedgeVector.emplace_back(hyperedge);
    countryMap[exporter].exporterHyperedgeVector.emplace_back(&hyperedgeVector.back());
    countryMap[exporter].countryHyperedgeVector.emplace_back(&hyperedgeVector.back());
    countryMap[importer].importerHyperedgeVector.emplace_back(&hyperedgeVector.back());
    countryMap[importer].countryHyperedgeVector.emplace_back(&hyperedgeVector.back());
    productMap[product].productHyperedgeVector.emplace_back(&hyperedgeVector.back());
  } scanFile.close();
  
  //Degree
  for (auto &hyperedge : hyperedgeVector) {
    hyperedge.degree=hyperedge.exporter->countryHyperedgeVector.size()+hyperedge.importer->countryHyperedgeVector.size()+hyperedge.product->productHyperedgeVector.size()-3;
    averageDegree+=hyperedge.degree; averageDegreeSquared+=(long double)hyperedge.degree*hyperedge.degree;
  } averageDegree/=hyperedgeVector.size(); averageDegreeSquared/=hyperedgeVector.size();
  
  //Correlation
  if (correlation!="Original") {
    vector<pair<long double, Hyperedge*>> degreeHyperedgeVector; degreeHyperedgeVector.reserve(hyperedgeVector.size());
    vector<long double> weightVector; weightVector.reserve(hyperedgeVector.size());
    for (auto &hyperedge : hyperedgeVector) {
      degreeHyperedgeVector.emplace_back(pair<long double, Hyperedge*>(hyperedge.degree, &hyperedge));
      weightVector.emplace_back(hyperedge.weight);
    }
    random_device rd; mt19937 gen(rd());
    sort(degreeHyperedgeVector.begin(), degreeHyperedgeVector.end());
    if (correlation=="Randomized") {shuffle(weightVector.begin(), weightVector.end(), gen);}
    else if (correlation=="MaximallyPositive") {sort(weightVector.rbegin(), weightVector.rend());}
    else if (correlation=="MaximallyNegative") {sort(weightVector.begin(), weightVector.end());}
    for (auto &degreeHyperedge : degreeHyperedgeVector) {
      degreeHyperedge.second->weight=weightVector.back(); weightVector.pop_back();
    }
  }
  
  //Weight
  for (auto &hyperedge : hyperedgeVector) {
    hyperedge.alphaWeight=pow(hyperedge.weight, alpha);
    averageAlphaWeight+=hyperedge.alphaWeight;
  } averageAlphaWeight/=hyperedgeVector.size();
  
  //SIR Dynamics
  Data outbreakFractionData=Data();
  for (unsigned int ensembleIndex=1; ensembleIndex<=fEnsembleNumber; ++ensembleIndex) {
    unsigned int pointIndex=0;
    for (long double lambda=ilambda; lambda<=flambda; lambda+=dlambda, ++pointIndex) {
      SIRDynamics(lambda);
      outbreakFractionData.Add(pointIndex, outbreakFraction);
    }
    outbreakFractionData.Print("OutbreakFraction"+to_string(year)+"_"+correlation+"_alpha="+to_string(alpha).substr(0,6)+".txt");
  }
}

void Hypergraph::SIRDynamics(long double lambda) {
  random_device rd; mt19937 gen(rd()); uniform_real_distribution<long double> urd(0, 1);
  //Initialization
  for (auto &hyperedge : hyperedgeVector) {hyperedge.state=0; hyperedge.newState=0;}
  Hyperedge* seed=&hyperedgeVector[gen()%hyperedgeVector.size()]; seed->state=1; seed->newState=1;
  //SIR dynamics
  outbreakFraction=0;
  vector<Hyperedge*> infectedHyperedgeVector; infectedHyperedgeVector.reserve(hyperedgeVector.size());
  infectedHyperedgeVector.emplace_back(seed);
  while (infectedHyperedgeVector.size()) {
    for (auto &infectedHyperedge : infectedHyperedgeVector) {
      for (auto &neighborHyperedge : infectedHyperedge->exporter->countryHyperedgeVector) {
        if (neighborHyperedge->state==0 && lambda*recoveryRate*neighborHyperedge->alphaWeight*averageDegree/(averageDegreeSquared*averageAlphaWeight)>urd(gen)) {neighborHyperedge->newState=1;}
      }
      for (auto &neighborHyperedge : infectedHyperedge->importer->countryHyperedgeVector) {
        if (neighborHyperedge->state==0 && lambda*recoveryRate*neighborHyperedge->alphaWeight*averageDegree/(averageDegreeSquared*averageAlphaWeight)>urd(gen)) {neighborHyperedge->newState=1;}
      }
      for (auto &neighborHyperedge : infectedHyperedge->product->productHyperedgeVector) {
        if (neighborHyperedge->state==0 && lambda*recoveryRate*neighborHyperedge->alphaWeight*averageDegree/(averageDegreeSquared*averageAlphaWeight)>urd(gen)) {neighborHyperedge->newState=1;}
      }
    }
    infectedHyperedgeVector.clear();
    for (auto &hyperedge : hyperedgeVector) {
      if (hyperedge.state==1 && recoveryRate>urd(gen)) {hyperedge.newState=2; ++outbreakFraction;}
      hyperedge.state=hyperedge.newState; if (hyperedge.state==1) {infectedHyperedgeVector.emplace_back(&hyperedge);}
    }
  }
  outbreakFraction/=hyperedgeVector.size();
}
