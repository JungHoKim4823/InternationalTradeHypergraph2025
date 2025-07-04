#include<random>
#include<cmath>
#include<utility>
#include<vector>
#include<map>
#include<algorithm>
#include<fstream>
using namespace std;

extern unsigned int year, fEnsembleNumber;
extern long double alpha, lambda, recoveryRate, outbreakThreshold, binSize, scale;

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
  long double outbreakFraction, averageDegree=0, averageSquaredDegree=0, maximumWeight=0, averageWeight=0, averageAlphaWeight=0;
  
  vector<unsigned int> collapsedHyperedgeNumberVector;
  vector<vector<pair<unsigned int, unsigned int>>> collapsedFractionVSWeightVector;
  vector<pair<long double, unsigned int>> collapsedFractionOfNormalHyperedgeNeighborsVector, collapsedFractionOfNormalHyperedgeOneNeighborsVector, collapsedFractionOfNormalHyperedgeTwoNeighborsVector, collapsedFractionOfNormalHyperedgeThreeNeighborsVector;
  vector<pair<long double, unsigned int>> collapsedFractionOfCollapsedHyperedgeNeighborsVector, collapsedFractionOfCollapsedHyperedgeOneNeighborsVector, collapsedFractionOfCollapsedHyperedgeTwoNeighborsVector, collapsedFractionOfCollapsedHyperedgeThreeNeighborsVector;
  vector<vector<pair<long double, unsigned int>>> collapsedFractionOfNormalHyperedgeNeighborsVSWeightVector, collapsedFractionOfNormalHyperedgeOneNeighborsVSWeightVector, collapsedFractionOfNormalHyperedgeTwoNeighborsVSWeightVector, collapsedFractionOfNormalHyperedgeThreeNeighborsVSWeightVector;
  vector<vector<pair<long double, unsigned int>>> collapsedFractionOfCollapsedHyperedgeNeighborsVSWeightVector, collapsedFractionOfCollapsedHyperedgeOneNeighborsVSWeightVector, collapsedFractionOfCollapsedHyperedgeTwoNeighborsVSWeightVector, collapsedFractionOfCollapsedHyperedgeThreeNeighborsVSWeightVector;
  
  Hypergraph();
  void SIRDynamics(long double lambda);
  void Statistics();
  void Print(unsigned int ensembleIndex);
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
    averageDegree+=hyperedge.degree; averageSquaredDegree+=(long double)hyperedge.degree*hyperedge.degree;
  } averageDegree/=hyperedgeVector.size(); averageSquaredDegree/=hyperedgeVector.size();
  
  //Weight
  for (auto &hyperedge : hyperedgeVector) {
    if (maximumWeight<hyperedge.weight) {maximumWeight=hyperedge.weight;}
    averageWeight+=hyperedge.weight;
    hyperedge.alphaWeight=pow(hyperedge.weight, alpha); averageAlphaWeight+=hyperedge.alphaWeight;
  } averageWeight/=hyperedgeVector.size(); averageAlphaWeight/=hyperedgeVector.size();
  
  //Initialization
  collapsedHyperedgeNumberVector.reserve(fEnsembleNumber);
  collapsedFractionVSWeightVector.reserve(fEnsembleNumber);
  collapsedFractionOfNormalHyperedgeNeighborsVector.reserve(fEnsembleNumber); collapsedFractionOfNormalHyperedgeOneNeighborsVector.reserve(fEnsembleNumber); collapsedFractionOfNormalHyperedgeTwoNeighborsVector.reserve(fEnsembleNumber); collapsedFractionOfNormalHyperedgeThreeNeighborsVector.reserve(fEnsembleNumber);
  collapsedFractionOfCollapsedHyperedgeNeighborsVector.reserve(fEnsembleNumber); collapsedFractionOfCollapsedHyperedgeOneNeighborsVector.reserve(fEnsembleNumber); collapsedFractionOfCollapsedHyperedgeTwoNeighborsVector.reserve(fEnsembleNumber); collapsedFractionOfCollapsedHyperedgeThreeNeighborsVector.reserve(fEnsembleNumber);
  collapsedFractionOfNormalHyperedgeNeighborsVSWeightVector.reserve(fEnsembleNumber); collapsedFractionOfNormalHyperedgeOneNeighborsVSWeightVector.reserve(fEnsembleNumber); collapsedFractionOfNormalHyperedgeTwoNeighborsVSWeightVector.reserve(fEnsembleNumber); collapsedFractionOfNormalHyperedgeThreeNeighborsVSWeightVector.reserve(fEnsembleNumber);
  collapsedFractionOfCollapsedHyperedgeNeighborsVSWeightVector.reserve(fEnsembleNumber); collapsedFractionOfCollapsedHyperedgeOneNeighborsVSWeightVector.reserve(fEnsembleNumber); collapsedFractionOfCollapsedHyperedgeTwoNeighborsVSWeightVector.reserve(fEnsembleNumber); collapsedFractionOfCollapsedHyperedgeThreeNeighborsVSWeightVector.reserve(fEnsembleNumber);
  
  //SIR Dynamics
  for (unsigned int ensembleIndex=1; ensembleIndex<=fEnsembleNumber; ++ensembleIndex) {
    SIRDynamics(lambda);
    if (outbreakFraction>outbreakThreshold) {Statistics(); Print(ensembleIndex);}
  }
  Print(fEnsembleNumber);
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
        if (neighborHyperedge->state==0 && lambda*recoveryRate*neighborHyperedge->alphaWeight*averageDegree/(averageSquaredDegree*averageAlphaWeight)>urd(gen)) {neighborHyperedge->newState=1;}
      }
      for (auto &neighborHyperedge : infectedHyperedge->importer->countryHyperedgeVector) {
        if (neighborHyperedge->state==0 && lambda*recoveryRate*neighborHyperedge->alphaWeight*averageDegree/(averageSquaredDegree*averageAlphaWeight)>urd(gen)) {neighborHyperedge->newState=1;}
      }
      for (auto &neighborHyperedge : infectedHyperedge->product->productHyperedgeVector) {
        if (neighborHyperedge->state==0 && lambda*recoveryRate*neighborHyperedge->alphaWeight*averageDegree/(averageSquaredDegree*averageAlphaWeight)>urd(gen)) {neighborHyperedge->newState=1;}
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

void Hypergraph::Statistics() {
  //The number of collapsed hyperedges
  unsigned int collapsedHyperedgeNumber=0;
  for (auto &hyperedge : hyperedgeVector) {if (hyperedge.state==2) {++collapsedHyperedgeNumber;}}
  collapsedHyperedgeNumberVector.emplace_back(collapsedHyperedgeNumber);
  
  //Collapsed fraction VS Weight
  unsigned int pointNumber=log(maximumWeight/(scale*averageWeight))/log(binSize)+1;
  collapsedFractionVSWeightVector.emplace_back(vector<pair<unsigned int, unsigned int>>(pointNumber));
  for (auto &hyperedge : hyperedgeVector) {
    unsigned int binIndex=log(hyperedge.weight/(scale*averageWeight))/log(binSize);
    if (hyperedge.state==2) {++collapsedFractionVSWeightVector.back()[binIndex].first;}
    ++collapsedFractionVSWeightVector.back()[binIndex].second;
  }
  
  //Neighbor related quantities
  collapsedFractionOfNormalHyperedgeNeighborsVector.emplace_back(pair<long double, unsigned int>()); collapsedFractionOfNormalHyperedgeOneNeighborsVector.emplace_back(pair<long double, unsigned int>()); collapsedFractionOfNormalHyperedgeTwoNeighborsVector.emplace_back(pair<long double, unsigned int>()); collapsedFractionOfNormalHyperedgeThreeNeighborsVector.emplace_back(pair<long double, unsigned int>());
  collapsedFractionOfCollapsedHyperedgeNeighborsVector.emplace_back(pair<long double, unsigned int>()); collapsedFractionOfCollapsedHyperedgeOneNeighborsVector.emplace_back(pair<long double, unsigned int>()); collapsedFractionOfCollapsedHyperedgeTwoNeighborsVector.emplace_back(pair<long double, unsigned int>()); collapsedFractionOfCollapsedHyperedgeThreeNeighborsVector.emplace_back(pair<long double, unsigned int>());
  collapsedFractionOfNormalHyperedgeNeighborsVSWeightVector.emplace_back(vector<pair<long double, unsigned int>>(pointNumber)); collapsedFractionOfNormalHyperedgeOneNeighborsVSWeightVector.emplace_back(vector<pair<long double, unsigned int>>(pointNumber)); collapsedFractionOfNormalHyperedgeTwoNeighborsVSWeightVector.emplace_back(vector<pair<long double, unsigned int>>(pointNumber)); collapsedFractionOfNormalHyperedgeThreeNeighborsVSWeightVector.emplace_back(vector<pair<long double, unsigned int>>(pointNumber));
  collapsedFractionOfCollapsedHyperedgeNeighborsVSWeightVector.emplace_back(vector<pair<long double, unsigned int>>(pointNumber)); collapsedFractionOfCollapsedHyperedgeOneNeighborsVSWeightVector.emplace_back(vector<pair<long double, unsigned int>>(pointNumber)); collapsedFractionOfCollapsedHyperedgeTwoNeighborsVSWeightVector.emplace_back(vector<pair<long double, unsigned int>>(pointNumber)); collapsedFractionOfCollapsedHyperedgeThreeNeighborsVSWeightVector.emplace_back(vector<pair<long double, unsigned int>>(pointNumber));
  for (auto &hyperedge : hyperedgeVector) {
    //Find neighbors
    vector<Hyperedge*> neighborVector, oneNeighborVector, twoNeighborVector, threeNeighborVector;
    neighborVector.reserve(1e5); oneNeighborVector.reserve(1e5); twoNeighborVector.reserve(1e5); threeNeighborVector.reserve(1);
    for (auto &neighbor : hyperedge.exporter->exporterHyperedgeVector) {
      if (neighbor->importer==hyperedge.importer && neighbor->product==hyperedge.product) {continue;}
      else if (neighbor->importer==hyperedge.importer) {neighborVector.emplace_back(neighbor); twoNeighborVector.emplace_back(neighbor);}
      else if (neighbor->product==hyperedge.product) {neighborVector.emplace_back(neighbor); twoNeighborVector.emplace_back(neighbor);}
      else {neighborVector.emplace_back(neighbor); oneNeighborVector.emplace_back(neighbor);}
    }
    for (auto &neighbor : hyperedge.exporter->importerHyperedgeVector) {
      if (neighbor->exporter==hyperedge.importer && neighbor->product==hyperedge.product) {neighborVector.emplace_back(neighbor); threeNeighborVector.emplace_back(neighbor);}
      else if (neighbor->exporter==hyperedge.importer) {neighborVector.emplace_back(neighbor); twoNeighborVector.emplace_back(neighbor);}
      else if (neighbor->product==hyperedge.product) {neighborVector.emplace_back(neighbor); twoNeighborVector.emplace_back(neighbor);}
      else {neighborVector.emplace_back(neighbor); oneNeighborVector.emplace_back(neighbor);}
    }
    for (auto &neighbor : hyperedge.importer->importerHyperedgeVector) {
      if (neighbor->exporter==hyperedge.exporter && neighbor->product==hyperedge.product) {continue;}
      else if (neighbor->exporter==hyperedge.exporter) {continue;}
      else if (neighbor->product==hyperedge.product) {neighborVector.emplace_back(neighbor); twoNeighborVector.emplace_back(neighbor);}
      else {neighborVector.emplace_back(neighbor); oneNeighborVector.emplace_back(neighbor);}
    }
    for (auto &neighbor : hyperedge.importer->exporterHyperedgeVector) {
      if (neighbor->importer==hyperedge.exporter && neighbor->product==hyperedge.product) {continue;}
      else if (neighbor->importer==hyperedge.exporter) {continue;}
      else if (neighbor->product==hyperedge.product) {neighborVector.emplace_back(neighbor); twoNeighborVector.emplace_back(neighbor);}
      else {neighborVector.emplace_back(neighbor); oneNeighborVector.emplace_back(neighbor);}
    }
    for (auto &neighbor : hyperedge.product->productHyperedgeVector) {
      if (neighbor->exporter==hyperedge.exporter && neighbor->importer==hyperedge.importer) {continue;}
      else if (neighbor->exporter==hyperedge.importer && neighbor->importer==hyperedge.exporter) {continue;}
      else if (neighbor->exporter==hyperedge.exporter) {continue;}
      else if (neighbor->exporter==hyperedge.importer) {continue;}
      else if (neighbor->importer==hyperedge.importer) {continue;}
      else if (neighbor->importer==hyperedge.exporter) {continue;}
      else {neighborVector.emplace_back(neighbor); oneNeighborVector.emplace_back(neighbor);}
    }
    
    //Collapsed fraction of hyperedge neighbors
    long double collapsedFractionOfHyperedgeNeighbors=0, collapsedFractionOfHyperedgeOneNeighbors=0, collapsedFractionOfHyperedgeTwoNeighbors=0, collapsedFractionOfHyperedgeThreeNeighbors=0;
    if (neighborVector.size()) {
      for (auto &neighbor : neighborVector) {
        if (neighbor->state==2) {++collapsedFractionOfHyperedgeNeighbors;}
      } collapsedFractionOfHyperedgeNeighbors/=neighborVector.size();
    }
    if (oneNeighborVector.size()) {
      for (auto &oneNeighbor : oneNeighborVector) {
        if (oneNeighbor->state==2) {++collapsedFractionOfHyperedgeOneNeighbors;}
      } collapsedFractionOfHyperedgeOneNeighbors/=oneNeighborVector.size();
    }
    if (twoNeighborVector.size()) {
      for (auto &twoNeighbor : twoNeighborVector) {
        if (twoNeighbor->state==2) {++collapsedFractionOfHyperedgeTwoNeighbors;}
      } collapsedFractionOfHyperedgeTwoNeighbors/=twoNeighborVector.size();
    }
    if (threeNeighborVector.size()) {
      for (auto &threeNeighbor : threeNeighborVector) {
        if (threeNeighbor->state==2) {++collapsedFractionOfHyperedgeThreeNeighbors;}
      } collapsedFractionOfHyperedgeThreeNeighbors/=threeNeighborVector.size();
    }
    if (hyperedge.state==2) {
      if (neighborVector.size()) {
        collapsedFractionOfCollapsedHyperedgeNeighborsVector.back().first+=collapsedFractionOfHyperedgeNeighbors;
        ++collapsedFractionOfCollapsedHyperedgeNeighborsVector.back().second;
      }
      if (oneNeighborVector.size()) {
        collapsedFractionOfCollapsedHyperedgeOneNeighborsVector.back().first+=collapsedFractionOfHyperedgeOneNeighbors;
        ++collapsedFractionOfCollapsedHyperedgeOneNeighborsVector.back().second;
      }
      if (twoNeighborVector.size()) {
        collapsedFractionOfCollapsedHyperedgeTwoNeighborsVector.back().first+=collapsedFractionOfHyperedgeTwoNeighbors;
        ++collapsedFractionOfCollapsedHyperedgeTwoNeighborsVector.back().second;
      }
      if (threeNeighborVector.size()) {
        collapsedFractionOfCollapsedHyperedgeThreeNeighborsVector.back().first+=collapsedFractionOfHyperedgeThreeNeighbors;
        ++collapsedFractionOfCollapsedHyperedgeThreeNeighborsVector.back().second;
      }
    } else {
      if (neighborVector.size()) {
        collapsedFractionOfNormalHyperedgeNeighborsVector.back().first+=collapsedFractionOfHyperedgeNeighbors;
        ++collapsedFractionOfNormalHyperedgeNeighborsVector.back().second;
      }
      if (oneNeighborVector.size()) {
        collapsedFractionOfNormalHyperedgeOneNeighborsVector.back().first+=collapsedFractionOfHyperedgeOneNeighbors;
        ++collapsedFractionOfNormalHyperedgeOneNeighborsVector.back().second;
      }
      if (twoNeighborVector.size()) {
        collapsedFractionOfNormalHyperedgeTwoNeighborsVector.back().first+=collapsedFractionOfHyperedgeTwoNeighbors;
        ++collapsedFractionOfNormalHyperedgeTwoNeighborsVector.back().second;
      }
      if (threeNeighborVector.size()) {
        collapsedFractionOfNormalHyperedgeThreeNeighborsVector.back().first+=collapsedFractionOfHyperedgeThreeNeighbors;
        ++collapsedFractionOfNormalHyperedgeThreeNeighborsVector.back().second;
      }
    }
    
    //Collapsed fraction of hyperedge neighbors VS weight
    unsigned int binIndex=log(hyperedge.weight/(scale*averageWeight))/log(binSize);
    if (hyperedge.state==2) {
      collapsedFractionOfCollapsedHyperedgeNeighborsVSWeightVector.back()[binIndex].first+=collapsedFractionOfHyperedgeNeighbors;
      ++collapsedFractionOfCollapsedHyperedgeNeighborsVSWeightVector.back()[binIndex].second;
      collapsedFractionOfCollapsedHyperedgeOneNeighborsVSWeightVector.back()[binIndex].first+=collapsedFractionOfHyperedgeOneNeighbors;
      ++collapsedFractionOfCollapsedHyperedgeOneNeighborsVSWeightVector.back()[binIndex].second;
      collapsedFractionOfCollapsedHyperedgeTwoNeighborsVSWeightVector.back()[binIndex].first+=collapsedFractionOfHyperedgeTwoNeighbors;
      ++collapsedFractionOfCollapsedHyperedgeTwoNeighborsVSWeightVector.back()[binIndex].second;
      collapsedFractionOfCollapsedHyperedgeThreeNeighborsVSWeightVector.back()[binIndex].first+=collapsedFractionOfHyperedgeThreeNeighbors;
      ++collapsedFractionOfCollapsedHyperedgeThreeNeighborsVSWeightVector.back()[binIndex].second;
    } else {
      collapsedFractionOfNormalHyperedgeNeighborsVSWeightVector.back()[binIndex].first+=collapsedFractionOfHyperedgeNeighbors;
      ++collapsedFractionOfNormalHyperedgeNeighborsVSWeightVector.back()[binIndex].second;
      collapsedFractionOfNormalHyperedgeOneNeighborsVSWeightVector.back()[binIndex].first+=collapsedFractionOfHyperedgeOneNeighbors;
      ++collapsedFractionOfNormalHyperedgeOneNeighborsVSWeightVector.back()[binIndex].second;
      collapsedFractionOfNormalHyperedgeTwoNeighborsVSWeightVector.back()[binIndex].first+=collapsedFractionOfHyperedgeTwoNeighbors;
      ++collapsedFractionOfNormalHyperedgeTwoNeighborsVSWeightVector.back()[binIndex].second;
      collapsedFractionOfNormalHyperedgeThreeNeighborsVSWeightVector.back()[binIndex].first+=collapsedFractionOfHyperedgeThreeNeighbors;
      ++collapsedFractionOfNormalHyperedgeThreeNeighborsVSWeightVector.back()[binIndex].second;
    }
  }
}

void Hypergraph::Print(unsigned int ensembleIndex) {
  //The number of collapsed hyperedges
  long double sumCollapsedHyperedgeNumber=0, sumSquaredCollapsedHyperedgeNumber=0;
  for (auto &collapsedHyperedgeNumber : collapsedHyperedgeNumberVector) {sumCollapsedHyperedgeNumber+=collapsedHyperedgeNumber; sumSquaredCollapsedHyperedgeNumber+=(long double)collapsedHyperedgeNumber*collapsedHyperedgeNumber;}
  ofstream printFileTheNumberOfCollapsedHyperedgesFirstLine("1. The number of collapsed hyperedges/TheNumberOfCollapsedHyperedges0000.txt");
  printFileTheNumberOfCollapsedHyperedgesFirstLine << "Year" << "\t" << "Ensembles" << "\t" << "Samples" << "\t" << "Hyperedges" << "\t" << "Collapsed" << "\t" << "<X^2>" << "\n";
  printFileTheNumberOfCollapsedHyperedgesFirstLine.close();
  ofstream printFileTheNumberOfCollapsedHyperedges("1. The number of collapsed hyperedges/TheNumberOfCollapsedHyperedges"+to_string(year)+".txt");
  printFileTheNumberOfCollapsedHyperedges << year << "\t" << ensembleIndex << "\t" << collapsedHyperedgeNumberVector.size() << "\t" << hyperedgeVector.size() << "\t" << sumCollapsedHyperedgeNumber/collapsedHyperedgeNumberVector.size() << "\t" << sumSquaredCollapsedHyperedgeNumber/collapsedHyperedgeNumberVector.size() << "\n";
  printFileTheNumberOfCollapsedHyperedges.close();
  
  //Collapsed fraction VS Weight
  ofstream printFileCollapsedFractionVSWeight("2. Collapsed fraction VS weight/CollapsedFractionVSWeight"+to_string(year)+".txt");
  ofstream printFileCollapsedFractionVSWeightFitting("2. Collapsed fraction VS weight/CollapsedFractionVSWeightFitting"+to_string(year)+".txt");
  printFileCollapsedFractionVSWeight << "Ensembles" << "\t" << "Samples" << "\t" <<  "w/<w>" << "\t" << "Collapsed fraction" << "\t" << "<X^2>" << "\n";
  printFileCollapsedFractionVSWeightFitting << "Ensembles" << "\t" << "Samples" << "\t" <<  "w/<w>" << "\t" << "Collapsed fraction" << "\t" << "<X^2>" << "\n";
  for (unsigned int binIndex=0; binIndex<collapsedFractionVSWeightVector.front().size(); ++binIndex) {
    long double sumCollapsedFractionVSWeight=0, sumSquaredCollapsedFractionVSWeight=0;
    for (auto &collapsedFractionVSWeight : collapsedFractionVSWeightVector) {
      sumCollapsedFractionVSWeight+=(long double)collapsedFractionVSWeight[binIndex].first/collapsedFractionVSWeight[binIndex].second;
      sumSquaredCollapsedFractionVSWeight+=((long double)collapsedFractionVSWeight[binIndex].first/collapsedFractionVSWeight[binIndex].second)*((long double)collapsedFractionVSWeight[binIndex].first/collapsedFractionVSWeight[binIndex].second);
    }
    printFileCollapsedFractionVSWeight << ensembleIndex << "\t" << collapsedFractionVSWeightVector.size() << "\t" << pow(binSize, binIndex+0.5)*scale << "\t" << sumCollapsedFractionVSWeight/collapsedFractionVSWeightVector.size() << "\t" << sumSquaredCollapsedFractionVSWeight/collapsedFractionVSWeightVector.size() << "\n";
    if (binIndex>=9 && binIndex<=16) {
      printFileCollapsedFractionVSWeightFitting << ensembleIndex << "\t" << collapsedFractionVSWeightVector.size() << "\t" << pow(binSize, binIndex+0.5)*scale << "\t" << sumCollapsedFractionVSWeight/collapsedFractionVSWeightVector.size() << "\t" << sumSquaredCollapsedFractionVSWeight/collapsedFractionVSWeightVector.size() << "\n";
    }
  } printFileCollapsedFractionVSWeight.close(); printFileCollapsedFractionVSWeightFitting.close();
  
  //Collapsed fraction of hyperedge neighbors
  ofstream printFileCollapsedFractionOfHyperedgeNeighborsFirstLine("3. Collapsed fraction of hyperedge neighbors/CollapsedFractionOfHyperedgeNeighbors0000.txt"), printFileCollapsedFractionOfHyperedgeOneNeighborsFirstLine("3. Collapsed fraction of hyperedge neighbors/CollapsedFractionOfHyperedgeOneNeighbors0000.txt"), printFileCollapsedFractionOfHyperedgeTwoNeighborsFirstLine("3. Collapsed fraction of hyperedge neighbors/CollapsedFractionOfHyperedgeTwoNeighbors0000.txt"), printFileCollapsedFractionOfHyperedgeThreeNeighborsFirstLine("3. Collapsed fraction of hyperedge neighbors/CollapsedFractionOfHyperedgeThreeNeighbors0000.txt");
  printFileCollapsedFractionOfHyperedgeNeighborsFirstLine << "Year" << "\t" << "Ensembles" << "\t" << "Samples" << "\t" << "Normal" << "\t" << "<X^2>" << "\t" << "Collapsed" << "\t" << "<X^2>" << "\n";
  printFileCollapsedFractionOfHyperedgeOneNeighborsFirstLine << "Year" << "\t" << "Ensembles" << "\t" << "Samples" << "\t" << "Normal" << "\t" << "<X^2>" << "\t" << "Collapsed" << "\t" << "<X^2>" << "\n";
  printFileCollapsedFractionOfHyperedgeTwoNeighborsFirstLine << "Year" << "\t" << "Ensembles" << "\t" << "Samples" << "\t" << "Normal" << "\t" << "<X^2>" << "\t" << "Collapsed" << "\t" << "<X^2>" << "\n";
  printFileCollapsedFractionOfHyperedgeThreeNeighborsFirstLine << "Year" << "\t" << "Ensembles" << "\t" << "Samples" << "\t" << "Normal" << "\t" << "<X^2>" << "\t" << "Collapsed" << "\t" << "<X^2>" << "\n";
  printFileCollapsedFractionOfHyperedgeNeighborsFirstLine.close(); printFileCollapsedFractionOfHyperedgeOneNeighborsFirstLine.close(); printFileCollapsedFractionOfHyperedgeTwoNeighborsFirstLine.close(); printFileCollapsedFractionOfHyperedgeThreeNeighborsFirstLine.close();
  ofstream printFileCollapsedFractionOfHyperedgeNeighbors("3. Collapsed fraction of hyperedge neighbors/CollapsedFractionOfHyperedgeNeighbors"+to_string(year)+".txt");
  long double sumCollapsedFractionOfNormalHyperedgeNeighbors=0, sumSquaredCollapsedFractionOfNormalHyperedgeNeighbors=0, sumCollapsedFractionOfCollapsedHyperedgeNeighbors=0, sumSquaredCollapsedFractionOfCollapsedHyperedgeNeighbors=0;
  for (auto &collapsedFractionOfNormalHyperedgeNeighbors : collapsedFractionOfNormalHyperedgeNeighborsVector) {
    sumCollapsedFractionOfNormalHyperedgeNeighbors+=collapsedFractionOfNormalHyperedgeNeighbors.first/collapsedFractionOfNormalHyperedgeNeighbors.second;
    sumSquaredCollapsedFractionOfNormalHyperedgeNeighbors+=(collapsedFractionOfNormalHyperedgeNeighbors.first/collapsedFractionOfNormalHyperedgeNeighbors.second)*(collapsedFractionOfNormalHyperedgeNeighbors.first/collapsedFractionOfNormalHyperedgeNeighbors.second);
  }
  for (auto &collapsedFractionOfCollapsedHyperedgeNeighbors : collapsedFractionOfCollapsedHyperedgeNeighborsVector) {
    sumCollapsedFractionOfCollapsedHyperedgeNeighbors+=collapsedFractionOfCollapsedHyperedgeNeighbors.first/collapsedFractionOfCollapsedHyperedgeNeighbors.second;
    sumSquaredCollapsedFractionOfCollapsedHyperedgeNeighbors+=(collapsedFractionOfCollapsedHyperedgeNeighbors.first/collapsedFractionOfCollapsedHyperedgeNeighbors.second)*(collapsedFractionOfCollapsedHyperedgeNeighbors.first/collapsedFractionOfCollapsedHyperedgeNeighbors.second);
  }
  printFileCollapsedFractionOfHyperedgeNeighbors << year << "\t" << ensembleIndex << "\t" << collapsedFractionOfNormalHyperedgeNeighborsVector.size() << "\t" << sumCollapsedFractionOfNormalHyperedgeNeighbors/collapsedFractionOfNormalHyperedgeNeighborsVector.size() << "\t" << sumSquaredCollapsedFractionOfNormalHyperedgeNeighbors/collapsedFractionOfNormalHyperedgeNeighborsVector.size() << "\t" << sumCollapsedFractionOfCollapsedHyperedgeNeighbors/collapsedFractionOfCollapsedHyperedgeNeighborsVector.size() << "\t" << sumSquaredCollapsedFractionOfCollapsedHyperedgeNeighbors/collapsedFractionOfCollapsedHyperedgeNeighborsVector.size() << "\n";
  printFileCollapsedFractionOfHyperedgeNeighbors.close();
  ofstream printFileCollapsedFractionOfHyperedgeOneNeighbors("3. Collapsed fraction of hyperedge neighbors/CollapsedFractionOfHyperedgeOneNeighbors"+to_string(year)+".txt");
  long double sumCollapsedFractionOfNormalHyperedgeOneNeighbors=0, sumSquaredCollapsedFractionOfNormalHyperedgeOneNeighbors=0, sumCollapsedFractionOfCollapsedHyperedgeOneNeighbors=0, sumSquaredCollapsedFractionOfCollapsedHyperedgeOneNeighbors=0;
  for (auto &collapsedFractionOfNormalHyperedgeOneNeighbors : collapsedFractionOfNormalHyperedgeOneNeighborsVector) {
    sumCollapsedFractionOfNormalHyperedgeOneNeighbors+=collapsedFractionOfNormalHyperedgeOneNeighbors.first/collapsedFractionOfNormalHyperedgeOneNeighbors.second;
    sumSquaredCollapsedFractionOfNormalHyperedgeOneNeighbors+=(collapsedFractionOfNormalHyperedgeOneNeighbors.first/collapsedFractionOfNormalHyperedgeOneNeighbors.second)*(collapsedFractionOfNormalHyperedgeOneNeighbors.first/collapsedFractionOfNormalHyperedgeOneNeighbors.second);
  }
  for (auto &collapsedFractionOfCollapsedHyperedgeOneNeighbors : collapsedFractionOfCollapsedHyperedgeOneNeighborsVector) {
    sumCollapsedFractionOfCollapsedHyperedgeOneNeighbors+=collapsedFractionOfCollapsedHyperedgeOneNeighbors.first/collapsedFractionOfCollapsedHyperedgeOneNeighbors.second;
    sumSquaredCollapsedFractionOfCollapsedHyperedgeOneNeighbors+=(collapsedFractionOfCollapsedHyperedgeOneNeighbors.first/collapsedFractionOfCollapsedHyperedgeOneNeighbors.second)*(collapsedFractionOfCollapsedHyperedgeOneNeighbors.first/collapsedFractionOfCollapsedHyperedgeOneNeighbors.second);
  }
  printFileCollapsedFractionOfHyperedgeOneNeighbors << year << "\t" << ensembleIndex << "\t" << collapsedFractionOfNormalHyperedgeOneNeighborsVector.size() << "\t" << sumCollapsedFractionOfNormalHyperedgeOneNeighbors/collapsedFractionOfNormalHyperedgeOneNeighborsVector.size() << "\t" << sumSquaredCollapsedFractionOfNormalHyperedgeOneNeighbors/collapsedFractionOfNormalHyperedgeOneNeighborsVector.size() << "\t" << sumCollapsedFractionOfCollapsedHyperedgeOneNeighbors/collapsedFractionOfCollapsedHyperedgeOneNeighborsVector.size() << "\t" << sumSquaredCollapsedFractionOfCollapsedHyperedgeOneNeighbors/collapsedFractionOfCollapsedHyperedgeOneNeighborsVector.size() << "\n";
  printFileCollapsedFractionOfHyperedgeOneNeighbors.close();
  ofstream printFileCollapsedFractionOfHyperedgeTwoNeighbors("3. Collapsed fraction of hyperedge neighbors/CollapsedFractionOfHyperedgeTwoNeighbors"+to_string(year)+".txt");
  long double sumCollapsedFractionOfNormalHyperedgeTwoNeighbors=0, sumSquaredCollapsedFractionOfNormalHyperedgeTwoNeighbors=0, sumCollapsedFractionOfCollapsedHyperedgeTwoNeighbors=0, sumSquaredCollapsedFractionOfCollapsedHyperedgeTwoNeighbors=0;
  for (auto &collapsedFractionOfNormalHyperedgeTwoNeighbors : collapsedFractionOfNormalHyperedgeTwoNeighborsVector) {
    sumCollapsedFractionOfNormalHyperedgeTwoNeighbors+=collapsedFractionOfNormalHyperedgeTwoNeighbors.first/collapsedFractionOfNormalHyperedgeTwoNeighbors.second;
    sumSquaredCollapsedFractionOfNormalHyperedgeTwoNeighbors+=(collapsedFractionOfNormalHyperedgeTwoNeighbors.first/collapsedFractionOfNormalHyperedgeTwoNeighbors.second)*(collapsedFractionOfNormalHyperedgeTwoNeighbors.first/collapsedFractionOfNormalHyperedgeTwoNeighbors.second);
  }
  for (auto &collapsedFractionOfCollapsedHyperedgeTwoNeighbors : collapsedFractionOfCollapsedHyperedgeTwoNeighborsVector) {
    sumCollapsedFractionOfCollapsedHyperedgeTwoNeighbors+=collapsedFractionOfCollapsedHyperedgeTwoNeighbors.first/collapsedFractionOfCollapsedHyperedgeTwoNeighbors.second;
    sumSquaredCollapsedFractionOfCollapsedHyperedgeTwoNeighbors+=(collapsedFractionOfCollapsedHyperedgeTwoNeighbors.first/collapsedFractionOfCollapsedHyperedgeTwoNeighbors.second)*(collapsedFractionOfCollapsedHyperedgeTwoNeighbors.first/collapsedFractionOfCollapsedHyperedgeTwoNeighbors.second);
  }
  printFileCollapsedFractionOfHyperedgeTwoNeighbors << year << "\t" << ensembleIndex << "\t" << collapsedFractionOfNormalHyperedgeTwoNeighborsVector.size() << "\t" << sumCollapsedFractionOfNormalHyperedgeTwoNeighbors/collapsedFractionOfNormalHyperedgeTwoNeighborsVector.size() << "\t" << sumSquaredCollapsedFractionOfNormalHyperedgeTwoNeighbors/collapsedFractionOfNormalHyperedgeTwoNeighborsVector.size() << "\t" << sumCollapsedFractionOfCollapsedHyperedgeTwoNeighbors/collapsedFractionOfCollapsedHyperedgeTwoNeighborsVector.size() << "\t" << sumSquaredCollapsedFractionOfCollapsedHyperedgeTwoNeighbors/collapsedFractionOfCollapsedHyperedgeTwoNeighborsVector.size() << "\n";
  printFileCollapsedFractionOfHyperedgeTwoNeighbors.close();
  ofstream printFileCollapsedFractionOfHyperedgeThreeNeighbors("3. Collapsed fraction of hyperedge neighbors/CollapsedFractionOfHyperedgeThreeNeighbors"+to_string(year)+".txt");
  long double sumCollapsedFractionOfNormalHyperedgeThreeNeighbors=0, sumSquaredCollapsedFractionOfNormalHyperedgeThreeNeighbors=0, sumCollapsedFractionOfCollapsedHyperedgeThreeNeighbors=0, sumSquaredCollapsedFractionOfCollapsedHyperedgeThreeNeighbors=0;
  for (auto &collapsedFractionOfNormalHyperedgeThreeNeighbors : collapsedFractionOfNormalHyperedgeThreeNeighborsVector) {
    sumCollapsedFractionOfNormalHyperedgeThreeNeighbors+=collapsedFractionOfNormalHyperedgeThreeNeighbors.first/collapsedFractionOfNormalHyperedgeThreeNeighbors.second;
    sumSquaredCollapsedFractionOfNormalHyperedgeThreeNeighbors+=(collapsedFractionOfNormalHyperedgeThreeNeighbors.first/collapsedFractionOfNormalHyperedgeThreeNeighbors.second)*(collapsedFractionOfNormalHyperedgeThreeNeighbors.first/collapsedFractionOfNormalHyperedgeThreeNeighbors.second);
  }
  for (auto &collapsedFractionOfCollapsedHyperedgeThreeNeighbors : collapsedFractionOfCollapsedHyperedgeThreeNeighborsVector) {
    sumCollapsedFractionOfCollapsedHyperedgeThreeNeighbors+=collapsedFractionOfCollapsedHyperedgeThreeNeighbors.first/collapsedFractionOfCollapsedHyperedgeThreeNeighbors.second;
    sumSquaredCollapsedFractionOfCollapsedHyperedgeThreeNeighbors+=(collapsedFractionOfCollapsedHyperedgeThreeNeighbors.first/collapsedFractionOfCollapsedHyperedgeThreeNeighbors.second)*(collapsedFractionOfCollapsedHyperedgeThreeNeighbors.first/collapsedFractionOfCollapsedHyperedgeThreeNeighbors.second);
  }
  printFileCollapsedFractionOfHyperedgeThreeNeighbors << year << "\t" << ensembleIndex << "\t" << collapsedFractionOfNormalHyperedgeThreeNeighborsVector.size() << "\t" << sumCollapsedFractionOfNormalHyperedgeThreeNeighbors/collapsedFractionOfNormalHyperedgeThreeNeighborsVector.size() << "\t" << sumSquaredCollapsedFractionOfNormalHyperedgeThreeNeighbors/collapsedFractionOfNormalHyperedgeThreeNeighborsVector.size() << "\t" << sumCollapsedFractionOfCollapsedHyperedgeThreeNeighbors/collapsedFractionOfCollapsedHyperedgeThreeNeighborsVector.size() << "\t" << sumSquaredCollapsedFractionOfCollapsedHyperedgeThreeNeighbors/collapsedFractionOfCollapsedHyperedgeThreeNeighborsVector.size() << "\n";
  printFileCollapsedFractionOfHyperedgeThreeNeighbors.close();
  
  //Collapsed fraction of hyperedge neighbors VS weight
  ofstream printFileCollapsedFractionOfHyperedgeNeighborsVSWeight("4. Collapsed fraction of hyperedge neighbors VS weight/CollapsedFractionOfHyperedgeNeighborsVSWeight"+to_string(year)+".txt");
  printFileCollapsedFractionOfHyperedgeNeighborsVSWeight << "Ensembles" << "\t" << "w/<w>" << "\t" << "Samples" << "\t" << "Normal" << "\t" << "<X^2>" << "\t" << "Samples" << "\t" << "Collapsed" << "\t" << "<X^2>" << "\n";
  for (unsigned int binIndex=0; binIndex<collapsedFractionOfNormalHyperedgeNeighborsVSWeightVector.front().size(); ++binIndex) {
    long double sumCollapsedFractionOfNormalHyperedgeNeighborsVSWeight=0, sumSquaredCollapsedFractionOfNormalHyperedgeNeighborsVSWeight=0, sumCollapsedFractionOfCollapsedHyperedgeNeighborsVSWeight=0, sumSquaredCollapsedFractionOfCollapsedHyperedgeNeighborsVSWeight=0;
    unsigned int normalSampleNumber=0, collapsedSampleNumber=0;
    for (auto &collapsedFractionOfNormalHyperedgeNeighborsVSWeight : collapsedFractionOfNormalHyperedgeNeighborsVSWeightVector) {
      if (collapsedFractionOfNormalHyperedgeNeighborsVSWeight[binIndex].second) {
        sumCollapsedFractionOfNormalHyperedgeNeighborsVSWeight+=collapsedFractionOfNormalHyperedgeNeighborsVSWeight[binIndex].first/collapsedFractionOfNormalHyperedgeNeighborsVSWeight[binIndex].second;
        sumSquaredCollapsedFractionOfNormalHyperedgeNeighborsVSWeight+=(collapsedFractionOfNormalHyperedgeNeighborsVSWeight[binIndex].first/collapsedFractionOfNormalHyperedgeNeighborsVSWeight[binIndex].second)*(collapsedFractionOfNormalHyperedgeNeighborsVSWeight[binIndex].first/collapsedFractionOfNormalHyperedgeNeighborsVSWeight[binIndex].second);
        ++normalSampleNumber;
      }
    }
    for (auto &collapsedFractionOfCollapsedHyperedgeNeighborsVSWeight : collapsedFractionOfCollapsedHyperedgeNeighborsVSWeightVector) {
      if (collapsedFractionOfCollapsedHyperedgeNeighborsVSWeight[binIndex].second) {
        sumCollapsedFractionOfCollapsedHyperedgeNeighborsVSWeight+=collapsedFractionOfCollapsedHyperedgeNeighborsVSWeight[binIndex].first/collapsedFractionOfCollapsedHyperedgeNeighborsVSWeight[binIndex].second;
        sumSquaredCollapsedFractionOfCollapsedHyperedgeNeighborsVSWeight+=(collapsedFractionOfCollapsedHyperedgeNeighborsVSWeight[binIndex].first/collapsedFractionOfCollapsedHyperedgeNeighborsVSWeight[binIndex].second)*(collapsedFractionOfCollapsedHyperedgeNeighborsVSWeight[binIndex].first/collapsedFractionOfCollapsedHyperedgeNeighborsVSWeight[binIndex].second);
        ++collapsedSampleNumber;
      }
    }
    printFileCollapsedFractionOfHyperedgeNeighborsVSWeight << ensembleIndex << "\t" << pow(binSize, binIndex+0.5)*scale << "\t" << normalSampleNumber << "\t" << sumCollapsedFractionOfNormalHyperedgeNeighborsVSWeight/normalSampleNumber << "\t" << sumSquaredCollapsedFractionOfNormalHyperedgeNeighborsVSWeight/normalSampleNumber << "\t" << collapsedSampleNumber << "\t" << sumCollapsedFractionOfCollapsedHyperedgeNeighborsVSWeight/collapsedSampleNumber << "\t" << sumSquaredCollapsedFractionOfCollapsedHyperedgeNeighborsVSWeight/collapsedSampleNumber << "\n";
  } printFileCollapsedFractionOfHyperedgeNeighborsVSWeight.close();
  ofstream printFileCollapsedFractionOfHyperedgeOneNeighborsVSWeight("4. Collapsed fraction of hyperedge neighbors VS weight/CollapsedFractionOfHyperedgeOneNeighborsVSWeight"+to_string(year)+".txt");
  printFileCollapsedFractionOfHyperedgeOneNeighborsVSWeight << "Ensembles" << "\t" << "w/<w>" << "\t" << "Samples" << "\t" << "Normal" << "\t" << "<X^2>" << "\t" << "Samples" << "\t" << "Collapsed" << "\t" << "<X^2>" << "\n";
  for (unsigned int binIndex=0; binIndex<collapsedFractionOfNormalHyperedgeOneNeighborsVSWeightVector.front().size(); ++binIndex) {
    long double sumCollapsedFractionOfNormalHyperedgeOneNeighborsVSWeight=0, sumSquaredCollapsedFractionOfNormalHyperedgeOneNeighborsVSWeight=0, sumCollapsedFractionOfCollapsedHyperedgeOneNeighborsVSWeight=0, sumSquaredCollapsedFractionOfCollapsedHyperedgeOneNeighborsVSWeight=0;
    unsigned int normalSampleNumber=0, collapsedSampleNumber=0;
    for (auto &collapsedFractionOfNormalHyperedgeOneNeighborsVSWeight : collapsedFractionOfNormalHyperedgeOneNeighborsVSWeightVector) {
      if (collapsedFractionOfNormalHyperedgeOneNeighborsVSWeight[binIndex].second) {
        sumCollapsedFractionOfNormalHyperedgeOneNeighborsVSWeight+=collapsedFractionOfNormalHyperedgeOneNeighborsVSWeight[binIndex].first/collapsedFractionOfNormalHyperedgeOneNeighborsVSWeight[binIndex].second;
        sumSquaredCollapsedFractionOfNormalHyperedgeOneNeighborsVSWeight+=(collapsedFractionOfNormalHyperedgeOneNeighborsVSWeight[binIndex].first/collapsedFractionOfNormalHyperedgeOneNeighborsVSWeight[binIndex].second)*(collapsedFractionOfNormalHyperedgeOneNeighborsVSWeight[binIndex].first/collapsedFractionOfNormalHyperedgeOneNeighborsVSWeight[binIndex].second);
        ++normalSampleNumber;
      }
    }
    for (auto &collapsedFractionOfCollapsedHyperedgeOneNeighborsVSWeight : collapsedFractionOfCollapsedHyperedgeOneNeighborsVSWeightVector) {
      if (collapsedFractionOfCollapsedHyperedgeOneNeighborsVSWeight[binIndex].second) {
        sumCollapsedFractionOfCollapsedHyperedgeOneNeighborsVSWeight+=collapsedFractionOfCollapsedHyperedgeOneNeighborsVSWeight[binIndex].first/collapsedFractionOfCollapsedHyperedgeOneNeighborsVSWeight[binIndex].second;
        sumSquaredCollapsedFractionOfCollapsedHyperedgeOneNeighborsVSWeight+=(collapsedFractionOfCollapsedHyperedgeOneNeighborsVSWeight[binIndex].first/collapsedFractionOfCollapsedHyperedgeOneNeighborsVSWeight[binIndex].second)*(collapsedFractionOfCollapsedHyperedgeOneNeighborsVSWeight[binIndex].first/collapsedFractionOfCollapsedHyperedgeOneNeighborsVSWeight[binIndex].second);
        ++collapsedSampleNumber;
      }
    }
    printFileCollapsedFractionOfHyperedgeOneNeighborsVSWeight << ensembleIndex << "\t" << pow(binSize, binIndex+0.5)*scale << "\t" << normalSampleNumber << "\t" << sumCollapsedFractionOfNormalHyperedgeOneNeighborsVSWeight/normalSampleNumber << "\t" << sumSquaredCollapsedFractionOfNormalHyperedgeOneNeighborsVSWeight/normalSampleNumber << "\t" << collapsedSampleNumber << "\t" << sumCollapsedFractionOfCollapsedHyperedgeOneNeighborsVSWeight/collapsedSampleNumber << "\t" << sumSquaredCollapsedFractionOfCollapsedHyperedgeOneNeighborsVSWeight/collapsedSampleNumber << "\n";
  } printFileCollapsedFractionOfHyperedgeOneNeighborsVSWeight.close();
  ofstream printFileCollapsedFractionOfHyperedgeTwoNeighborsVSWeight("4. Collapsed fraction of hyperedge neighbors VS weight/CollapsedFractionOfHyperedgeTwoNeighborsVSWeight"+to_string(year)+".txt");
  printFileCollapsedFractionOfHyperedgeTwoNeighborsVSWeight << "Ensembles" << "\t" << "w/<w>" << "\t" << "Samples" << "\t" << "Normal" << "\t" << "<X^2>" << "\t" << "Samples" << "\t" << "Collapsed" << "\t" << "<X^2>" << "\n";
  for (unsigned int binIndex=0; binIndex<collapsedFractionOfNormalHyperedgeTwoNeighborsVSWeightVector.front().size(); ++binIndex) {
    long double sumCollapsedFractionOfNormalHyperedgeTwoNeighborsVSWeight=0, sumSquaredCollapsedFractionOfNormalHyperedgeTwoNeighborsVSWeight=0, sumCollapsedFractionOfCollapsedHyperedgeTwoNeighborsVSWeight=0, sumSquaredCollapsedFractionOfCollapsedHyperedgeTwoNeighborsVSWeight=0;
    unsigned int normalSampleNumber=0, collapsedSampleNumber=0;
    for (auto &collapsedFractionOfNormalHyperedgeTwoNeighborsVSWeight : collapsedFractionOfNormalHyperedgeTwoNeighborsVSWeightVector) {
      if (collapsedFractionOfNormalHyperedgeTwoNeighborsVSWeight[binIndex].second) {
        sumCollapsedFractionOfNormalHyperedgeTwoNeighborsVSWeight+=collapsedFractionOfNormalHyperedgeTwoNeighborsVSWeight[binIndex].first/collapsedFractionOfNormalHyperedgeTwoNeighborsVSWeight[binIndex].second;
        sumSquaredCollapsedFractionOfNormalHyperedgeTwoNeighborsVSWeight+=(collapsedFractionOfNormalHyperedgeTwoNeighborsVSWeight[binIndex].first/collapsedFractionOfNormalHyperedgeTwoNeighborsVSWeight[binIndex].second)*(collapsedFractionOfNormalHyperedgeTwoNeighborsVSWeight[binIndex].first/collapsedFractionOfNormalHyperedgeTwoNeighborsVSWeight[binIndex].second);
        ++normalSampleNumber;
      }
    }
    for (auto &collapsedFractionOfCollapsedHyperedgeTwoNeighborsVSWeight : collapsedFractionOfCollapsedHyperedgeTwoNeighborsVSWeightVector) {
      if (collapsedFractionOfCollapsedHyperedgeTwoNeighborsVSWeight[binIndex].second) {
        sumCollapsedFractionOfCollapsedHyperedgeTwoNeighborsVSWeight+=collapsedFractionOfCollapsedHyperedgeTwoNeighborsVSWeight[binIndex].first/collapsedFractionOfCollapsedHyperedgeTwoNeighborsVSWeight[binIndex].second;
        sumSquaredCollapsedFractionOfCollapsedHyperedgeTwoNeighborsVSWeight+=(collapsedFractionOfCollapsedHyperedgeTwoNeighborsVSWeight[binIndex].first/collapsedFractionOfCollapsedHyperedgeTwoNeighborsVSWeight[binIndex].second)*(collapsedFractionOfCollapsedHyperedgeTwoNeighborsVSWeight[binIndex].first/collapsedFractionOfCollapsedHyperedgeTwoNeighborsVSWeight[binIndex].second);
        ++collapsedSampleNumber;
      }
    }
    printFileCollapsedFractionOfHyperedgeTwoNeighborsVSWeight << ensembleIndex << "\t" << pow(binSize, binIndex+0.5)*scale << "\t" << normalSampleNumber << "\t" << sumCollapsedFractionOfNormalHyperedgeTwoNeighborsVSWeight/normalSampleNumber << "\t" << sumSquaredCollapsedFractionOfNormalHyperedgeTwoNeighborsVSWeight/normalSampleNumber << "\t" << collapsedSampleNumber << "\t" << sumCollapsedFractionOfCollapsedHyperedgeTwoNeighborsVSWeight/collapsedSampleNumber << "\t" << sumSquaredCollapsedFractionOfCollapsedHyperedgeTwoNeighborsVSWeight/collapsedSampleNumber << "\n";
  } printFileCollapsedFractionOfHyperedgeTwoNeighborsVSWeight.close();
  ofstream printFileCollapsedFractionOfHyperedgeThreeNeighborsVSWeight("4. Collapsed fraction of hyperedge neighbors VS weight/CollapsedFractionOfHyperedgeThreeNeighborsVSWeight"+to_string(year)+".txt");
  printFileCollapsedFractionOfHyperedgeThreeNeighborsVSWeight << "Ensembles" << "\t" << "w/<w>" << "\t" << "Samples" << "\t" << "Normal" << "\t" << "<X^2>" << "\t" << "Samples" << "\t" << "Collapsed" << "\t" << "<X^2>" << "\n";
  for (unsigned int binIndex=0; binIndex<collapsedFractionOfNormalHyperedgeThreeNeighborsVSWeightVector.front().size(); ++binIndex) {
    long double sumCollapsedFractionOfNormalHyperedgeThreeNeighborsVSWeight=0, sumSquaredCollapsedFractionOfNormalHyperedgeThreeNeighborsVSWeight=0, sumCollapsedFractionOfCollapsedHyperedgeThreeNeighborsVSWeight=0, sumSquaredCollapsedFractionOfCollapsedHyperedgeThreeNeighborsVSWeight=0;
    unsigned int normalSampleNumber=0, collapsedSampleNumber=0;
    for (auto &collapsedFractionOfNormalHyperedgeThreeNeighborsVSWeight : collapsedFractionOfNormalHyperedgeThreeNeighborsVSWeightVector) {
      if (collapsedFractionOfNormalHyperedgeThreeNeighborsVSWeight[binIndex].second) {
        sumCollapsedFractionOfNormalHyperedgeThreeNeighborsVSWeight+=collapsedFractionOfNormalHyperedgeThreeNeighborsVSWeight[binIndex].first/collapsedFractionOfNormalHyperedgeThreeNeighborsVSWeight[binIndex].second;
        sumSquaredCollapsedFractionOfNormalHyperedgeThreeNeighborsVSWeight+=(collapsedFractionOfNormalHyperedgeThreeNeighborsVSWeight[binIndex].first/collapsedFractionOfNormalHyperedgeThreeNeighborsVSWeight[binIndex].second)*(collapsedFractionOfNormalHyperedgeThreeNeighborsVSWeight[binIndex].first/collapsedFractionOfNormalHyperedgeThreeNeighborsVSWeight[binIndex].second);
        ++normalSampleNumber;
      }
    }
    for (auto &collapsedFractionOfCollapsedHyperedgeThreeNeighborsVSWeight : collapsedFractionOfCollapsedHyperedgeThreeNeighborsVSWeightVector) {
      if (collapsedFractionOfCollapsedHyperedgeThreeNeighborsVSWeight[binIndex].second) {
        sumCollapsedFractionOfCollapsedHyperedgeThreeNeighborsVSWeight+=collapsedFractionOfCollapsedHyperedgeThreeNeighborsVSWeight[binIndex].first/collapsedFractionOfCollapsedHyperedgeThreeNeighborsVSWeight[binIndex].second;
        sumSquaredCollapsedFractionOfCollapsedHyperedgeThreeNeighborsVSWeight+=(collapsedFractionOfCollapsedHyperedgeThreeNeighborsVSWeight[binIndex].first/collapsedFractionOfCollapsedHyperedgeThreeNeighborsVSWeight[binIndex].second)*(collapsedFractionOfCollapsedHyperedgeThreeNeighborsVSWeight[binIndex].first/collapsedFractionOfCollapsedHyperedgeThreeNeighborsVSWeight[binIndex].second);
        ++collapsedSampleNumber;
      }
    }
    printFileCollapsedFractionOfHyperedgeThreeNeighborsVSWeight << ensembleIndex << "\t" << pow(binSize, binIndex+0.5)*scale << "\t" << normalSampleNumber << "\t" << sumCollapsedFractionOfNormalHyperedgeThreeNeighborsVSWeight/normalSampleNumber << "\t" << sumSquaredCollapsedFractionOfNormalHyperedgeThreeNeighborsVSWeight/normalSampleNumber << "\t" << collapsedSampleNumber << "\t" << sumCollapsedFractionOfCollapsedHyperedgeThreeNeighborsVSWeight/collapsedSampleNumber << "\t" << sumSquaredCollapsedFractionOfCollapsedHyperedgeThreeNeighborsVSWeight/collapsedSampleNumber << "\n";
  } printFileCollapsedFractionOfHyperedgeThreeNeighborsVSWeight.close();
}
