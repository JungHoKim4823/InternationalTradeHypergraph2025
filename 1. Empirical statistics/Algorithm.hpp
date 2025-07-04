#include<cmath>
#include<utility>
#include<tuple>
#include<vector>
#include<set>
#include<map>
#include<iostream>
#include<fstream>
using namespace std;

extern long double logRatioThreshold;
extern string stringLogRatioThreshold;

class Country; class Product; class Hyperedge;
class Country {public: vector<Hyperedge*> exporterHyperedgeVector, importerHyperedgeVector, countryHyperedgeVector;};
class Product {public: vector<Hyperedge*> productHyperedgeVector;};
class Hyperedge {public: Country *exporter, *importer; Product *product; long double weight, logRatio;};

class Hypergraph {
public:
  map<unsigned int, Country> countryMap; map<unsigned int, Product> productMap; vector<Hyperedge> hyperedgeVector;
  Hypergraph(unsigned int year);
};

Hypergraph::Hypergraph(unsigned int year) {
  cout << year << "\n";
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
    if (logRatio=="-inf") {hyperedge.logRatio=-INFINITY;} else {hyperedge.logRatio=stold(logRatio);}
    hyperedgeVector.emplace_back(hyperedge);
    countryMap[exporter].exporterHyperedgeVector.emplace_back(&hyperedgeVector.back());
    countryMap[exporter].countryHyperedgeVector.emplace_back(&hyperedgeVector.back());
    countryMap[importer].importerHyperedgeVector.emplace_back(&hyperedgeVector.back());
    countryMap[importer].countryHyperedgeVector.emplace_back(&hyperedgeVector.back());
    productMap[product].productHyperedgeVector.emplace_back(&hyperedgeVector.back());
  } scanFile.close();

  //The number of nodes and hyperedges
  ofstream printFileTheNumberOfNodesAndHyperedges("01. The number of nodes and hyperedges/TheNumberOfNodesAndHyperedges.txt", ios::app);
  printFileTheNumberOfNodesAndHyperedges << year << "\t" << countryMap.size() << "\t" << productMap.size() << "\t" << hyperedgeVector.size() << "\n";
  printFileTheNumberOfNodesAndHyperedges.close();
  
  //Weight distribution
  long double averageWeight=0, averageSquaredWeight=0, maximumWeight=0;
  for (auto &hyperedge : hyperedgeVector) {
    averageWeight+=hyperedge.weight; averageSquaredWeight+=hyperedge.weight*hyperedge.weight; if (maximumWeight<hyperedge.weight) {maximumWeight=hyperedge.weight;}
  } averageWeight/=hyperedgeVector.size(); averageSquaredWeight/=hyperedgeVector.size();
  long double binSize=4, scale=1e-8; unsigned int pointNumber=log(maximumWeight/(scale*averageWeight))/log(binSize)+1;
  vector<unsigned int> numberVSWeight(pointNumber);
  for (auto &hyperedge : hyperedgeVector) {
    unsigned int binIndex=log(hyperedge.weight/(scale*averageWeight))/log(binSize);
    ++numberVSWeight[binIndex];
  }
  ofstream printFileWeightDistribution("02. Weight distibutions/WeightDistribution"+to_string(year)+".txt");
  printFileWeightDistribution << "w/<w>" << "\t" << "P(w/<w>)" << "\t" << "Samples" << "\n";
  for (unsigned int binIndex=0; binIndex<pointNumber; ++binIndex) {
    printFileWeightDistribution << pow(binSize, binIndex+0.5)*scale << "\t" << numberVSWeight[binIndex]/(hyperedgeVector.size()*(pow(binSize, binIndex+1)-pow(binSize, binIndex))*scale) << "\t" << numberVSWeight[binIndex] << "\n";
  } printFileWeightDistribution.close();
  ofstream printFileAverageWeight("02. Weight distibutions/AverageWeight.txt", ios::app);
  printFileAverageWeight << year << "\t" << averageWeight << "\t" << averageSquaredWeight << "\n";
  printFileAverageWeight.close();
  
  //Degree distribution
  long double averageDegree=0, averageSquaredDegree=0, maximumDegree=0;
  for (auto &hyperedge : hyperedgeVector) {
    unsigned int degree=hyperedge.exporter->countryHyperedgeVector.size()+hyperedge.importer->countryHyperedgeVector.size()+hyperedge.product->productHyperedgeVector.size()-3;
    averageDegree+=degree; averageSquaredDegree+=(long double)degree*degree; if (maximumDegree<degree) {maximumDegree=degree;}
  } averageDegree/=hyperedgeVector.size(); averageSquaredDegree/=hyperedgeVector.size();
  binSize=0.05; pointNumber=(maximumDegree/averageDegree)/binSize+1;
  vector<unsigned int> numberVSDegree(pointNumber);
  for (auto &hyperedge : hyperedgeVector) {
    unsigned int degree=hyperedge.exporter->countryHyperedgeVector.size()+hyperedge.importer->countryHyperedgeVector.size()+hyperedge.product->productHyperedgeVector.size()-3;
    unsigned int binIndex=(degree/averageDegree)/binSize;
    ++numberVSDegree[binIndex];
  }
  ofstream printFileDegreeDistribution("03. Degree distributions/DegreeDistribution"+to_string(year)+".txt");
  printFileDegreeDistribution << "k/<k>" << "\t" << "P(k/<k>)" << "\t" << "Samples" << "\n";
  for (unsigned int binIndex=0; binIndex<pointNumber; ++binIndex) {
    printFileDegreeDistribution << binSize*binIndex+binSize/2 << "\t" << numberVSDegree[binIndex]/(hyperedgeVector.size()*binSize) << "\t" << numberVSDegree[binIndex] << "\n";
  } printFileDegreeDistribution.close();
  ofstream printFileAverageDegree("03. Degree distributions/AverageDegree.txt", ios::app);
  printFileAverageDegree << year << "\t" << averageDegree << "\t" << averageSquaredDegree << "\n";
  printFileAverageDegree.close();
  
  //Finite log ratio distribution
  long double averageFiniteLogRatio=0, averageSquaredFiniteLogRatio=0; unsigned int finiteLogRatioHyperedgeNumber=0;
  for (auto &hyperedge : hyperedgeVector) {
    if (isinf(hyperedge.logRatio)) {continue;}
    averageFiniteLogRatio+=hyperedge.logRatio; averageSquaredFiniteLogRatio+=hyperedge.logRatio*hyperedge.logRatio; ++finiteLogRatioHyperedgeNumber;
  } averageFiniteLogRatio/=finiteLogRatioHyperedgeNumber; averageSquaredFiniteLogRatio/=finiteLogRatioHyperedgeNumber;
  binSize=1; pointNumber=(20*2)/binSize;
  vector<unsigned int> numberVSFiniteLogRatio(pointNumber);
  for (auto &hyperedge : hyperedgeVector) {
    if (isinf(hyperedge.logRatio)) {continue;}
    unsigned int binIndex=((hyperedge.logRatio-averageFiniteLogRatio)+20)/binSize;
    ++numberVSFiniteLogRatio[binIndex];
  }
  ofstream printFileFiniteLogRatioDistribution("04. Finite log ratio distributions/FiniteLogRatioDistribution"+to_string(year)+".txt");
  printFileFiniteLogRatioDistribution << "g-<g>_f" << "\t" << "P(g-<g>_f)" << "\t" << "Samples" << "\n";
  for (unsigned int binIndex=0; binIndex<pointNumber; ++binIndex) {
    printFileFiniteLogRatioDistribution << binSize*binIndex-20+binSize/2 << "\t" << numberVSFiniteLogRatio[binIndex]/(finiteLogRatioHyperedgeNumber*binSize) << "\t" << numberVSFiniteLogRatio[binIndex] << "\n";
  } printFileFiniteLogRatioDistribution.close();
  ofstream printFileAverageFiniteLogRatio("04. Finite log ratio distributions/AverageFiniteLogRatio.txt", ios::app);
  printFileAverageFiniteLogRatio << year << "\t" << averageFiniteLogRatio << "\t" << averageSquaredFiniteLogRatio << "\n";
  printFileAverageFiniteLogRatio.close();
  
  //Weight-degree correlation
  binSize=4, scale=1e-8; pointNumber=log(maximumWeight/(scale*averageWeight))/log(binSize)+2;
  vector<tuple<long double, long double, unsigned int>> weightDegreeCorrelation(pointNumber);
  for (auto &hyperedge : hyperedgeVector) {
    unsigned int degree=hyperedge.exporter->countryHyperedgeVector.size()+hyperedge.importer->countryHyperedgeVector.size()+hyperedge.product->productHyperedgeVector.size()-3;
    unsigned int binIndex=log(hyperedge.weight/(scale*averageWeight))/log(binSize);
    get<0>(weightDegreeCorrelation[binIndex])+=degree/averageDegree;
    get<1>(weightDegreeCorrelation[binIndex])+=(degree/averageDegree)*(degree/averageDegree);
    ++get<2>(weightDegreeCorrelation[binIndex]);
  }
  ofstream printFileWeightDegreeCorrelation("05. Weight-degree correlation/WeightDegreeCorrelation"+to_string(year)+".txt");
  printFileWeightDegreeCorrelation << "w/<w>" << "\t" << "k/<k>" << "\t" << "<X^2>" << "\n";
  for (unsigned int binIndex=0; binIndex<pointNumber; ++binIndex) {
    printFileWeightDegreeCorrelation << pow(binSize, binIndex+0.5)*scale << "\t" << get<0>(weightDegreeCorrelation[binIndex])/get<2>(weightDegreeCorrelation[binIndex]) << "\t" << get<1>(weightDegreeCorrelation[binIndex])/get<2>(weightDegreeCorrelation[binIndex]) << "\n";
  } printFileWeightDegreeCorrelation.close();
  
  //Degree-weight correlation
  binSize=0.05; pointNumber=(maximumDegree/averageDegree)/binSize+1;
  vector<tuple<long double, long double, unsigned int>> degreeWeightCorrelation(pointNumber);
  for (auto &hyperedge : hyperedgeVector) {
    unsigned int degree=hyperedge.exporter->countryHyperedgeVector.size()+hyperedge.importer->countryHyperedgeVector.size()+hyperedge.product->productHyperedgeVector.size()-3;
    unsigned int binIndex=(degree/averageDegree)/binSize;
    get<0>(degreeWeightCorrelation[binIndex])+=hyperedge.weight/averageWeight;
    get<1>(degreeWeightCorrelation[binIndex])+=(hyperedge.weight/averageWeight)*(hyperedge.weight/averageWeight);
    ++get<2>(degreeWeightCorrelation[binIndex]);
  }
  ofstream printFileDegreeWeightCorrelation("05. Weight-degree correlation/DegreeWeightCorrelation"+to_string(year)+".txt");
  printFileDegreeWeightCorrelation << "k/<k>" << "\t" << "w/<w>" << "\t" << "<X^2>" << "\n";
  for (unsigned int binIndex=0; binIndex<pointNumber; ++binIndex) {
    printFileDegreeWeightCorrelation << binSize*binIndex+binSize/2 << "\t" << get<0>(degreeWeightCorrelation[binIndex])/get<2>(degreeWeightCorrelation[binIndex]) << "\t" << get<1>(degreeWeightCorrelation[binIndex])/get<2>(degreeWeightCorrelation[binIndex]) << "\n";
  } printFileDegreeWeightCorrelation.close();
  
  //The number of collapsed hyperedges
  unsigned int infiniteCollapsedHyperedgeNumber=0, finiteCollapsedHyperedgeNumber=0;
  for (auto &hyperedge : hyperedgeVector) {
    if (isinf(hyperedge.logRatio)) {++infiniteCollapsedHyperedgeNumber;}
    else if (hyperedge.logRatio<=logRatioThreshold) {++finiteCollapsedHyperedgeNumber;}
  }
  ofstream printFileTheNumberOfCollapsedHyperedges("07. The number of collapsed hyperedges/"+stringLogRatioThreshold+"/TheNumberOfCollapsedHyperedges.txt", ios::app);
  printFileTheNumberOfCollapsedHyperedges << year << "\t" << hyperedgeVector.size() << "\t" << infiniteCollapsedHyperedgeNumber << "\t" << finiteCollapsedHyperedgeNumber << "\n";
  printFileTheNumberOfCollapsedHyperedges.close();
  
  //Collapsed fraction VS Weight
  binSize=4; scale=1e-8; pointNumber=log(maximumWeight/(scale*averageWeight))/log(binSize)+1;
  vector<pair<unsigned int, unsigned int>> collapsedFractionVSWeight(pointNumber);
  for (auto &hyperedge : hyperedgeVector) {
    unsigned int binIndex=log(hyperedge.weight/(scale*averageWeight))/log(binSize);
    if (hyperedge.logRatio<=logRatioThreshold) {++collapsedFractionVSWeight[binIndex].first;}
    ++collapsedFractionVSWeight[binIndex].second;
  }
  ofstream printFileCollapsedFractionVSWeight("08. Collapsed fraction VS weight/"+stringLogRatioThreshold+"/CollapsedFractionVSWeight"+to_string(year)+".txt");
  ofstream printFileCollapsedFractionVSWeightFitting("08. Collapsed fraction VS weight/"+stringLogRatioThreshold+"/CollapsedFractionVSWeightFitting"+to_string(year)+".txt");
  printFileCollapsedFractionVSWeight << "w/<w>" << "\t" << "Collapsed fraction" << "\n";
  printFileCollapsedFractionVSWeightFitting << "w/<w>" << "\t" << "Collapsed fraction" << "\n";
  for (unsigned int binIndex=0; binIndex<pointNumber; ++binIndex) {
    printFileCollapsedFractionVSWeight << pow(binSize, binIndex+0.5)*scale << "\t" << (long double)collapsedFractionVSWeight[binIndex].first/collapsedFractionVSWeight[binIndex].second << "\n";
    if (binIndex>=9 && binIndex<=16) {
      printFileCollapsedFractionVSWeightFitting << pow(binSize, binIndex+0.5)*scale << "\t" << (long double)collapsedFractionVSWeight[binIndex].first/collapsedFractionVSWeight[binIndex].second << "\n";
    }
  } printFileCollapsedFractionVSWeight.close(); printFileCollapsedFractionVSWeightFitting.close();
  
  //Baseline of collapsed fraction VS Weight
  for (unsigned int exponent=1; exponent<=5; exponent+=2) {
    long double blockWeight=pow(10, exponent);
    long double averageDeltaWeight=0, averageSquaredDeltaWeight=0, varianceDeltaWeight; unsigned int blockHyperedgeNumber=0;
    for (auto &hyperedge : hyperedgeVector) {
      if (hyperedge.weight>=0.95*blockWeight && hyperedge.weight<=1.05*blockWeight) {
        long double deltaWeight;
        if (isinf(hyperedge.logRatio)) {deltaWeight=-hyperedge.weight;}
        else {deltaWeight=hyperedge.weight*(exp(hyperedge.logRatio)-1);}
        averageDeltaWeight+=deltaWeight; averageSquaredDeltaWeight+=deltaWeight*deltaWeight; ++blockHyperedgeNumber;
      }
    } averageDeltaWeight/=blockHyperedgeNumber; averageSquaredDeltaWeight/=blockHyperedgeNumber;
    varianceDeltaWeight=averageSquaredDeltaWeight-averageDeltaWeight*averageDeltaWeight;
    binSize=4; scale=1e-8; pointNumber=log(maximumWeight/(scale*averageWeight))/log(binSize)+1;
    ofstream printFileCollapsedFractionVSWeightBaseline("08. Collapsed fraction VS weight/"+stringLogRatioThreshold+"/CollapsedFractionVSWeight"+to_string(exponent)+"Baseline"+to_string(year)+".txt");
    printFileCollapsedFractionVSWeightBaseline << "w/<w>" << "\t" << "Collapsed fraction" << "\n";
    for (unsigned int binIndex=0; binIndex<pointNumber; ++binIndex) {
      long double pointWeight=pow(binSize, binIndex+0.5)*scale*averageWeight;
      long double blockNumber=(unsigned int)(pointWeight/blockWeight)+1;
      long double average=averageDeltaWeight*blockNumber, variance=varianceDeltaWeight*blockNumber;
      long double threshold=pointWeight*(exp(logRatioThreshold)-1);
      long double collapsedFraction=0.5*(erfc(-(threshold-average)/sqrt(2*variance)));
      printFileCollapsedFractionVSWeightBaseline << pow(binSize, binIndex+0.5)*scale << "\t" << collapsedFraction << "\n";
    } printFileCollapsedFractionVSWeightBaseline.close();
  }
  
  //Neighbor related quantities
  binSize=4; scale=1e-8; pointNumber=log(maximumWeight/(scale*averageWeight))/log(binSize)+1;
  vector<tuple<long double, long double, unsigned int>> weightNeighborsWeightCorrelation(pointNumber), weightOneNeighborsWeightCorrelation(pointNumber), weightTwoNeighborsWeightCorrelation(pointNumber), weightThreeNeighborsWeightCorrelation(pointNumber);
  tuple<long double, long double, unsigned int> collapsedFractionOfNormalHyperedgeNeighbors, collapsedFractionOfNormalHyperedgeOneNeighbors, collapsedFractionOfNormalHyperedgeTwoNeighbors, collapsedFractionOfNormalHyperedgeThreeNeighbors;
  tuple<long double, long double, unsigned int> collapsedFractionOfCollapsedHyperedgeNeighbors, collapsedFractionOfCollapsedHyperedgeOneNeighbors, collapsedFractionOfCollapsedHyperedgeTwoNeighbors, collapsedFractionOfCollapsedHyperedgeThreeNeighbors;
  vector<tuple<long double, long double, unsigned int>> collapsedFractionOfNormalHyperedgeNeighborsVSWeight(pointNumber), collapsedFractionOfNormalHyperedgeOneNeighborsVSWeight(pointNumber), collapsedFractionOfNormalHyperedgeTwoNeighborsVSWeight(pointNumber), collapsedFractionOfNormalHyperedgeThreeNeighborsVSWeight(pointNumber);
  vector<tuple<long double, long double, unsigned int>> collapsedFractionOfCollapsedHyperedgeNeighborsVSWeight(pointNumber), collapsedFractionOfCollapsedHyperedgeOneNeighborsVSWeight(pointNumber), collapsedFractionOfCollapsedHyperedgeTwoNeighborsVSWeight(pointNumber), collapsedFractionOfCollapsedHyperedgeThreeNeighborsVSWeight(pointNumber);
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
    
    //Weight-weight correlation
    long double averageNeighborWeight=0, averageOneNeighborWeight=0, averageTwoNeighborWeight=0, averageThreeNeighborWeight=0;
    unsigned int binIndex=log(hyperedge.weight/(scale*averageWeight))/log(binSize);
    if (neighborVector.size()) {
      for (auto &neighbor : neighborVector) {averageNeighborWeight+=neighbor->weight;} averageNeighborWeight/=neighborVector.size();
      get<0>(weightNeighborsWeightCorrelation[binIndex])+=averageNeighborWeight/averageWeight;
      get<1>(weightNeighborsWeightCorrelation[binIndex])+=(averageNeighborWeight/averageWeight)*(averageNeighborWeight/averageWeight);
      ++get<2>(weightNeighborsWeightCorrelation[binIndex]);
    }
    if (oneNeighborVector.size()) {
      for (auto &oneNeighbor : oneNeighborVector) {averageOneNeighborWeight+=oneNeighbor->weight;} averageOneNeighborWeight/=oneNeighborVector.size();
      get<0>(weightOneNeighborsWeightCorrelation[binIndex])+=averageOneNeighborWeight/averageWeight;
      get<1>(weightOneNeighborsWeightCorrelation[binIndex])+=(averageOneNeighborWeight/averageWeight)*(averageOneNeighborWeight/averageWeight);
      ++get<2>(weightOneNeighborsWeightCorrelation[binIndex]);
    }
    if (twoNeighborVector.size()) {
      for (auto &twoNeighbor : twoNeighborVector) {averageTwoNeighborWeight+=twoNeighbor->weight;} averageTwoNeighborWeight/=twoNeighborVector.size();
      get<0>(weightTwoNeighborsWeightCorrelation[binIndex])+=averageTwoNeighborWeight/averageWeight;
      get<1>(weightTwoNeighborsWeightCorrelation[binIndex])+=(averageTwoNeighborWeight/averageWeight)*(averageTwoNeighborWeight/averageWeight);
      ++get<2>(weightTwoNeighborsWeightCorrelation[binIndex]);
    }
    if (threeNeighborVector.size()) {
      for (auto &threeNeighbor : threeNeighborVector) {averageThreeNeighborWeight+=threeNeighbor->weight;} averageThreeNeighborWeight/=threeNeighborVector.size();
      get<0>(weightThreeNeighborsWeightCorrelation[binIndex])+=averageThreeNeighborWeight/averageWeight;
      get<1>(weightThreeNeighborsWeightCorrelation[binIndex])+=(averageThreeNeighborWeight/averageWeight)*(averageThreeNeighborWeight/averageWeight);
      ++get<2>(weightThreeNeighborsWeightCorrelation[binIndex]);
    }
    
    //Collapsed fraction of hyperedge neighbors
    long double collapsedFractionOfHyperedgeNeighbors=0, collapsedFractionOfHyperedgeOneNeighbors=0, collapsedFractionOfHyperedgeTwoNeighbors=0, collapsedFractionOfHyperedgeThreeNeighbors=0;
    if (neighborVector.size()) {
      for (auto &neighbor : neighborVector) {
        if (neighbor->logRatio<=logRatioThreshold) {++collapsedFractionOfHyperedgeNeighbors;}
      } collapsedFractionOfHyperedgeNeighbors/=neighborVector.size();
    }
    if (oneNeighborVector.size()) {
      for (auto &oneNeighbor : oneNeighborVector) {
        if (oneNeighbor->logRatio<=logRatioThreshold) {++collapsedFractionOfHyperedgeOneNeighbors;}
      } collapsedFractionOfHyperedgeOneNeighbors/=oneNeighborVector.size();
    }
    if (twoNeighborVector.size()) {
      for (auto &twoNeighbor : twoNeighborVector) {
        if (twoNeighbor->logRatio<=logRatioThreshold) {++collapsedFractionOfHyperedgeTwoNeighbors;}
      } collapsedFractionOfHyperedgeTwoNeighbors/=twoNeighborVector.size();
    }
    if (threeNeighborVector.size()) {
      for (auto &threeNeighbor : threeNeighborVector) {
        if (threeNeighbor->logRatio<=logRatioThreshold) {++collapsedFractionOfHyperedgeThreeNeighbors;}
      } collapsedFractionOfHyperedgeThreeNeighbors/=threeNeighborVector.size();
    }
    if (hyperedge.logRatio<=logRatioThreshold) {
      if (neighborVector.size()) {
        get<0>(collapsedFractionOfCollapsedHyperedgeNeighbors)+=collapsedFractionOfHyperedgeNeighbors;
        get<1>(collapsedFractionOfCollapsedHyperedgeNeighbors)+=collapsedFractionOfHyperedgeNeighbors*collapsedFractionOfHyperedgeNeighbors;
        ++get<2>(collapsedFractionOfCollapsedHyperedgeNeighbors);
      }
      if (oneNeighborVector.size()) {
        get<0>(collapsedFractionOfCollapsedHyperedgeOneNeighbors)+=collapsedFractionOfHyperedgeOneNeighbors;
        get<1>(collapsedFractionOfCollapsedHyperedgeOneNeighbors)+=collapsedFractionOfHyperedgeOneNeighbors*collapsedFractionOfHyperedgeOneNeighbors;
        ++get<2>(collapsedFractionOfCollapsedHyperedgeOneNeighbors);
      }
      if (twoNeighborVector.size()) {
        get<0>(collapsedFractionOfCollapsedHyperedgeTwoNeighbors)+=collapsedFractionOfHyperedgeTwoNeighbors;
        get<1>(collapsedFractionOfCollapsedHyperedgeTwoNeighbors)+=collapsedFractionOfHyperedgeTwoNeighbors*collapsedFractionOfHyperedgeTwoNeighbors;
        ++get<2>(collapsedFractionOfCollapsedHyperedgeTwoNeighbors);
      }
      if (threeNeighborVector.size()) {
        get<0>(collapsedFractionOfCollapsedHyperedgeThreeNeighbors)+=collapsedFractionOfHyperedgeThreeNeighbors;
        get<1>(collapsedFractionOfCollapsedHyperedgeThreeNeighbors)+=collapsedFractionOfHyperedgeThreeNeighbors*collapsedFractionOfHyperedgeThreeNeighbors;
        ++get<2>(collapsedFractionOfCollapsedHyperedgeThreeNeighbors);
      }
    } else {
      if (neighborVector.size()) {
        get<0>(collapsedFractionOfNormalHyperedgeNeighbors)+=collapsedFractionOfHyperedgeNeighbors;
        get<1>(collapsedFractionOfNormalHyperedgeNeighbors)+=collapsedFractionOfHyperedgeNeighbors*collapsedFractionOfHyperedgeNeighbors;
        ++get<2>(collapsedFractionOfNormalHyperedgeNeighbors);
      }
      if (oneNeighborVector.size()) {
        get<0>(collapsedFractionOfNormalHyperedgeOneNeighbors)+=collapsedFractionOfHyperedgeOneNeighbors;
        get<1>(collapsedFractionOfNormalHyperedgeOneNeighbors)+=collapsedFractionOfHyperedgeOneNeighbors*collapsedFractionOfHyperedgeOneNeighbors;
        ++get<2>(collapsedFractionOfNormalHyperedgeOneNeighbors);
      }
      if (twoNeighborVector.size()) {
        get<0>(collapsedFractionOfNormalHyperedgeTwoNeighbors)+=collapsedFractionOfHyperedgeTwoNeighbors;
        get<1>(collapsedFractionOfNormalHyperedgeTwoNeighbors)+=collapsedFractionOfHyperedgeTwoNeighbors*collapsedFractionOfHyperedgeTwoNeighbors;
        ++get<2>(collapsedFractionOfNormalHyperedgeTwoNeighbors);
      }
      if (threeNeighborVector.size()) {
        get<0>(collapsedFractionOfNormalHyperedgeThreeNeighbors)+=collapsedFractionOfHyperedgeThreeNeighbors;
        get<1>(collapsedFractionOfNormalHyperedgeThreeNeighbors)+=collapsedFractionOfHyperedgeThreeNeighbors*collapsedFractionOfHyperedgeThreeNeighbors;
        ++get<2>(collapsedFractionOfNormalHyperedgeThreeNeighbors);
      }
    }
    
    //Collapsed fraction of hyperedge neighbors VS weight
    if (hyperedge.logRatio<=logRatioThreshold) {
      get<0>(collapsedFractionOfCollapsedHyperedgeNeighborsVSWeight[binIndex])+=collapsedFractionOfHyperedgeNeighbors;
      get<1>(collapsedFractionOfCollapsedHyperedgeNeighborsVSWeight[binIndex])+=collapsedFractionOfHyperedgeNeighbors*collapsedFractionOfHyperedgeNeighbors;
      ++get<2>(collapsedFractionOfCollapsedHyperedgeNeighborsVSWeight[binIndex]);
      get<0>(collapsedFractionOfCollapsedHyperedgeOneNeighborsVSWeight[binIndex])+=collapsedFractionOfHyperedgeOneNeighbors;
      get<1>(collapsedFractionOfCollapsedHyperedgeOneNeighborsVSWeight[binIndex])+=collapsedFractionOfHyperedgeOneNeighbors*collapsedFractionOfHyperedgeOneNeighbors;
      ++get<2>(collapsedFractionOfCollapsedHyperedgeOneNeighborsVSWeight[binIndex]);
      get<0>(collapsedFractionOfCollapsedHyperedgeTwoNeighborsVSWeight[binIndex])+=collapsedFractionOfHyperedgeTwoNeighbors;
      get<1>(collapsedFractionOfCollapsedHyperedgeTwoNeighborsVSWeight[binIndex])+=collapsedFractionOfHyperedgeTwoNeighbors*collapsedFractionOfHyperedgeTwoNeighbors;
      ++get<2>(collapsedFractionOfCollapsedHyperedgeTwoNeighborsVSWeight[binIndex]);
      get<0>(collapsedFractionOfCollapsedHyperedgeThreeNeighborsVSWeight[binIndex])+=collapsedFractionOfHyperedgeThreeNeighbors;
      get<1>(collapsedFractionOfCollapsedHyperedgeThreeNeighborsVSWeight[binIndex])+=collapsedFractionOfHyperedgeThreeNeighbors*collapsedFractionOfHyperedgeThreeNeighbors;
      ++get<2>(collapsedFractionOfCollapsedHyperedgeThreeNeighborsVSWeight[binIndex]);
    } else {
      get<0>(collapsedFractionOfNormalHyperedgeNeighborsVSWeight[binIndex])+=collapsedFractionOfHyperedgeNeighbors;
      get<1>(collapsedFractionOfNormalHyperedgeNeighborsVSWeight[binIndex])+=collapsedFractionOfHyperedgeNeighbors*collapsedFractionOfHyperedgeNeighbors;
      ++get<2>(collapsedFractionOfNormalHyperedgeNeighborsVSWeight[binIndex]);
      get<0>(collapsedFractionOfNormalHyperedgeOneNeighborsVSWeight[binIndex])+=collapsedFractionOfHyperedgeOneNeighbors;
      get<1>(collapsedFractionOfNormalHyperedgeOneNeighborsVSWeight[binIndex])+=collapsedFractionOfHyperedgeOneNeighbors*collapsedFractionOfHyperedgeOneNeighbors;
      ++get<2>(collapsedFractionOfNormalHyperedgeOneNeighborsVSWeight[binIndex]);
      get<0>(collapsedFractionOfNormalHyperedgeTwoNeighborsVSWeight[binIndex])+=collapsedFractionOfHyperedgeTwoNeighbors;
      get<1>(collapsedFractionOfNormalHyperedgeTwoNeighborsVSWeight[binIndex])+=collapsedFractionOfHyperedgeTwoNeighbors*collapsedFractionOfHyperedgeTwoNeighbors;
      ++get<2>(collapsedFractionOfNormalHyperedgeTwoNeighborsVSWeight[binIndex]);
      get<0>(collapsedFractionOfNormalHyperedgeThreeNeighborsVSWeight[binIndex])+=collapsedFractionOfHyperedgeThreeNeighbors;
      get<1>(collapsedFractionOfNormalHyperedgeThreeNeighborsVSWeight[binIndex])+=collapsedFractionOfHyperedgeThreeNeighbors*collapsedFractionOfHyperedgeThreeNeighbors;
      ++get<2>(collapsedFractionOfNormalHyperedgeThreeNeighborsVSWeight[binIndex]);
    }
  }
  
  //Weight-weight correlation
  ofstream printFileWeightNeighborWeightCorrelation("06. Weight-weight correlation/WeightNeighborsWeightCorrelation"+to_string(year)+".txt");
  printFileWeightNeighborWeightCorrelation << "w/<w>" << "\t" << "w_nei/<w>" << "\t" << "<X^2>" << "\n";
  for (unsigned int binIndex=0; binIndex<pointNumber; ++binIndex) {
    printFileWeightNeighborWeightCorrelation << pow(binSize, binIndex+0.5)*scale << "\t" << get<0>(weightNeighborsWeightCorrelation[binIndex])/get<2>(weightNeighborsWeightCorrelation[binIndex]) << "\t" << get<1>(weightNeighborsWeightCorrelation[binIndex])/get<2>(weightNeighborsWeightCorrelation[binIndex]) << "\n";
  } printFileWeightNeighborWeightCorrelation.close();
  ofstream printFileWeightOneNeighborWeightCorrelation("06. Weight-weight correlation/WeightOneNeighborsWeightCorrelation"+to_string(year)+".txt");
  printFileWeightOneNeighborWeightCorrelation << "w/<w>" << "\t" << "w_nei/<w>" << "\t" << "<X^2>" << "\n";
  for (unsigned int binIndex=0; binIndex<pointNumber; ++binIndex) {
    printFileWeightOneNeighborWeightCorrelation << pow(binSize, binIndex+0.5)*scale << "\t" << get<0>(weightOneNeighborsWeightCorrelation[binIndex])/get<2>(weightOneNeighborsWeightCorrelation[binIndex]) << "\t" << get<1>(weightOneNeighborsWeightCorrelation[binIndex])/get<2>(weightOneNeighborsWeightCorrelation[binIndex]) << "\n";
  } printFileWeightOneNeighborWeightCorrelation.close();
  ofstream printFileWeightTwoNeighborWeightCorrelation("06. Weight-weight correlation/WeightTwoNeighborsWeightCorrelation"+to_string(year)+".txt");
  printFileWeightTwoNeighborWeightCorrelation << "w/<w>" << "\t" << "w_nei/<w>" << "\t" << "<X^2>" << "\n";
  for (unsigned int binIndex=0; binIndex<pointNumber; ++binIndex) {
    printFileWeightTwoNeighborWeightCorrelation << pow(binSize, binIndex+0.5)*scale << "\t" << get<0>(weightTwoNeighborsWeightCorrelation[binIndex])/get<2>(weightTwoNeighborsWeightCorrelation[binIndex]) << "\t" << get<1>(weightTwoNeighborsWeightCorrelation[binIndex])/get<2>(weightTwoNeighborsWeightCorrelation[binIndex]) << "\n";
  } printFileWeightTwoNeighborWeightCorrelation.close();
  ofstream printFileWeightThreeNeighborWeightCorrelation("06. Weight-weight correlation/WeightThreeNeighborsWeightCorrelation"+to_string(year)+".txt");
  printFileWeightThreeNeighborWeightCorrelation << "w/<w>" << "\t" << "w_nei/<w>" << "\t" << "<X^2>" << "\n";
  for (unsigned int binIndex=0; binIndex<pointNumber; ++binIndex) {
    printFileWeightThreeNeighborWeightCorrelation << pow(binSize, binIndex+0.5)*scale << "\t" << get<0>(weightThreeNeighborsWeightCorrelation[binIndex])/get<2>(weightThreeNeighborsWeightCorrelation[binIndex]) << "\t" << get<1>(weightThreeNeighborsWeightCorrelation[binIndex])/get<2>(weightThreeNeighborsWeightCorrelation[binIndex]) << "\n";
  } printFileWeightThreeNeighborWeightCorrelation.close();
  
  //Collapsed fraction of hyperedge neighbors
  ofstream printFileCollapsedFractionOfHyperedgeNeighbors("09. Collapsed fraction of hyperedge neighbors/"+stringLogRatioThreshold+"/CollapsedFractionOfHyperedgeNeighbors.txt", ios::app);
  printFileCollapsedFractionOfHyperedgeNeighbors << year << "\t" << get<0>(collapsedFractionOfNormalHyperedgeNeighbors)/get<2>(collapsedFractionOfNormalHyperedgeNeighbors) << "\t" << get<1>(collapsedFractionOfNormalHyperedgeNeighbors)/get<2>(collapsedFractionOfNormalHyperedgeNeighbors) << "\t";
  printFileCollapsedFractionOfHyperedgeNeighbors << get<0>(collapsedFractionOfCollapsedHyperedgeNeighbors)/get<2>(collapsedFractionOfCollapsedHyperedgeNeighbors) << "\t" << get<1>(collapsedFractionOfCollapsedHyperedgeNeighbors)/get<2>(collapsedFractionOfCollapsedHyperedgeNeighbors) << "\n";
  printFileCollapsedFractionOfHyperedgeNeighbors.close();
  ofstream printFileCollapsedFractionOfHyperedgeOneNeighbors("09. Collapsed fraction of hyperedge neighbors/"+stringLogRatioThreshold+"/CollapsedFractionOfHyperedgeOneNeighbors.txt", ios::app);
  printFileCollapsedFractionOfHyperedgeOneNeighbors << year << "\t" <<  get<0>(collapsedFractionOfNormalHyperedgeOneNeighbors)/get<2>(collapsedFractionOfNormalHyperedgeOneNeighbors) << "\t" << get<1>(collapsedFractionOfNormalHyperedgeOneNeighbors)/get<2>(collapsedFractionOfNormalHyperedgeOneNeighbors) << "\t";
  printFileCollapsedFractionOfHyperedgeOneNeighbors << get<0>(collapsedFractionOfCollapsedHyperedgeOneNeighbors)/get<2>(collapsedFractionOfCollapsedHyperedgeOneNeighbors) << "\t" << get<1>(collapsedFractionOfCollapsedHyperedgeOneNeighbors)/get<2>(collapsedFractionOfCollapsedHyperedgeOneNeighbors) << "\n";
  printFileCollapsedFractionOfHyperedgeOneNeighbors.close();
  ofstream printFileCollapsedFractionOfHyperedgeTwoNeighbors("09. Collapsed fraction of hyperedge neighbors/"+stringLogRatioThreshold+"/CollapsedFractionOfHyperedgeTwoNeighbors.txt", ios::app);
  printFileCollapsedFractionOfHyperedgeTwoNeighbors << year << "\t" <<  get<0>(collapsedFractionOfNormalHyperedgeTwoNeighbors)/get<2>(collapsedFractionOfNormalHyperedgeTwoNeighbors) << "\t" << get<1>(collapsedFractionOfNormalHyperedgeTwoNeighbors)/get<2>(collapsedFractionOfNormalHyperedgeTwoNeighbors) << "\t";
  printFileCollapsedFractionOfHyperedgeTwoNeighbors << get<0>(collapsedFractionOfCollapsedHyperedgeTwoNeighbors)/get<2>(collapsedFractionOfCollapsedHyperedgeTwoNeighbors) << "\t" << get<1>(collapsedFractionOfCollapsedHyperedgeTwoNeighbors)/get<2>(collapsedFractionOfCollapsedHyperedgeTwoNeighbors) << "\n";
  printFileCollapsedFractionOfHyperedgeTwoNeighbors.close();
  ofstream printFileCollapsedFractionOfHyperedgeThreeNeighbors("09. Collapsed fraction of hyperedge neighbors/"+stringLogRatioThreshold+"/CollapsedFractionOfHyperedgeThreeNeighbors.txt", ios::app);
  printFileCollapsedFractionOfHyperedgeThreeNeighbors << year << "\t" <<  get<0>(collapsedFractionOfNormalHyperedgeThreeNeighbors)/get<2>(collapsedFractionOfNormalHyperedgeThreeNeighbors) << "\t" << get<1>(collapsedFractionOfNormalHyperedgeThreeNeighbors)/get<2>(collapsedFractionOfNormalHyperedgeThreeNeighbors) << "\t";
  printFileCollapsedFractionOfHyperedgeThreeNeighbors << get<0>(collapsedFractionOfCollapsedHyperedgeThreeNeighbors)/get<2>(collapsedFractionOfCollapsedHyperedgeThreeNeighbors) << "\t" << get<1>(collapsedFractionOfCollapsedHyperedgeThreeNeighbors)/get<2>(collapsedFractionOfCollapsedHyperedgeThreeNeighbors) << "\n";
  printFileCollapsedFractionOfHyperedgeThreeNeighbors.close();
  
  //Collapsed fraction of hyperedge neighbors VS weight
  ofstream printFileCollapsedFractionOfHyperedgeNeighborsVSWeight("10. Collapsed fraction of hyperedge neighbors VS weight/"+stringLogRatioThreshold+"/CollapsedFractionOfHyperedgeNeighborsVSWeight"+to_string(year)+".txt");
  printFileCollapsedFractionOfHyperedgeNeighborsVSWeight << "w/<w>" << "\t" << "Normal" << "\t" << "<X^2>" << "\t" << "Collapsed" << "\t" << "<X^2>" << "\n";
  for (unsigned int binIndex=0; binIndex<pointNumber; ++binIndex) {
    printFileCollapsedFractionOfHyperedgeNeighborsVSWeight << pow(binSize, binIndex+0.5)*scale << "\t";
    printFileCollapsedFractionOfHyperedgeNeighborsVSWeight << get<0>(collapsedFractionOfNormalHyperedgeNeighborsVSWeight[binIndex])/get<2>(collapsedFractionOfNormalHyperedgeNeighborsVSWeight[binIndex]) << "\t" << get<1>(collapsedFractionOfNormalHyperedgeNeighborsVSWeight[binIndex])/get<2>(collapsedFractionOfNormalHyperedgeNeighborsVSWeight[binIndex]) << "\t";
    printFileCollapsedFractionOfHyperedgeNeighborsVSWeight << get<0>(collapsedFractionOfCollapsedHyperedgeNeighborsVSWeight[binIndex])/get<2>(collapsedFractionOfCollapsedHyperedgeNeighborsVSWeight[binIndex]) << "\t" << get<1>(collapsedFractionOfCollapsedHyperedgeNeighborsVSWeight[binIndex])/get<2>(collapsedFractionOfCollapsedHyperedgeNeighborsVSWeight[binIndex]) << "\n";
  } printFileCollapsedFractionOfHyperedgeNeighborsVSWeight.close();
  ofstream printFileCollapsedFractionOfHyperedgeOneNeighborsVSWeight("10. Collapsed fraction of hyperedge neighbors VS weight/"+stringLogRatioThreshold+"/CollapsedFractionOfHyperedgeOneNeighborsVSWeight"+to_string(year)+".txt");
  printFileCollapsedFractionOfHyperedgeOneNeighborsVSWeight << "w/<w>" << "\t" << "Normal" << "\t" << "<X^2>" << "\t" << "Collapsed" << "\t" << "<X^2>" << "\n";
  for (unsigned int binIndex=0; binIndex<pointNumber; ++binIndex) {
    printFileCollapsedFractionOfHyperedgeOneNeighborsVSWeight << pow(binSize, binIndex+0.5)*scale << "\t";
    printFileCollapsedFractionOfHyperedgeOneNeighborsVSWeight << get<0>(collapsedFractionOfNormalHyperedgeOneNeighborsVSWeight[binIndex])/get<2>(collapsedFractionOfNormalHyperedgeOneNeighborsVSWeight[binIndex]) << "\t" << get<1>(collapsedFractionOfNormalHyperedgeOneNeighborsVSWeight[binIndex])/get<2>(collapsedFractionOfNormalHyperedgeOneNeighborsVSWeight[binIndex]) << "\t";
    printFileCollapsedFractionOfHyperedgeOneNeighborsVSWeight << get<0>(collapsedFractionOfCollapsedHyperedgeOneNeighborsVSWeight[binIndex])/get<2>(collapsedFractionOfCollapsedHyperedgeOneNeighborsVSWeight[binIndex]) << "\t" << get<1>(collapsedFractionOfCollapsedHyperedgeOneNeighborsVSWeight[binIndex])/get<2>(collapsedFractionOfCollapsedHyperedgeOneNeighborsVSWeight[binIndex]) << "\n";
  } printFileCollapsedFractionOfHyperedgeOneNeighborsVSWeight.close();
  ofstream printFileCollapsedFractionOfHyperedgeTwoNeighborsVSWeight("10. Collapsed fraction of hyperedge neighbors VS weight/"+stringLogRatioThreshold+"/CollapsedFractionOfHyperedgeTwoNeighborsVSWeight"+to_string(year)+".txt");
  printFileCollapsedFractionOfHyperedgeTwoNeighborsVSWeight << "w/<w>" << "\t" << "Normal" << "\t" << "<X^2>" << "\t" << "Collapsed" << "\t" << "<X^2>" << "\n";
  for (unsigned int binIndex=0; binIndex<pointNumber; ++binIndex) {
    printFileCollapsedFractionOfHyperedgeTwoNeighborsVSWeight << pow(binSize, binIndex+0.5)*scale << "\t";
    printFileCollapsedFractionOfHyperedgeTwoNeighborsVSWeight << get<0>(collapsedFractionOfNormalHyperedgeTwoNeighborsVSWeight[binIndex])/get<2>(collapsedFractionOfNormalHyperedgeTwoNeighborsVSWeight[binIndex]) << "\t" << get<1>(collapsedFractionOfNormalHyperedgeTwoNeighborsVSWeight[binIndex])/get<2>(collapsedFractionOfNormalHyperedgeTwoNeighborsVSWeight[binIndex]) << "\t";
    printFileCollapsedFractionOfHyperedgeTwoNeighborsVSWeight << get<0>(collapsedFractionOfCollapsedHyperedgeTwoNeighborsVSWeight[binIndex])/get<2>(collapsedFractionOfCollapsedHyperedgeTwoNeighborsVSWeight[binIndex]) << "\t" << get<1>(collapsedFractionOfCollapsedHyperedgeTwoNeighborsVSWeight[binIndex])/get<2>(collapsedFractionOfCollapsedHyperedgeTwoNeighborsVSWeight[binIndex]) << "\n";
  } printFileCollapsedFractionOfHyperedgeTwoNeighborsVSWeight.close();
  ofstream printFileCollapsedFractionOfHyperedgeThreeNeighborsVSWeight("10. Collapsed fraction of hyperedge neighbors VS weight/"+stringLogRatioThreshold+"/CollapsedFractionOfHyperedgeThreeNeighborsVSWeight"+to_string(year)+".txt");
  printFileCollapsedFractionOfHyperedgeThreeNeighborsVSWeight << "w/<w>" << "\t" << "Normal" << "\t" << "<X^2>" << "\t" << "Collapsed" << "\t" << "<X^2>" << "\n";
  for (unsigned int binIndex=0; binIndex<pointNumber; ++binIndex) {
    printFileCollapsedFractionOfHyperedgeThreeNeighborsVSWeight << pow(binSize, binIndex+0.5)*scale << "\t";
    printFileCollapsedFractionOfHyperedgeThreeNeighborsVSWeight << get<0>(collapsedFractionOfNormalHyperedgeThreeNeighborsVSWeight[binIndex])/get<2>(collapsedFractionOfNormalHyperedgeThreeNeighborsVSWeight[binIndex]) << "\t" << get<1>(collapsedFractionOfNormalHyperedgeThreeNeighborsVSWeight[binIndex])/get<2>(collapsedFractionOfNormalHyperedgeThreeNeighborsVSWeight[binIndex]) << "\t";
    printFileCollapsedFractionOfHyperedgeThreeNeighborsVSWeight << get<0>(collapsedFractionOfCollapsedHyperedgeThreeNeighborsVSWeight[binIndex])/get<2>(collapsedFractionOfCollapsedHyperedgeThreeNeighborsVSWeight[binIndex]) << "\t" << get<1>(collapsedFractionOfCollapsedHyperedgeThreeNeighborsVSWeight[binIndex])/get<2>(collapsedFractionOfCollapsedHyperedgeThreeNeighborsVSWeight[binIndex]) << "\n";
  } printFileCollapsedFractionOfHyperedgeThreeNeighborsVSWeight.close();
}
