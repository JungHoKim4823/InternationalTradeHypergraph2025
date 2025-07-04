#include<random>
#include<cmath>
#include<utility>
#include<vector>
#include<map>
#include<fstream>
#include<iostream>
using namespace std;

extern unsigned int year;
extern const string correlation;
extern long double alpha;
extern const long double ilambda, dlambda, flambda;

class Country; class Product; class Hyperedge;
class Country {public: vector<Hyperedge*> exporterHyperedgeVector, importerHyperedgeVector, countryHyperedgeVector;};
class Product {public: vector<Hyperedge*> productHyperedgeVector;};
class Hyperedge {
public:
  Country *exporter, *importer; Product *product;
  long double weight, alphaWeight, logRatio; unsigned int degree;
};

class Hypergraph {
public:
  map<unsigned int, Country> countryMap; map<unsigned int, Product> productMap; vector<Hyperedge> hyperedgeVector;
  long double outbreakFraction, averageDegree=0, averageDegreeSquared=0, averageAlphaWeight=0;
  
  Hypergraph();
};

Hypergraph::Hypergraph() {
  //Construct a trade hypergraph
  ifstream scanFile("UN"+to_string(year)+".txt");
  hyperedgeVector.reserve(1e6);
  unsigned int exporter, importer, product; long double weight, logRatio;
  while (scanFile >> exporter >> importer >> product >> weight >> logRatio) {
    Hyperedge hyperedge=Hyperedge();
    hyperedge.exporter=&countryMap[exporter];
    hyperedge.importer=&countryMap[importer];
    hyperedge.product=&productMap[product];
    hyperedge.weight=weight;
    hyperedge.logRatio=logRatio;
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
  
  //Weight to the power of alpha
  for (auto &hyperedge : hyperedgeVector) {
    hyperedge.alphaWeight=pow(hyperedge.weight, alpha);
    averageAlphaWeight+=hyperedge.alphaWeight;
  } averageAlphaWeight/=hyperedgeVector.size();
  
  //Weight-dependent mean field approximation
  ofstream printFile("Solution"+to_string(year)+"_"+correlation+"_alpha="+to_string(alpha).substr(0,6)+".txt"); printFile.precision(9);
  printFile << "lambda" << "\t" << "Outbreak Fraction" << "\n";
  for (long double lambda=ilambda; lambda<=flambda; lambda+=dlambda) {
    long double cInfTilde=1, oldcInfTilde=0;
    while (abs(cInfTilde-oldcInfTilde)>1e-12) {
      oldcInfTilde=cInfTilde; cInfTilde=0;
      for (auto &hyperedge : hyperedgeVector) {
        const long double numerator=lambda*averageDegree*hyperedge.degree*hyperedge.alphaWeight*oldcInfTilde;
        const long double denominator=averageDegreeSquared*averageAlphaWeight;
        cInfTilde+=hyperedge.degree*(1-exp(-(numerator/denominator)))/(averageDegree*hyperedgeVector.size());
      }
    }
    long double cInf=0;
    for (auto &hyperedge : hyperedgeVector) {
      const long double numerator=lambda*averageDegree*hyperedge.degree*hyperedge.alphaWeight*cInfTilde;
      const long double denominator=averageDegreeSquared*averageAlphaWeight;
      cInf+=(1-exp(-(numerator/denominator)))/hyperedgeVector.size();
    }
    printFile << lambda << "\t" << cInf << "\n";
    cout << year << "\t" << lambda << "\t" << cInf << "\n";
  } printFile.close();
}
