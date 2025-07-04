#include<cmath>
#include<utility>
#include<vector>
#include<fstream>
using namespace std;

extern unsigned int year, fEnsembleNumber;
extern string correlation;
extern long double alpha, ilambda, dlambda, flambda, recoveryRate, outbreakThreshold;

class Point {
public:
  vector<long double> datumVector;
  Point() {datumVector.reserve(fEnsembleNumber);}
  void Add(long double datum) {datumVector.emplace_back(datum);}
  unsigned int SurvivalNum() {
    unsigned int count=0;
    for (auto &datum : datumVector) {if (datum>outbreakThreshold) {++count;}}
    return count;
  }
  long double SurvivalAvg(unsigned int exponent) {
    long double sum=0; unsigned int count=0;
    for (auto &datum : datumVector) {if (datum>outbreakThreshold) {sum+=pow(datum, exponent); ++count;}}
    return count ? sum/count : 0;
  }
};

class Data {
public:
  vector<Point> pointVector=vector<Point>(((flambda-ilambda)/dlambda)+1.5);
  void Add(unsigned int pointIndex, long double datum) {pointVector[pointIndex].Add(datum);}
  void Print(string filename) {
    ofstream printFile(filename); printFile.precision(9);
    printFile << "Year=" << year << "\tCorrelation=" << correlation << "\n";
    printFile << "alpha=" << alpha << "\tgamma=" << recoveryRate << "\tOutbreak threshold=" << outbreakThreshold << "\n";
    printFile << "Ensembles" << "\t" << "Samples" << "\t" << "lambda" << "\t" << "<X>" << "\t" << "<X^2>" << "\n";
    for (auto &point : pointVector) {
      unsigned int pointIndex=&point-&pointVector.front();
      printFile << point.datumVector.size() << "\t" << point.SurvivalNum() << "\t" << ilambda+dlambda*pointIndex << "\t";
      printFile << point.SurvivalAvg(1) << "\t"  << point.SurvivalAvg(2) << "\n";
    } printFile.close();
  }
};
