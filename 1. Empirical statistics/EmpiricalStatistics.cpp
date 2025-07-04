#include"Algorithm.hpp"
using namespace std;

long double logRatioThreshold=-1; //-0.5, -1, -2, -100
string stringLogRatioThreshold="Threshold=-1";

int main() {
  ofstream printFileTheNumberOfNodesAndHyperedges("01. The number of nodes and hyperedges/TheNumberOfNodesAndHyperedges.txt");
  printFileTheNumberOfNodesAndHyperedges << "Year" << "\t" << "Countries" << "\t" << "Products" << "\t" << "Hyperedges" << "\n";
  printFileTheNumberOfNodesAndHyperedges.close();
  ofstream printFileAverageWeight("02. Weight distibutions/AverageWeight.txt");
  printFileAverageWeight << "Year" << "\t" << "<w>" << "\t" << "<w^2>" << "\n";
  printFileAverageWeight.close();
  ofstream printFileAverageDegree("03. Degree distributions/AverageDegree.txt");
  printFileAverageDegree << "Year" << "\t" << "<k>" << "\t" << "<k^2>" << "\n";
  printFileAverageDegree.close();
  ofstream printFileAverageFiniteLogRatio("04. Finite log ratio distributions/AverageFiniteLogRatio.txt");
  printFileAverageFiniteLogRatio << "Year" << "\t" << "<g>_f" << "\t" << "<g^2>_f" << "\n";
  printFileAverageFiniteLogRatio.close();
  ofstream printFileTheNumberOfCollapsedHyperedge("07. The number of collapsed hyperedges/"+stringLogRatioThreshold+"/TheNumberOfCollapsedHyperedges.txt");
  printFileTheNumberOfCollapsedHyperedge << "Year" << "\t" << "Hyperedges" << "\t" << "Infinitely collapsed" << "\t" << "finitely collapsed" << "\n";
  printFileTheNumberOfCollapsedHyperedge.close();
  ofstream printFileCollapsedFractionOfHyperedgeNeighbors("09. Collapsed fraction of hyperedge neighbors/"+stringLogRatioThreshold+"/CollapsedFractionOfHyperedgeNeighbors.txt");
  printFileCollapsedFractionOfHyperedgeNeighbors << "Year" << "\t" << "Normal" << "\t" << "<X^2>" << "\t" << "Collapsed" << "\t" << "<X^2>" << "\n";
  printFileCollapsedFractionOfHyperedgeNeighbors.close();
  ofstream printFileCollapsedFractionOfHyperedgeOneNeighbors("09. Collapsed fraction of hyperedge neighbors/"+stringLogRatioThreshold+"/CollapsedFractionOfHyperedgeOneNeighbors.txt");
  printFileCollapsedFractionOfHyperedgeOneNeighbors << "Year" << "\t" << "Normal" << "\t" << "<X^2>" << "\t" << "Collapsed" << "\t" << "<X^2>" << "\n";
  printFileCollapsedFractionOfHyperedgeOneNeighbors.close();
  ofstream printFileCollapsedFractionOfHyperedgeTwoNeighbors("09. Collapsed fraction of hyperedge neighbors/"+stringLogRatioThreshold+"/CollapsedFractionOfHyperedgeTwoNeighbors.txt");
  printFileCollapsedFractionOfHyperedgeTwoNeighbors << "Year" << "\t" << "Normal" << "\t" << "<X^2>" << "\t" << "Collapsed" << "\t" << "<X^2>" << "\n";
  printFileCollapsedFractionOfHyperedgeTwoNeighbors.close();
  ofstream printFileCollapsedFractionOfHyperedgeThreeNeighbors("09. Collapsed fraction of hyperedge neighbors/"+stringLogRatioThreshold+"/CollapsedFractionOfHyperedgeThreeNeighbors.txt");
  printFileCollapsedFractionOfHyperedgeThreeNeighbors << "Year" << "\t" << "Normal" << "\t" << "<X^2>" << "\t" << "Collapsed" << "\t" << "<X^2>" << "\n";
  printFileCollapsedFractionOfHyperedgeThreeNeighbors.close();
  for (unsigned int year=2004; year<=2018; ++year) {Hypergraph hypergraph(year);}
}
