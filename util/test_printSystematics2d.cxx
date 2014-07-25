#include "DileptonMatrixMethod/DileptonMatrixMethod.h"
#include "DileptonMatrixMethod/Lepton.h"

#include "TArrayD.h"

#include <iostream>
#include <iterator>
#include <string>
#include <cassert>
#include <cmath>
#include <vector>

/*
  Test the matrix method with the recommended, most up-to-date matrix.

  davide.gerbaudo@gmail.com
  April 2014
*/

using namespace std;
namespace sf = susy::fake;

//----------------------------------------------------------
vector<double>  bin_edges2bin_centers(const TArrayD* binEdges)
{
    vector<double> binCenters;
    if(binEdges){
        for(int i=0; i<binEdges->GetSize()-1; ++i)
            binCenters.push_back(0.5*(binEdges->At(i) + binEdges->At(i+1)));
    }
    return binCenters;
}
//----------------------------------------------------------
void printSystematics(sf::DileptonMatrixMethod & matrix, bool isEl, bool isMu)
{
    if(isEl==isMu) cout<<"pick either el or mu"<<endl;
    assert(isEl!=isMu);

    vector<double> ptBins (bin_edges2bin_centers(matrix.getPtBins()));
    vector<double> etaBins(bin_edges2bin_centers(matrix.getEtaBins()));
    cout<<"pt bin centers: ";
    copy(ptBins.begin(), ptBins.end(), ostream_iterator<double>(cout, ", "));
    cout<<endl;
    cout<<"eta bin centers: ";
    copy(etaBins.begin(), etaBins.end(), ostream_iterator<double>(cout, ", "));
    cout<<endl;
    bool isTight(false), isElectron(true);
    float gev2mev(1.0e3); // the matrix is stored in GeV, but the inputs are supposed to be in MeV...how nice
    const std::string regionName = "CR_SSInc1j";
    for(size_t i=0; i<ptBins.size(); ++i)
        for(size_t j=0; j<etaBins.size(); ++j){
            float pt(ptBins[i]), eta(etaBins[j]);
            sf::Lepton l(isTight, isElectron, pt*gev2mev, eta);
            sf::DileptonMatrixMethod::RATE_TYPE rt = sf::DileptonMatrixMethod::FAKE;
            matrix.printRateSystematics(l, rt, matrix.getIndexRegion(regionName));
        }
}
//----------------------------------------------------------
int main(int argc, char **argv)
{

  string inputFilename="data/FinalFakeHist_May_20.root";

  bool isEl(false), isMu(false);

  const sf::Parametrization::Value rp = sf::Parametrization::PT_ETA;
  std::vector<std::string> regions;
  regions.push_back("CR_SSInc1j");
  sf::DileptonMatrixMethod matrix;
  if(matrix.configure(inputFilename, regions, rp, rp, rp, rp)) {
      cout<<endl<<" --- Fake systematic uncertainties --- "<<endl;

      cout<<endl<<" --- Fake electron  --- "<<endl;
      isEl = true; isMu = false;
      printSystematics(matrix, isEl, isMu);

      isEl = false; isMu = true;
      cout<<endl<<" --- Fake muon      --- "<<endl;
      printSystematics(matrix, isEl, isMu);
      cout<<endl;
  }
  return 0;
}
