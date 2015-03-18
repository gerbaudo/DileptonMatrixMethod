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
void printSystematics1d(sf::DileptonMatrixMethod & matrix, bool isEl, bool isMu, const string &region)
{
    if(isEl==isMu) cout<<"pick either el or mu"<<endl;
    assert(isEl!=isMu);

    vector<double> ptBins (bin_edges2bin_centers(matrix.getPtBins()));
    cout<<"pt bin centers: ";
    copy(ptBins.begin(), ptBins.end(), ostream_iterator<double>(cout, ", "));
    cout<<endl;
    bool isTight(false), isElectron(isEl);
    float gev(1.0);
    float eta(1.5); // just a dummy value for the 1D param
    for(size_t i=0; i<ptBins.size(); ++i){
        float pt(ptBins[i]);
        sf::Lepton l(isTight, isElectron, pt*gev, eta);
        sf::DileptonMatrixMethod::RATE_TYPE rt = sf::DileptonMatrixMethod::FAKE;
        matrix.printRateSystematics(l, rt, matrix.getIndexRegion(region));
    }
}
//----------------------------------------------------------
void printSystematics2d(sf::DileptonMatrixMethod & matrix, bool isEl, bool isMu, const string &region)
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
    bool isTight(false), isElectron(isEl);
    float gev(1.0);
    for(size_t i=0; i<ptBins.size(); ++i)
        for(size_t j=0; j<etaBins.size(); ++j){
            float pt(ptBins[i]), eta(etaBins[j]);
            sf::Lepton l(isTight, isElectron, pt*gev, eta);
            sf::DileptonMatrixMethod::RATE_TYPE rt = sf::DileptonMatrixMethod::FAKE;
            matrix.printRateSystematics(l, rt, matrix.getIndexRegion(region));
        }
}
//----------------------------------------------------------
int main(int argc, char **argv)
{
    if(argc!=4){
        cout<<"usage : "<<argv[0]<<" FakeMatrix.root region_name pt"<<endl;
        cout<<"Example:"<<endl
            <<argv[0]<<" ${ROOTCOREBIN}/data/DileptonMatrixMethod/FakeMatrix_Dec_03_with_syst.root razor  pt"
            <<argv[0]<<" ${ROOTCOREBIN}/data/DileptonMatrixMethod/FakeMatrix_Mar_18.root           emuInc pteta"
            <<endl;
        return 0;
    }
    string inputFilename   = argv[1];
    string region          = argv[2];
    string parametrization = argv[3];
    sf::Parametrization::Value rp = sf::Parametrization::PT;
    if(parametrization=="pt") rp = sf::Parametrization::PT;
    else if(parametrization=="pteta") rp = sf::Parametrization::PT_ETA;
    else{
        cout<<"Invalid parametrization '"<<parametrization<<"', must be either pt or pteta"<<endl;
        return 0;
    }

  bool isEl(false), isMu(false);

  std::vector<std::string> regions;
  regions.push_back(region);
  sf::DileptonMatrixMethod matrix;
  if(matrix.configure(inputFilename, regions, rp, rp, rp, rp)) {
      cout<<endl<<" --- Fake systematic uncertainties --- "<<endl;

      cout<<endl<<" --- Fake electron  --- "<<endl;
      isEl = true; isMu = false;
      if(rp==sf::Parametrization::PT)
          printSystematics1d(matrix, isEl, isMu, region);
      else
          printSystematics2d(matrix, isEl, isMu, region);

      isEl = false; isMu = true;
      cout<<endl<<" --- Fake muon      --- "<<endl;
      if(rp==sf::Parametrization::PT)
          printSystematics1d(matrix, isEl, isMu, region);
      else
          printSystematics2d(matrix, isEl, isMu, region);
      cout<<endl;
  }
  return 0;
}
