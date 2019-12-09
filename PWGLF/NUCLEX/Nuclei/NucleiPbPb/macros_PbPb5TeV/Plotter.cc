#include "src/Common.h"
#include "src/Utils.h"
using namespace utils;
#include <TFile.h>
#include <RooPlot.h>
#include <TCanvas.h>
#include <TDirectory.h>
#include <iostream>

#include <string>
using std::string;

using namespace RooFit;

string output_name = "../results/plots.root";

void Plotter(bool short_mode = true){

  TFile file_in(kSignalOutput.data());
  TFile file_out(output_name.data(),"RECREATE");

  TCanvas* fCanvas= nullptr;
  int iPad = 0;
  int nPads = 0;

  int counter = 0;

  for (auto list_key : *file_in.GetListOfKeys()) {
    /// Preliminary operation to read the list and create an output dir
    if(short_mode)
      if(counter>0) return;
    counter++;
    if (string(list_key->GetName()).find(kFilterListNames.data()) == string::npos) continue;

    TTList* list = (TTList*)file_in.Get(list_key->GetName());
    TDirectory* base_dir = file_out.mkdir(list_key->GetName());
    file_out.cd(list_key->GetName());
    gSystem->Exec(Form("mkdir -p ../results/images/%s/",list_key->GetName()));
    for(int iS=0; iS<2; iS++){
      TDirectory* dir = base_dir->mkdir(kNames[iS].data());
      dir->cd();
      for(int iC=0; iC<kCentLength; iC++){
        dir->mkdir(Form("cent_%d",iC));
        dir->cd(Form("cent_%d",iC));
        iPad = 0;
        nPads = kPtBinLimit[iC] - 5 + 1;
        for(int iB=5; iB<=kPtBinLimit[iC]; iB++){
          std::cout << "iPad: " << iPad << ", iB: " << iB << std::endl;
          if(iPad%6 == 0){
            if(fCanvas) delete fCanvas;
            fCanvas = new TCanvas(Form("Canvas_%d",iPad/6),Form("Canvas_%d",iPad/6),kCanvasW,kCanvasH);
            fCanvas->Divide(3,2);
          }
          //if(iPad==0) fCanvas->Print(Form("cent_%d_%c.pdf[",iC,kLetter[iS]));
          fCanvas->cd(iPad%6+1);
          string path = string(list_key->GetName()) + "/" + kNames[iS] + "/TailTail/C_" + to_string(iC) ;
          RooPlot* fPlot = (RooPlot*)file_in.Get(Form("%s/d%i_%i",path.data(),iC,iB));
          // std::cout << "path " << Form("%s/d%i_%i",path.data(),iC,iB) << std::endl;
          Requires(fPlot,"RooPlot");
          fPlot->Draw();
          iPad++;
          if(iPad==nPads){
            fCanvas->Print(Form("../results/images/%s/cent_%d_%c_%d.pdf",list_key->GetName(),iC,kLetter[iS],iB/6));
            fCanvas->Write();
          }
          else if((iPad)%6 == 0){
            fCanvas->Print(Form("../results/images/%s/cent_%d_%c_%d.pdf",list_key->GetName(),iC,kLetter[iS],iB/6));
            fCanvas->Write();
          }
        }
      }
    }
    if (kTest) break; // for the test
  }
}
