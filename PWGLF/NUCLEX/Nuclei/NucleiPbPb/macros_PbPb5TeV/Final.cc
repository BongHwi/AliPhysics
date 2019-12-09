#include "src/Common.h"
#include "src/Plotting.h"
#include "src/Utils.h"
using namespace utils;

#include <map>
#include <vector>
#include <array>
using std::array;
#include <memory>

#include <Riostream.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TList.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TLine.h>
#include <TF1.h>

#include <AliPID.h>
#include <AliPWGFunc.h>

#include <RooRealVar.h>
#include <RooArgSet.h>
#include <RooDataSet.h>
#include <RooPlot.h>
#include <RooGaussian.h>
#include <TVirtualFitter.h>
#include <TPaveText.h>
#include <TLatex.h>

const char* kPrefix[2] = {"","anti"};

const char* kParticleNames[2] = {"Deuterons", "Antideuterons"};

void Final(const char* fitFunctionName = "Boltzmann") {
  TFile spectra_file(kSpectraOutput.data());
  TFile TOF_systematics_file(kSystematicsOutput.data());
  TFile TPC_systematics_file(kSystematicsOutputTPC.data());
  TFile final_file(kFinalOutput.data(),"recreate");

  // const double pt_bin_limits[16] = {0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.4,1.6,1.8,2.0,2.2,2.6,3.0,3.4,3.8};//,4.4};
  // const int kNPtBins = 15;

  const int n_centralities = kCentLength;

  TH1F* stat[2][n_centralities];
  TH1F* syst[2][n_centralities];
  TH1F* syst_corr[2][n_centralities];
  TH1F* syst_uncorr[2][n_centralities];
  TH1F* stat_tpc[2][n_centralities];
  TH1F* syst_tpc[2][n_centralities];
  // TH1F* syst_tpc_corr[2][kCentLength];
  // TH1F* syst_tpc_uncorr[2][kCentLength];

  for (int iS = 0; iS < 2; ++iS) {
    TDirectory* s_dir = final_file.mkdir(kNames[iS].data());
    for (int iC = 0; iC < n_centralities; ++iC) {
      TDirectory *c_dir = s_dir->mkdir(to_string(iC).data());
      c_dir->cd();
      string tof_basepath = kFilterListNames + "/" + kNames[iS] + "/TOFspectra" + to_string(iC);
      TH1F* spectra_tof_tmp  = (TH1F*)spectra_file.Get(tof_basepath.data());
      Requires(spectra_tof_tmp,tof_basepath.data());
      stat[iS][iC] = (TH1F*)spectra_tof_tmp->Clone(Form("stat_%d_%d",iS,iC));
      auto ptAxis = stat[iS][iC]->GetXaxis();

      TH1F* totsyst_tmp = (TH1F*)TOF_systematics_file.Get((to_string(iC) + "/" + kNames[iS] + "/totsyst").data());
      Requires(totsyst_tmp, "Missing totsyt");
      TH1F* totsyst = (TH1F*)totsyst_tmp->Clone(Form("totsyst_%d_%d",iS,iC));

      string corr_sytpath = (to_string(iC) + "/" + kNames[iS] + "/syst_corr").data();
      TH1F* corr_syst_rel = (TH1F*)TOF_systematics_file.Get(corr_sytpath.data());
      Requires(corr_syst_rel,corr_sytpath.data());
      syst_corr[iS][iC]  = (TH1F*)stat[iS][iC]->Clone(("syst_corr" + std::to_string(iC)).data());

      string uncorr_sytpath = (to_string(iC) + "/" + kNames[iS] + "/syst_uncorr").data();
      TH1F* uncorr_syst_rel = (TH1F*)TOF_systematics_file.Get(uncorr_sytpath.data());
      Requires(uncorr_syst_rel,uncorr_sytpath.data());
      syst_uncorr[iS][iC]  = (TH1F*)stat[iS][iC]->Clone(("syst_uncorr" + std::to_string(iC)).data());

      TH1F* totsyst_tpc_tmp = (TH1F*)TPC_systematics_file.Get((to_string(iC) + "/" + kNames[iS] + "/totsyst_tpc").data());
      Requires(totsyst_tpc_tmp, "Missing totsyt_tpc");
      TH1F* totsyst_tpc = (TH1F*)totsyst_tpc_tmp->Clone(Form("totsyst_tpc_%d_%d",iS,iC));
      
      string tpc_basepath = kFilterListNames + "/" + kNames[iS] + "/TPCspectra" + to_string(iC);
      TH1F* spectra_tpc_tmp  = (TH1F*)spectra_file.Get(tpc_basepath.data());
      Requires(spectra_tpc_tmp,tpc_basepath.data());
      stat_tpc[iS][iC] = (TH1F*)spectra_tpc_tmp->Clone(Form("stat_tpc_%d_%d",iS,iC));
      syst[iS][iC]  = (TH1F*)totsyst->Clone(("syst" + to_string(iC)).data());
      syst[iS][iC]->Reset();
      syst_tpc[iS][iC]  = (TH1F*)totsyst->Clone(("syst_tpc" + to_string(iC)).data());
      syst_tpc[iS][iC]->Reset();

      for (int iB = 1; iB <= kPtBinLimit[iC]; ++iB) {
        std::cout << "bin: " << ptAxis->GetBinUpEdge(iB) << ", limit: " << kCentPtLimits[iC] << std::endl;
        if (ptAxis->GetBinCenter(iB) < kPtRange[0]|| ptAxis->GetBinCenter(iB) > kCentPtLimits[iC]){
          continue;
        }
        else{
          syst[iS][iC]->SetBinContent(iB,stat[iS][iC]->GetBinContent(iB));
          syst[iS][iC]->SetBinError(iB,totsyst->GetBinContent(iB) * stat[iS][iC]->GetBinContent(iB));
          syst_corr[iS][iC]->SetBinError(iB,corr_syst_rel->GetBinContent(iB) * stat[iS][iC]->GetBinContent(iB));
          syst_uncorr[iS][iC]->SetBinError(iB,uncorr_syst_rel->GetBinContent(iB) * stat[iS][iC]->GetBinContent(iB));
        }
        // if(ptAxis->GetBinCenter(iB)<1.4){
        //   syst_tpc[iS][iC]->SetBinContent(iB,stat_tpc[iS][iC]->GetBinContent(iB));
        //   syst_tpc[iS][iC]->SetBinError(iB,totsyst_tpc->GetBinContent(iB) * stat_tpc[iS][iC]->GetBinContent(iB));
        // }
        // else{
        //   stat_tpc[iS][iC]->SetBinContent(iB,0.);
        //   stat_tpc[iS][iC]->SetBinError(iB,0.);
        //   syst_tpc[iS][iC]->SetBinContent(iB,0.);
        //   syst_tpc[iS][iC]->SetBinError(iB,0.);
        // }
      }
      stat[iS][iC]->Write("stat_tof_orig");
      syst[iS][iC]->Write("syst_tof_orig");
      syst_corr[iS][iC]->Write("syst_corr_orig");
      syst_uncorr[iS][iC]->Write("syst_uncorr_orig");
      plotting::SetHistStyle(stat[iS][iC],plotting::MaterialColors[iC]);
      plotting::SetHistStyle(syst[iS][iC],plotting::MaterialColors[iC]);
      plotting::SetHistStyle(syst_corr[iS][iC],plotting::MaterialColors[iC]);
      plotting::SetHistStyle(syst_uncorr[iS][iC],plotting::MaterialColors[iC]);
      stat[iS][iC]->Scale(kScaleFactor[iC]);
      stat[iS][iC]->Write("stat_tof");
      syst[iS][iC]->Scale(kScaleFactor[iC]);
      syst[iS][iC]->Write("syst_tof");

      syst_corr[iS][iC]->Scale(kScaleFactor[iC]);
      syst_corr[iS][iC]->Write("syst_corr");
      syst_uncorr[iS][iC]->Scale(kScaleFactor[iC]);
      syst_uncorr[iS][iC]->Write("syst_uncorr");
      // plotting::SetHistStyle(stat_tpc[iS][iC],plotting::MaterialColors[iC],21);
      // plotting::SetHistStyle(syst_tpc[iS][iC],plotting::MaterialColors[iC],21);
      // stat_tpc[iS][iC]->Scale(kScaleFactor[iC]);
      // stat_tpc[iS][iC]->Write("stat_tpc");
      // syst_tpc[iS][iC]->Scale(kScaleFactor[iC]);
      // syst_tpc[iS][iC]->Write("syst_tpc");

    }

    TCanvas spectra("spectra","spectra",800,600);
    TH1* hFrame = spectra.DrawFrame(
        kPtRange[0] - 0.1,
        // 0.9 * syst[iS][n_centralities-2]->GetMinimum(),
        1.1e-7,
        kPtBins[kPtBinLimit[0]]+0.2,
        // 1e2 * syst[iS][0]->GetMaximum(),
        0.9e2,
        ";#it{p}_{T} (GeV/#it{c});#frac{1}{#it{N}_{ev}} #frac{d^{2}#it{N}}{d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-1}"
        );
    spectra.SetLeftMargin(0.15);
    spectra.SetRightMargin(0.03);
    spectra.SetTopMargin(0.1);
    spectra.SetBottomMargin(0.14);
    hFrame->GetYaxis()->SetTitleOffset(1.3);
    TLatex text;
    text.SetTextFont(63);
    text.SetTextSize(22);
    //text.DrawText(0.5,0.75,"ALICE Preliminary");
    float name_position = (iS==0) ? 3.3 : 3.0;
    // text.DrawLatex(name_position,7.5,Form("#bf{%s, pp, #sqrt{#it{s}} = 13 TeV}",kNamePlot[iS].data()));
    // text.DrawLatex(3.35461,0.41997,Form("#bf{#color[%d]{#LTd#it{N}_{ch} / d#it{#eta}#GT = 26.02}}",plotting::kSpectraColors[0]));
    // text.DrawLatex(0.641832,1.02914e-06,Form("#bf{#color[%d]{#LTd#it{N}_{ch} / d#it{#eta}#GT = 2.55}}",plotting::kSpectraColors[8]));
    // TLegend final_leg(0.70,0.60,0.90,0.88);
    // TLegend final_leg(0.699248,0.233043,0.943609,0.558261);
    TLegend final_leg(0.17,0.15,0.9,0.35);
    final_leg.SetBorderSize(0);
    final_leg.SetFillColor(0);
    final_leg.SetTextSize(0.027);
    final_leg.SetHeader("Pb-Pb #sqrt{s_{NN}=5.02 TeV}");
    TLegendEntry *header = (TLegendEntry*)final_leg.GetListOfPrimitives()->First();
    header->SetTextAlign(22);
    final_leg.SetNColumns(3);
    TFile bwfile(Form("%s%sfits.root",kBaseOutputDir.data(),kPrefix[iS]),"read");
    if (bwfile.IsOpen()) {
      TF1* bw = nullptr;
      TH1F* scaledbw = nullptr;
      for (int iC = 0; iC < kCentLength; ++iC) {
        bw = (TF1*)bwfile.Get(Form("%s/%i/%s%i",fitFunctionName,iC,fitFunctionName,iC));
        Requires(bw, Form("%s/%i/%s%i",fitFunctionName,iC,fitFunctionName,iC));
        if (!bw) continue;
        scaledbw = new TH1F(Form("scaledbw%i",iC),"",1000,0.5,1.05*kCentPtLimits[iC]);
        scaledbw->Add(bw);
        scaledbw->Scale(kScaleFactor[iC]);
        scaledbw->SetLineStyle(kDashed);
        scaledbw->SetLineColor(kBlack);
        spectra.cd();
        scaledbw->Draw("lsame");
      }
      final_leg.AddEntry(scaledbw,"Individual fit","l");
    }
    for (int iC = 0; iC < kCentLength; ++iC) {
      stat[iS][iC]->Draw("esamex0");
      syst[iS][iC]->Draw("e2same");
      // stat_tpc[iS][iC]->Draw("esamex0");
      // syst_tpc[iS][iC]->Draw("e2same");
      final_leg.AddEntry(syst[iS][iC],Form("%4.0f - %2.0f %% (#times %.1f)",kCentLabels[iC][0],kCentLabels[iC][1],kScaleFactor[iC]),"fp");
    }
    final_leg.AddEntry(syst[iS][7],"TOF","p");
    // final_leg.AddEntry(syst_tpc[iS][7],"TPC","p");
    final_leg.Draw();
    spectra.SetLogy();
    s_dir->cd();
    spectra.Write();
     if (kPrintFigures) {
      spectra.SaveAs((kFiguresFolder + "spectraTOF" + kLetter[iS] + ".eps").data());
      spectra.SaveAs((kMacrosFolder + "spectraTOF" + kLetter[iS] + ".C").data());
    }
  }
  TDirectory* r_dir = final_file.mkdir("ratio");
  for (int iC = 0; iC < n_centralities; ++iC) {
    r_dir->mkdir(to_string(iC).data())->cd();
    stat[1][iC]->Divide(stat[0][iC]);
    syst[1][iC]->Divide(syst[0][iC]);
    for(int iB=1; iB<=kNPtBins; iB++){
      if(stat[1][iC]->GetBinCenter(iB)<1.){
        stat[1][iC]->SetBinContent(iB, 0.);
        stat[1][iC]->SetBinError(iB, 0.);
        syst[1][iC]->SetBinContent(iB, 0.);
        syst[1][iC]->SetBinError(iB, 0.);
      }
    }
    stat[1][iC]->Write("stat_tof");
    syst[1][iC]->Write("syst_tof");

    // stat_tpc[1][iC]->Divide(stat_tpc[0][iC]);
    // syst_tpc[1][iC]->Divide(syst_tpc[0][iC]);
    // for(int iB=1; iB<=kNPtBins; iB++){
    //   if(stat_tpc[1][iC]->GetBinCenter(iB)>1.2){
    //     stat_tpc[1][iC]->SetBinContent(iB, 0.);
    //     stat_tpc[1][iC]->SetBinError(iB, 0.);
    //     syst_tpc[1][iC]->SetBinContent(iB, 0.);
    //     syst_tpc[1][iC]->SetBinError(iB, 0.);
    //   }
    // }
    TCanvas ratio("ratio","ratio");
    ratio.DrawFrame(
        0.5 * kPtRange[0],
        0.1,
        1.05 * kCentPtLimits[iC],
        1.9,
        ";#it{p}_{T} (GeV/#it{c});#bar{d}/d"
        );
    stat[1][iC]->Draw("esamex0");
    syst[1][iC]->Draw("e2same");
    // stat_tpc[1][iC]->Draw("esamex0");
    // syst_tpc[1][iC]->Draw("e2same");
    TLine *line = new TLine(0.5 * kPtRange[0],1,1.05 * kPtRange[1],1);
    line->SetLineColor(kBlack);
    line->Draw();
    if (kPrintFigures) ratio.SaveAs((kFiguresFolder + "ratio.eps").data());
    ratio.Write();
  }
}