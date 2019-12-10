#include "src/Common.h"
#include "src/Utils.h"
#include "src/Plotting.h"
using namespace utils;

#include <TDirectoryFile.h>
#include <TEfficiency.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TH2F.h>
using namespace std;

void Compare(TFile* oldreult = new TFile("../final.root","read")){
    bool antimatter_analysys = false;
    const char* kind_of_particle = (antimatter_analysys) ? "anti" : "";
    TFile *input_file = TFile::Open(kFinalOutput.data());
    TFile output_file(kCompareOutput.data(),"recreate");

    TH1D *stat[kCentLength],*syst[kCentLength],*stat_old[kCentLength],*syst_old[kCentLength],*stat_ratio[kCentLength],*syst_ratio[kCentLength];
    
    for (int iC = 0; iC < kCentLength; ++iC) {
        stat[iC] = (TH1D*)input_file->Get(Form("%sdeuterons/%i/stat_tof_orig",kind_of_particle,iC));
        stat_old[iC] = (TH1D*)oldreult->Get(Form("%sdeuterons/%i/stat",kind_of_particle,iC));

        syst[iC] = (TH1D*)input_file->Get(Form("%sdeuterons/%i/syst_tof_orig",kind_of_particle,iC));
        syst_old[iC] = (TH1D*)oldreult->Get(Form("%sdeuterons/%i/syst",kind_of_particle,iC));
        
        stat_ratio[iC] = (TH1D*)stat_old[iC]->Rebin(20,"test",kPtBins_old);
        TH1* h2 = (TH1*)stat[iC]->Rebin(20,"test",kPtBins_old);

        syst_ratio[iC] = (TH1D*)syst_old[iC]->Rebin(20,"test",kPtBins_old);
        TH1* h2_sys = (TH1*)syst[iC]->Rebin(20,"test",kPtBins_old);

        stat_ratio[iC]->Divide(h2);
        syst_ratio[iC]->Divide(h2_sys);
        output_file.cd();
        // std::cout << stat_ratio[iC]->GetBinContent(9) << std::endl;
        plotting::SetHistStyle(stat_ratio[iC],plotting::MaterialColors[iC]);
        plotting::SetHistStyle(syst_ratio[iC],plotting::MaterialColors[iC]);

        stat_ratio[iC]->Write(Form("%sdeuterons_stat_ratio_%d",kind_of_particle,iC));
        syst_ratio[iC]->Write(Form("%sdeuterons_syst_ratio_%d",kind_of_particle,iC));
    }

    TLatex* t = new TLatex();
        t->SetNDC();
        t->SetTextSize(0.035);


    TDirectory* r_dir = output_file.mkdir("ratio");
    TCanvas* fCanvas = new TCanvas("Canvas","Canvas",kCanvasW,kCanvasH);
    fCanvas->Divide(5,2);
    for (int iC = 0; iC < kCentLength; ++iC) {
        r_dir->cd();
        TCanvas ratio("ratio","ratio");
        ratio.DrawFrame(
            0.5 * kPtRange[0],
            0.1,
            1.05 * kCentPtLimits[iC],
            1.9,
            ";#it{p}_{T} (GeV/#it{c});new/old"
            );
        fCanvas->DrawFrame(
            0.5 * kPtRange[0],
            0.1,
            1.05 * kCentPtLimits[iC],
            1.9,
            ";#it{p}_{T} (GeV/#it{c});new/old"
            );
        stat_ratio[iC]->Draw("esamex0");
        syst_ratio[iC]->Draw("e2same");
        TLine *line = new TLine(0.5 * kPtRange[0],1,1.05 * kPtRange[1],1);
        line->SetLineColor(kBlack);
        line->Draw();
        TF1* pol0fit = new TF1("basicline", "[0]", 0, 8);
        syst_ratio[iC]->Fit(pol0fit,"","",1.5,8);
        // std::cout << pol0fit->GetParameter(0) << std::endl;
        auto diff = abs(1-pol0fit->GetParameter(0))*1e2;
        std::cout << diff << std::endl;

        t->DrawLatex(0.15, 0.8, Form("#bf{Difference: %.1f}",diff));


        ratio.Write(Form("%sdeuterons_ratio_%d",kind_of_particle,iC));

        fCanvas->cd(iC+1);
        stat_ratio[iC]->GetXaxis()->SetRangeUser(1.0,6.0);
        stat_ratio[iC]->SetMaximum(1.5);
        stat_ratio[iC]->SetMinimum(0.5);
        stat_ratio[iC]->Draw("esamex0");
        syst_ratio[iC]->Draw("e2same");
        line->Draw();
    }
    fCanvas->Write(Form("%sdeuterons_ratio",kind_of_particle));

}