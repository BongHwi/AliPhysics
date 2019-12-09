#ifndef PLOTTING_H
#define PLOTTING_H

#include <array>

#include <TCanvas.h>
#include <TColor.h>
#include <TH1.h>
#include <TStyle.h>

namespace plotting {
  constexpr std::array<int,7> kHighContrastColors{
    kRed,
    kOrange,
    kGreen+1,
    kAzure+1,
    kBlue,
    kViolet+2,
    kBlack
  };
  vector<int> MaterialColors = {kOrange + 10, kOrange - 3, kOrange,
                              kTeal + 2,    kCyan + 2,    kAzure + 7,  kAzure - 1,
                              kViolet - 1,  kPink + 8, kBlack};
                              
  const std::array<int,11> kSpectraColors{
    TColor::GetColor("#ff3300"),
    TColor::GetColor("#ec6e0a"),
    TColor::GetColor("#daaa14"),
    TColor::GetColor("#c7e51e"),
    TColor::GetColor("#85dd69"),
    TColor::GetColor("#42d6b4"),
    TColor::GetColor("#00ceff"),
    kBlack,
    TColor::GetColor("#009adf"),
    TColor::GetColor("#0067c0"),
    TColor::GetColor("#0033a1")
  };

  const Style_t kChargeMarker[2] = {20,24};

  void SetHistStyle(TH1* h, int color, int marker = 20, int linew = 1, int fillstyle = 0) {
    h->SetMarkerColor(color);
    h->SetLineColor(color);
    h->SetFillStyle(fillstyle);
    h->SetLineWidth(linew);
    h->SetMarkerStyle(marker);
  }
}

#endif
