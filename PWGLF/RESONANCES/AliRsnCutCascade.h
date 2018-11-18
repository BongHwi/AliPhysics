#ifndef ALIRSNCUTCASCADE_H
#define ALIRSNCUTCASCADE_H

// Class AliRsnCutCascade
//
// General Cascade cut for Xi/Omega related Analysis.
// Based on AliRsnCutV0, it contains every cuts for the CutV0 as prefix "V0"
//
// authors: Bong-Hwi Lim (bong-hwi.lim@cern.ch)
//          advised by Beomkyu Kim(kimb@cern.ch)
//

#include <TMath.h>
#include <TString.h>

#include "AliRsnCut.h"
#include "AliPIDResponse.h"
#include "AliRsnCutTrackQuality.h"

class AliESDtrack;
class AliAODTrack;
class AliRsnCutCascade : public AliRsnCut {
    public:
    
    AliRsnCutCascade(const char *name = "AliRsnCutCascade", Int_t hypothesis = kXiMinus, AliPID::EParticleType pid = AliPID::kProton, AliPID::EParticleType pid2 = AliPID::kPion, AliPID::EParticleType pid3 = AliPID::kPion);
    AliRsnCutCascade(const AliRsnCutCascade &copy);
    AliRsnCutCascade &operator=(const AliRsnCutCascade &copy);
    virtual ~AliRsnCutCascade() { }
    
    void           SetESDtrackCuts(AliESDtrackCuts *cuts)   {fESDtrackCuts = cuts;}
    void           SetHypothesis(Int_t code);
    void           SetXiTolerance(Double_t value)             {fXiTolerance = value;}
    void           SetXiToleranceVeto(Double_t value)         {fXiToleranceVeto = value;}
    void           SetSwitch(Bool_t value)                  {fSwitch = value;}
    
    void           SetV0fLowRadius(Double_t value)            {fV0LowRadius = value;}
    void           SetV0fHighRadius(Double_t value)           {fV0HighRadius = value;}
    
    void           SetXifLowRadius(Double_t value)            {fXiLowRadius = value;}
    void           SetXifHighRadius(Double_t value)           {fXiHighRadius = value;}
    
    void           SetV0MinDCAVertex(Double_t value)          {fV0MinDCAVertex = value;}
    void           SetXiMinDCAVertex(Double_t value)          {fXiMinDCAVertex = value;}
    void           SetV0MaxDCAVertex(Double_t value)          {fV0MaxDCAVertex = value;}
    void           SetXiMaxDCAVertex(Double_t value)          {fXiMaxDCAVertex = value;}
    
    void           SetV0MinCosPointingAngle(Double_t value)   {fV0MinCosPointAngle = value;}
    void           SetXiMinCosPointingAngle(Double_t value)   {fXiMinCosPointAngle = value;}
    
    void           SetV0MaxDaughtersDCA(Double_t value)       {fV0MaxDaughtersDCA = value;}
    void           SetXiMaxDaughtersDCA(Double_t value)       {fXiMaxDaughtersDCA = value;}
    
    void           SetMinTPCcluster(Int_t value)            {fMinTPCcluster = value;}
    void           SetMaxRapidity(Double_t value)           {fMaxRapidity = value;}
    
    void           SetPIDCutV0Proton(Double_t value)        {fPIDCutV0Proton = value;}
    void           SetPIDCutV0Pion(Double_t value)          {fPIDCutV0Proton = value;}
    void           SetPIDCutBPion(Double_t value)           {fPIDCutBPion = value;}
    
    
    AliRsnCutTrackQuality *CutQuality()                     {return &fCutQuality;}
    void           SetAODTestFilterBit(Int_t value)         {fAODTestFilterBit = value;}
    Int_t          GetAODTestFilterBit()                    {return fAODTestFilterBit;}
    
    virtual Bool_t IsSelected(TObject *obj);
    virtual void   Print(const Option_t *option = "") const;
    
protected:
    
    Bool_t      CheckESD(AliESDcascade *track);
    Bool_t      CheckAOD(AliAODcascade *track);
    
    
    Int_t            fHypothesis;       // PDG code corresponding to expected V0 hypothesis
    Double_t         fMass;             // mass corresponding to hypothesis
    Double_t         fXiTolerance;        // tolerance in the difference between computed and expected mass for Cascade
    Double_t         fXiToleranceVeto;    // Competing Cascade Rejection. Read the note in AliRsnCutCascade.cxx for more info.
    Bool_t           fSwitch;           // Switch for using Competing V0 Rejection
    Double_t         fV0LowRadius;        // Lower Limit on Fiducial Volume for V0
    Double_t         fV0HighRadius;       // Higher Limit on Fiducial Volume for V0
    Double_t         fXiLowRadius;        // Lower Limit on Fiducial Volume for Cascade
    Double_t         fXiHighRadius;       // Higher Limit on Fiducial Volume for Cascade
    Double_t         fV0MinDCAVertex;     // min allowed DCA from primary vertex of V0
    Double_t         fXiMinDCAVertex;     // min allowed DCA from primary vertex of Cascade
    Double_t         fV0MaxDCAVertex;     // max allowed DCA from primary vertex of V0
    Double_t         fXiMaxDCAVertex;     // max allowed DCA from primary vertex of Cascade
    Double_t         fV0MinCosPointAngle; // min allowed cosine of pointing angle of V0
    Double_t         fXiMinCosPointAngle; // min allowed cosine of pointing angle of Cascade
    Double_t         fV0MaxDaughtersDCA;  // max allowed DCA between the two daughers
    Double_t         fXiMaxDaughtersDCA;  // max allowed DCA between the two daughers
    Int_t            fMinTPCcluster;    // min allowed TOC cluster
    Double_t         fMaxRapidity;      // max allowed V0 rapidity
    
    AliPID::EParticleType fPID;         // PID for track
    AliPID::EParticleType fPID2;        // PID for track
    AliPID::EParticleType fPID3;        // PID for track
    
    Double_t         fPIDCutV0Proton;       // nsigmas for V0 proton
    Double_t         fPIDCutV0Pion;         // nsigmas for V0 pion
    Double_t         fPIDCutBPion;          // nsigmas for Bachelor pion
    
    AliESDtrackCuts *fESDtrackCuts;     // quality cuts for v0 daughters
    
    AliRsnCutTrackQuality fCutQuality;  // track quality cut
    
    Int_t            fAODTestFilterBit; // test filter bit for AODs
    
    ClassDef(AliRsnCutCascade, 1)
};

//__________________________________________________________________________________________________
inline void AliRsnCutCascade::SetHypothesis(Int_t code)
{
    //
    // Assign a Cascade species hypothesis, which also assign the expected mass
    //
    
    fHypothesis = code;
    
    switch (fHypothesis) {
        case kXiMinus:
            fMass = 1.3217;
            break;
        case kXiPlusBar:
            fMass = 1.3217;
            break;
        case kOmegaMinus:
            fMass = 1.6725;
            break;
        case kOmegaPlusBar:
            fMass = 1.6725;
            break;
        default:
            AliError(Form("You are setting an unexpected hypothesis: %d", fHypothesis));
            fMass = 1E20;
    }
}

#endif
