#ifndef SELECTOR_3G_H
#define SELECTOR_3G_H

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"

namespace Constants {
    constexpr double kLightVelocity_cm_ns = 29.9792458;
}

struct CutStats {
    int totalEvents = 0;
    int signalEvents = 0;
};

class Selector3g {
public:
    Selector3g(const std::string& inputFile,
               const std::string& outputFile);
    Selector3g(const std::string& inputFile); 
    ~Selector3g();

    void runAnalysis();
    void processEvent(Long64_t entry);

private:
    std::string fInput;
    std::string fOutput;
    bool fisMC; 
    bool fIsSignalEvent; 
    bool fIsSignalTriplet;
    bool fIsSignalWithPrompt;
    bool isSignal;
    TFile* fIn = nullptr;
    TTree* fTree = nullptr;

    TFile* fOut = nullptr;
    TTree* fOutputTree = nullptr;
    Long64_t fEventID;
    int fNAnnihilationTriplets;
    int fNPromptHits;
    bool fHasPromptHit;

    bool fTripletDOPOK;
    bool fTripletGeometryOK;
    bool fTripletTimeOK;

    bool fTripletGeometryOK_global;
    bool fTripletTimeOK_global;
    bool fTripletDOPOK_global;

    bool fTripletGeometryOK_signal;
    bool fTripletTimeOK_signal;
    bool fTripletDOPOK_signal;

    bool hasGoodTriplet;
    bool hasGoodTriplet_signal;

    Int_t nhits;
    std::vector<double>* times = nullptr;
    std::vector<TVector3>* pos = nullptr;
    std::vector<double>* tots = nullptr;
    std::vector<int>* scins = nullptr;

    std::vector<int>* gammaTags = nullptr;
    std::vector<int>* MCVtxIndices = nullptr;
    std::vector<int>* MCIndices = nullptr;
    std::vector<double>* MCTimes = nullptr;
    std::vector<double>* MCEnergies = nullptr;
    std::vector<double>* MCZ = nullptr;


    float fHit1X, fHit1Y, fHit1Z;
    float fHit1Time, fHit1ToT, fHit1Energy;

    float fHit2X, fHit2Y, fHit2Z;
    float fHit2Time, fHit2ToT, fHit2Energy;

    float fHit3X, fHit3Y, fHit3Z;
    float fHit3Time, fHit3ToT, fHit3Energy;

    double fHit1MCTime, fHit2MCTime, fHit3MCTime;
    double fHit1MCZ, fHit2MCZ, fHit3MCZ;

    float fAngle12, fAngle23, fAngle31;
    float fAngleSum, fAngleDiff;

    float fAnnihilationPointX, fAnnihilationPointY, fAnnihilationPointZ;

    std::vector<float> fTimeDiff;

    float fPromptX = -1000000;
    float fPromptY = -1000000;
    float fPromptZ = -1000000;
    float fPromptTime = -1000000;
    float fPromptToT = -1000000;
    float fPromptMCZ = -1000000;
    float fPromptMCTime = -1000000;
    float fPromptMCEnergy = -1000000;

    double fLifetime = -1000000;

    double fTripletDOP = 0.;

    double fMaxZ = 22.;
    //double fTotCutAnniMin = 0., fTotCutAnniMax = 350.;
    //double fTotCutDeexMin = 400., fTotCutDeexMax = 1500.;
    double fTotCutAnniMin = 0., fTotCutAnniMax = 77000000.;
    double fTotCutDeexMin = 8500000., fTotCutDeexMax = 150000000;
    double fMinRelAngleCut = 190.;
    double fMaxTimeDiff = 3.;
    double fMaxDOP = 8.;
    TVector3 fSourcePos{0,0,0};

    std::map<std::string,CutStats> fCutStats;

    bool checkZ(int idx,double maxZ) const;
    bool checkToT(int idx,double min,double max) const;
    double calculateTOF(int idx,const TVector3& source) const;
    double calculateAngle(int i1,int i2);
    TVector3 calculateAnnihilationPoint(int i1,int i2,int i3) const;

    void updateCutStatistics(const std::string& name,bool passed);
    void printCutAnalysis();

    void identifyPromptHits(const std::vector<int>& goodHits,std::vector<int>& prompts) const;
    void identifyAnnihilationHits(const std::vector<int>& goodHits,std::vector<std::vector<int>>& triplets) const;
    double getTimeDiff(int i,int j,int k) const;       
    double getRelativeAngle(int i,int j,int k) const;
    double getDOP(int i,int j,int k) const;    
    bool isTrueSignalTriplet(int i1, int i2, int i3) const;
    bool isTrueSignalWithPrompt(int i1, int i2, int i3, const std::vector<int>& prompts) const;
    bool isTrueSignalEvent() const;   
    void fillTripletInfo(int i1, int i2, int i3);
    void fillPromptInfo(int idx);
};

#endif // SELECTOR_3G_H
