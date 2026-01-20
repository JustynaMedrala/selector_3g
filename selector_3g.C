#include "selector_3g.h"
#include "TMath.h"
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <random>
#include <map>


Selector3g::Selector3g(const std::string& inputFile,
                       const std::string& outputFile,
                       bool isMC)
    : fInput(inputFile), fOutput(outputFile), fisMC(isMC), fIsSignalEvent(false), fIsSignalTriplet(false), fIsSignalWithPrompt(false)
{
    fIn = TFile::Open(fInput.c_str(),"READ");
    if (!fIn || !fIn->IsOpen()) {
        std::cerr << "Nie udało się otworzyć pliku\n";
        return;
    }

    fTree = (TTree*)fIn->Get("T");
    fTree->SetBranchAddress("nhits",&nhits);
    fTree->SetBranchAddress("times",&times);
    fTree->SetBranchAddress("pos",&pos);
    fTree->SetBranchAddress("tots",&tots);
    fTree->SetBranchAddress("scins",&scins);

    const char* mcBranches[] = {"gammaTags","VtxIndices","MCIndices","MCTimes","MCEnergy","MCZ"};
    fisMC = true;
    for (const auto& name : mcBranches) {
        if (!fTree->GetBranch(name)) { fisMC = false; break; }
    }

    if (fisMC) {
        fTree->SetBranchAddress("gammaTags",&gammaTags);
        fTree->SetBranchAddress("VtxIndices",&MCVtxIndices);
        fTree->SetBranchAddress("MCIndices",&MCIndices);
        fTree->SetBranchAddress("MCTimes",&MCTimes);
        fTree->SetBranchAddress("MCEnergy",&MCEnergies);
        fTree->SetBranchAddress("MCZ",&MCZ);
    }

    fOut = TFile::Open(fOutput.c_str(),"RECREATE");
    fOutputTree = new TTree("HitTree","3γ reconstructed data");

    fOutputTree->Branch("eventID",&fEventID);
    fOutputTree->Branch("multiplicity",&nhits);
    fOutputTree->Branch("nAnnihilationTriplets",&fNAnnihilationTriplets);
    fOutputTree->Branch("nPromptHits",&fNPromptHits);
    fOutputTree->Branch("hasPromptHit",&fHasPromptHit);

    fOutputTree->Branch("hit1_x", &fHit1X); fOutputTree->Branch("hit1_y", &fHit1Y); fOutputTree->Branch("hit1_z", &fHit1Z);
    fOutputTree->Branch("hit1_time", &fHit1Time); fOutputTree->Branch("hit1_tot", &fHit1ToT); fOutputTree->Branch("hit1_energy", &fHit1Energy);
    
    fOutputTree->Branch("hit2_x", &fHit2X); fOutputTree->Branch("hit2_y", &fHit2Y); fOutputTree->Branch("hit2_z", &fHit2Z);
    fOutputTree->Branch("hit2_time", &fHit2Time); fOutputTree->Branch("hit2_tot", &fHit2ToT); fOutputTree->Branch("hit2_energy", &fHit2Energy);
    
    fOutputTree->Branch("hit3_x", &fHit3X); fOutputTree->Branch("hit3_y", &fHit3Y); fOutputTree->Branch("hit3_z", &fHit3Z);
    fOutputTree->Branch("hit3_time", &fHit3Time); fOutputTree->Branch("hit3_tot", &fHit3ToT); fOutputTree->Branch("hit3_energy", &fHit3Energy);
    
    fOutputTree->Branch("hit1_mcTime", &fHit1MCTime); fOutputTree->Branch("hit2_mcTime", &fHit2MCTime); fOutputTree->Branch("hit3_mcTime", &fHit3MCTime);
    fOutputTree->Branch("hit1_mcZ", &fHit1MCZ); fOutputTree->Branch("hit2_mcZ", &fHit2MCZ); fOutputTree->Branch("hit3_mcZ", &fHit3MCZ);

    fOutputTree->Branch("prompt_x", &fPromptX); fOutputTree->Branch("prompt_y", &fPromptY); fOutputTree->Branch("prompt_z", &fPromptZ);
    fOutputTree->Branch("prompt_time", &fPromptTime); fOutputTree->Branch("prompt_tot", &fPromptToT);

    fOutputTree->Branch("prompt_mcZ", &fPromptMCZ); fOutputTree->Branch("prompt_mcEnergy", &fPromptMCEnergy); fOutputTree->Branch("prompt_mcTime", &fPromptMCTime);

    fOutputTree->Branch("lifetime", &fLifetime); 

    fOutputTree->Branch("angle12", &fAngle12); fOutputTree->Branch("angle23", &fAngle23); fOutputTree->Branch("angle31", &fAngle31);
    fOutputTree->Branch("angleDiff", &fAngleDiff); fOutputTree->Branch("angleSum", &fAngleSum);
    fOutputTree->Branch("timeDiff", &fTimeDiff);
}

Selector3g::Selector3g(const std::string& inputFile)
    : Selector3g(inputFile,"T",false) {}

Selector3g::~Selector3g() {
    if(fIn) fIn->Close();
    if(fOut) fOut->Close();
}

void Selector3g::runAnalysis() {
    Long64_t entries = fTree->GetEntries();
    for(Long64_t i=0;i<entries;i++)
        processEvent(i);

    //printCutAnalysis();
    fOut->cd();
    fOutputTree->Write();
    fOut->Close();
}

void Selector3g::processEvent(Long64_t entry) {
    fTree->GetEntry(entry);
    fEventID = entry;

    fNAnnihilationTriplets = 0;
    fNPromptHits = 0;
    fHasPromptHit = false;

    if (nhits < 4) return; 

    std::vector<int> goodHits;
    for (Int_t i = 0; i < nhits; ++i)
        if (checkZ(i, fMaxZ)) goodHits.push_back(i);

    if (goodHits.size() < 4) return;

    bool isSignalEvent = false;
    if (fisMC) isSignalEvent = isTrueSignalEvent();
    updateCutStatistics("02_Multiplicity_ge_4", isSignalEvent);
    updateCutStatistics("03_Z_cut", isSignalEvent);

    std::vector<int> prompts;
    identifyPromptHits(goodHits, prompts);

    std::vector<std::vector<int>> triplets;
    identifyAnnihilationHits(goodHits, triplets);
    if (triplets.empty()) return;

    fTripletGeometryOK_global = false;
    fTripletTimeOK_global = false;
    fTripletDOPOK_global = false;
    fTripletGeometryOK_signal = false;
    fTripletTimeOK_signal = false;
    fTripletDOPOK_signal = false;
    hasGoodTriplet = false;
    hasGoodTriplet_signal = false;

    struct TripletInfo {
        int i1, i2, i3;
        double DOP;
        bool passGeometry;
        bool passTime;
        bool passDOP;
    };

    std::vector<TripletInfo> passingTriplets;

    for (auto& t : triplets) {
        int i1 = t[0], i2 = t[1], i3 = t[2];

        bool isSignalTriplet = false;
        if (fisMC)
            isSignalTriplet = isTrueSignalWithPrompt(i1, i2, i3, prompts);

        fTripletGeometryOK = getRelativeAngle(i1, i2, i3)>fMinRelAngleCut;
        fTripletTimeOK = getTimeDiff(i1, i2, i3)<fMaxTimeDiff;
        fTripletDOPOK = getDOP(i1, i2, i3)<=fMaxDOP;

        if (fTripletGeometryOK){
            if(isSignalTriplet) fTripletGeometryOK_signal = true; //signal triplet that passed geometry cut 
            fTripletGeometryOK_global = true; //at least one triplet per event passed geometry cut
            }
        if (fTripletTimeOK){
            if(isSignalTriplet) fTripletTimeOK_signal = true;
             fTripletTimeOK_global = true;}
        if (fTripletDOPOK){
            if(isSignalTriplet) fTripletDOPOK_signal = true; 
            fTripletDOPOK_global = true;}

        if (fTripletGeometryOK && fTripletTimeOK && fTripletDOPOK){
            ++fNAnnihilationTriplets;
            hasGoodTriplet = true;
            if(isSignalTriplet) hasGoodTriplet_signal = true;
            passingTriplets.push_back({i1,i2,i3,fTripletDOP,fTripletGeometryOK,fTripletTimeOK,fTripletDOPOK});
        }

        // Fill tree for this triplet
        //fillTripletInfo(i1, i2, i3);
        //fOutputTree->Fill();
    }

    if (!passingTriplets.empty()&&!prompts.empty()) {
        auto best = *std::min_element(passingTriplets.begin(), passingTriplets.end(),
                                      [](const TripletInfo& a, const TripletInfo& b){ return a.DOP < b.DOP; });

        int i1 = best.i1, i2 = best.i2, i3 = best.i3;
        fTripletGeometryOK = best.passGeometry;
        fTripletTimeOK     = best.passTime;
        fTripletDOPOK      = best.passDOP;
        fTripletDOP        = best.DOP;

        fNPromptHits = prompts.size();
        fLifetime = -1000000;
        fHasPromptHit = false;


        isSignal = false;
        if (fisMC) isSignal = isTrueSignalWithPrompt(i1,i2,i3,prompts);
        updateCutStatistics("08_At_least_one_prompt", isSignal);

        if (prompts.size() == 1) {
            int p = prompts[0];
            if (fisMC) isSignal = isTrueSignalWithPrompt(i1,i2,i3,{p});
            fillPromptInfo(p);
            fillTripletInfo(i1,i2,i3);
            fHasPromptHit = true;
            fLifetime = ((*times)[i1] + (*times)[i2] + (*times)[i3])/3. - (*times)[p];
            updateCutStatistics("09_Prompt", isSignal);
        } else {
            fPromptX = fPromptY = fPromptZ = fPromptTime = fPromptToT =
            fPromptMCEnergy = fPromptMCTime = fPromptMCZ = -1000000;
        }

        fOutputTree->Fill();
    }

    if(fTripletGeometryOK_global) updateCutStatistics("04_Geometry", fTripletGeometryOK_signal);
    if(fTripletTimeOK_global) updateCutStatistics("05_Time", fTripletTimeOK_signal);
    if(fTripletDOPOK_global) updateCutStatistics("06_DOP", fTripletDOPOK_signal);
    if(hasGoodTriplet) updateCutStatistics("07_Good_Triplet", hasGoodTriplet_signal);
}

bool Selector3g::checkZ(int idx,double maxZ) const {
    return ((*pos)[idx].Z() >= -maxZ && (*pos)[idx].Z() <= maxZ);
}

bool Selector3g::checkToT(int idx,double min,double max) const {
    double tot = (*tots)[idx];
    return (tot>=min && tot<=max);
}

double Selector3g::calculateTOF(int idx,const TVector3& source) const {
    double dist = ((*pos)[idx]-source).Mag();
    return dist/Constants::kLightVelocity_cm_ns;
}

double Selector3g::calculateAngle(int i1,int i2) {
    return (*pos)[i1].Angle((*pos)[i2])*180./TMath::Pi();
}

TVector3 Selector3g::calculateAnnihilationPoint(int i1,int i2,int i3) const {
    return ((*pos)[i1]+(*pos)[i2]+(*pos)[i3])*(1./3.);
}

void Selector3g::identifyPromptHits(const std::vector<int>& goodHits,std::vector<int>& prompts) const {
    for (auto i : goodHits) {
        if (checkToT(i, fTotCutDeexMin, fTotCutDeexMax))
            prompts.push_back(i);
    }
}

void Selector3g::identifyAnnihilationHits(const std::vector<int>& goodHits,std::vector<std::vector<int>>& triplets) const {
    for (size_t a=0;a<goodHits.size();a++)
        for (size_t b=a+1;b<goodHits.size();b++)
            for (size_t c=b+1;c<goodHits.size();c++) {
                int i=goodHits[a], j=goodHits[b], k=goodHits[c];
                if (checkToT(i,fTotCutAnniMin,fTotCutAnniMax) &&
                    checkToT(j,fTotCutAnniMin,fTotCutAnniMax) &&
                    checkToT(k,fTotCutAnniMin,fTotCutAnniMax)) {
                    triplets.push_back({i,j,k});
                }
            }
}

void Selector3g::updateCutStatistics(const std::string& name,bool passed){
    fCutStats[name].totalEvents++;
    if(passed) fCutStats[name].signalEvents++;
}

void Selector3g::printCutAnalysis(){
    std::cout<<"\nCUT SUMMARY:\n";
    int signal_total = fCutStats.begin()->second.signalEvents;
    for(auto& c : fCutStats)
        std::cout<<std::setw(30)<<c.first
                 <<"  passed="<<c.second.totalEvents
                 <<"  signal="<<c.second.signalEvents
                 <<"  efficiency="<<std::fixed
                 <<c.second.signalEvents/(double)signal_total
                 <<"  purity="<<std::fixed
                 <<c.second.signalEvents/(double)c.second.totalEvents<<"\n";
}

double Selector3g::getTimeDiff(int i1, int i2, int i3) const {
    double t1 = (*times)[i1] - calculateTOF(i1, fSourcePos);
    double t2 = (*times)[i2] - calculateTOF(i2, fSourcePos);
    double t3 = (*times)[i3] - calculateTOF(i3, fSourcePos);

    double dt12 = fabs(t1 - t2);
    double dt23 = fabs(t2 - t3);
    double dt31 = fabs(t3 - t1);

    return std::max({dt12, dt23, dt31});
}

double Selector3g::getRelativeAngle(int i1, int i2, int i3) const {
    TVector3 p1 = (*pos)[i1];
    TVector3 p2 = (*pos)[i2];
    TVector3 p3 = (*pos)[i3];

    double a12 = p1.Angle(p2) * 180.0 / TMath::Pi();
    double a23 = p2.Angle(p3) * 180.0 / TMath::Pi();
    double a31 = p3.Angle(p1) * 180.0 / TMath::Pi();

    std::vector<double> angles = {a12, a23, a31};
    std::sort(angles.begin(), angles.end());
    return angles[0] + angles[1];
}

double Selector3g::getDOP(int i1, int i2, int i3) const {
    TVector3 p1 = (*pos)[i1];
    TVector3 p2 = (*pos)[i2];
    TVector3 p3 = (*pos)[i3];

    TVector3 ap = calculateAnnihilationPoint(i1, i2, i3);
    TVector3 n = (p2 - p1).Cross(p3 - p1);

    if (n.Mag() < 1e-9) return 1e9;
    return fabs(n.Dot(ap - p1)) / n.Mag();
}

bool Selector3g::isTrueSignalTriplet(int i1,int i2,int i3) const {
    bool isAnnihilation = (*gammaTags)[i1]==3 && (*gammaTags)[i2]==3 && (*gammaTags)[i3]==3;
    int v1 = (*MCVtxIndices)[i1], v2 = (*MCVtxIndices)[i2], v3 = (*MCVtxIndices)[i3];
    bool sameVertex = (v1==v2) && (v2==v3);
    //if(isAnnihilation) std::cout<<v1<<" "<<v2<<" "<<v3<<std::endl;
    return isAnnihilation && sameVertex;
}

bool Selector3g::isTrueSignalWithPrompt(int i1,int i2,int i3,const std::vector<int>& prompts) const {
    if (!isTrueSignalTriplet(i1,i2,i3)) return false;
    int vtx = (*MCVtxIndices)[i1];
    for (auto p : prompts)
        if ((*gammaTags)[p]==1 && (*MCVtxIndices)[p]==vtx)
            return true;
    return false;
}

bool Selector3g::isTrueSignalEvent() const {
    if (!fisMC) return false;
    std::map<int,std::pair<int,int>> vtxCounts;
    for (Int_t i=0;i<nhits;i++){
        int tag = (*gammaTags)[i];
        int vtx = (*MCVtxIndices)[i];
        if(tag==3) vtxCounts[vtx].first++;
        else if(tag==1) vtxCounts[vtx].second++;
    }
    for (auto& [vtx,counts] : vtxCounts)
        if(counts.first>=3 && counts.second>=1) return true;
    return false;
}

void Selector3g::fillTripletInfo(int i1, int i2, int i3)
{
    TVector3 p1 = (*pos)[i1];
    TVector3 p2 = (*pos)[i2];
    TVector3 p3 = (*pos)[i3];

    fHit1X = p1.X(); fHit1Y = p1.Y(); fHit1Z = p1.Z(); 
    fHit1Time = (*times)[i1]; 
    fHit1ToT = (*tots)[i1]; 

    fHit2X = p2.X(); fHit2Y = p2.Y(); fHit2Z = p2.Z(); 
    fHit2Time = (*times)[i2]; 
    fHit2ToT = (*tots)[i2]; 

    fHit3X = p3.X(); fHit3Y = p3.Y(); fHit3Z = p3.Z(); 
    fHit3Time = (*times)[i3]; 
    fHit3ToT = (*tots)[i3]; 

    if(fisMC) {
        fHit1MCZ = (*MCZ)[i1]; fHit1MCTime = (*MCTimes)[i1]; fHit1Energy = (*MCEnergies)[i1];
        fHit2MCZ = (*MCZ)[i2]; fHit2MCTime = (*MCTimes)[i2]; fHit2Energy = (*MCEnergies)[i2];
        fHit3MCZ = (*MCZ)[i3]; fHit3MCTime = (*MCTimes)[i3]; fHit3Energy = (*MCEnergies)[i3];}


    fAngle12 = calculateAngle(i1, i2);
    fAngle23 = calculateAngle(i2, i3);
    fAngle31 = calculateAngle(i3, i1);

    std::vector<double> ang = {fAngle12, fAngle23, fAngle31};
    std::sort(ang.begin(), ang.end());
    fAngleSum = ang[0] + ang[1];
    fAngleDiff = ang[1] - ang[0];

    TVector3 ap = calculateAnnihilationPoint(i1, i2, i3);
    fAnnihilationPointX = ap.X();
    fAnnihilationPointY = ap.Y();
    fAnnihilationPointZ = ap.Z();

    std::vector<double> tripletCopy = {(*times)[i1], (*times)[i2], (*times)[i3]};  // skopiuj dane
    std::shuffle(tripletCopy.begin(), tripletCopy.end(), std::mt19937{std::random_device{}()});

    fTimeDiff.clear();
    fTimeDiff.push_back(tripletCopy[0]- tripletCopy[1]);
    fTimeDiff.push_back(tripletCopy[1] - tripletCopy[2]);
    fTimeDiff.push_back(tripletCopy[2] - tripletCopy[0]);   
}

void Selector3g::fillPromptInfo(int idx)
{
    TVector3 p = (*pos)[idx];

    fPromptX = p.X();
    fPromptY = p.Y();
    fPromptZ = p.Z();

    fPromptTime = (*times)[idx];
    fPromptToT = (*tots)[idx];
    if(fisMC) {
        fPromptMCTime = (*MCTimes)[idx];
        fPromptMCEnergy = (*MCEnergies)[idx];
        fPromptMCZ = (*MCZ)[idx];}
}