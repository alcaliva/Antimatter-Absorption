#include "AliAnalysisTaskJulia.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisUtils.h"
#include "AliMultSelection.h"
#include "AliMultEstimator.h"
#include "AliAnalysisTask.h"
#include "AliPIDResponse.h"
#include "TLorentzVector.h"
#include "AliCentrality.h"
#include "AliEventplane.h"
#include "AliEventCuts.h"
#include "TDatabasePDG.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODEvent.h"
#include "TObjArray.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TVectorD.h"
#include "AliAODv0.h"
#include "TRandom.h"
#include "TChain.h"
#include "TMath.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TF1.h"

ClassImp(AliAnalysisTaskJulia)

//_________________________________________________________________________________________________________________________________________________________________________________________________
AliAnalysisTaskJulia::AliAnalysisTaskJulia():
AliAnalysisTaskSE(),
fAODevent(NULL),
fPIDResponse(NULL),
fAODeventCuts(),
fUtils(NULL),
fOutputList(NULL),
fQAList(NULL),
fCentralityMin(0),
fCentralityMax(0),
fVertexZmin(0),
fVertexZmax(0),
fNumberOfVertexContributorsMin(0),
fMagFieldSign(0),
fDCAxyMax(0),
fDCAzMax(0),
fPtMin(0),
fPtMax(0),
fEtaMax(0),
fNumberOfClustersITSMin(0),
fNumberOfClustersTPCMin(0),
fNumberOfClustersTRDMin(0),
fNumberOfClustersTPCdEdxMin(0),
fNumberOfCrossedRowsMin(0),
fCrossedRowsOverFindableClsMin(0),
fChiSquarePerNDFMax(0),
fChiSquareTRDMax(0),
fNumberOfTRDtrackletsPIDMin(0),
fIsTRDrefit(0),
fTPCnsigmaDeuteronMax(0),
fTOFnsigmaDeuteronMax(0)
{}
//_________________________________________________________________________________________________________________________________________________________________________________________________
AliAnalysisTaskJulia::AliAnalysisTaskJulia(const char *name):
AliAnalysisTaskSE(name),
fAODevent(NULL),
fPIDResponse(NULL),
fAODeventCuts(),
fUtils(NULL),
fOutputList(NULL),
fQAList(NULL),
fCentralityMin(0),
fCentralityMax(0),
fVertexZmin(0),
fVertexZmax(0),
fNumberOfVertexContributorsMin(0),
fMagFieldSign(0),
fDCAxyMax(0),
fDCAzMax(0),
fPtMin(0),
fPtMax(0),
fEtaMax(0),
fNumberOfClustersITSMin(0),
fNumberOfClustersTPCMin(0),
fNumberOfClustersTRDMin(0),
fNumberOfClustersTPCdEdxMin(0),
fNumberOfCrossedRowsMin(0),
fCrossedRowsOverFindableClsMin(0),
fChiSquarePerNDFMax(0),
fChiSquareTRDMax(0),
fNumberOfTRDtrackletsPIDMin(0),
fIsTRDrefit(0),
fTPCnsigmaDeuteronMax(0),
fTOFnsigmaDeuteronMax(0)
{
    fUtils = new AliAnalysisUtils();
    DefineInput (0, TChain::Class());
    DefineOutput(1, TList::Class());
    DefineOutput(2, TList::Class());
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
AliAnalysisTaskJulia::~AliAnalysisTaskJulia()  {
    
    fOutputList->Clear();
    delete fAODevent;
    delete fPIDResponse;
    delete fUtils;
    delete fOutputList;
    delete fQAList;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskJulia::UserCreateOutputObjects()  {
    
    fOutputList = new TList();
    fOutputList -> SetOwner();
    
    fQAList = new TList();
    fQAList -> SetOwner();
    
    //QA Plots of Event Selection
    fAODeventCuts.AddQAplotsToList(fQAList);
    
    
    //Number of Events
    hEvents = new TH1F ("hEvents","",10,0,10);
    fOutputList -> Add(hEvents);
    
    
    //dE/dx vs. momentum
    hdEdx_vs_momentum_positive = new TH2F ("hdEdx_vs_momentum_positive","",1000,0.1,10.0,500,0,1000.0);
    hdEdx_vs_momentum_negative = new TH2F ("hdEdx_vs_momentum_negative","",1000,0.1,10.0,500,0,1000.0);
    fOutputList -> Add (hdEdx_vs_momentum_positive);
    fOutputList -> Add (hdEdx_vs_momentum_negative);

    
    
    PostData(1, fOutputList);
    PostData(2, fQAList);
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskJulia::UserExec(Option_t *)  {
    

    //Get Input Event
    if ( !GetEvent ()) return;
    

    //Load PID Response
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler *inputHandler = (AliInputEventHandler*) (mgr->GetInputEventHandler());
    fPIDResponse = inputHandler->GetPIDResponse();
    
    
    //Loop over Reconstructed Tracks
    for (Int_t i=0 ; i<fAODevent->GetNumberOfTracks() ; i++)  {
        
        //Track Selection
        AliAODTrack *track = (AliAODTrack*) fAODevent -> GetTrack(i);
        if ( !track ) continue;
        if ( !IsGoodQualityTrack (track)) continue;
        
        Double_t dEdx = track -> GetTPCsignal();
        Double_t p    = track -> GetTPCmomentum();
        
        if (track -> Charge() > 0) hdEdx_vs_momentum_positive -> Fill (p,dEdx);
        if (track -> Charge() < 0) hdEdx_vs_momentum_negative -> Fill (p,dEdx);
    }
    
    
    PostData(1, fOutputList);
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskJulia::GetEvent ()  {
    
    //Get Input Event
    fAODevent = dynamic_cast <AliAODEvent*>(InputEvent());
    if (!fAODevent) return false;
    hEvents -> Fill(0.5);
    
    //Standard Event Cuts
    if (!fAODeventCuts.AcceptEvent(fAODevent)) {
        PostData(2, fQAList);
        return false;
    }
    hEvents -> Fill(1.5);

    //Centrality
    AliMultSelection *multiplicitySelection = (AliMultSelection*) fAODevent->FindListObject("MultSelection");
    if( !multiplicitySelection) return false;
    hEvents -> Fill(2.5);
    Double_t centrality = multiplicitySelection->GetMultiplicityPercentile("V0M");
    
    //Selection of Centrality Range
    if (centrality<fCentralityMin || centrality>=fCentralityMax ) return false;
    hEvents -> Fill(3.5);

    //Primary Vertex
    AliAODVertex *vertex = (AliAODVertex*) fAODevent->GetPrimaryVertex();
    if ( !vertex ) return false;
    hEvents -> Fill(4.5);
    
    //Primary Vertex Selection
    if ( vertex->GetZ() < fVertexZmin ) return false;
    if ( vertex->GetZ() > fVertexZmax ) return false;
    hEvents -> Fill(5.5);
    
    if ( vertex->GetNContributors() < fNumberOfVertexContributorsMin ) return false;
    hEvents -> Fill(6.5);

    //Event Plane
    AliEventplane *eventPlane =  fAODevent->GetEventplane();
    if (!eventPlane) return false;
    hEvents -> Fill(7.5);

    //Magnetic Field Orientation
    Double_t signMagneticField = fAODevent->GetMagneticField()/TMath::Abs(fAODevent->GetMagneticField());
    if (fMagFieldSign>0 && signMagneticField<0) return false;
    if (fMagFieldSign<0 && signMagneticField>0) return false;
    hEvents -> Fill(8.5);
    
    return true;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskJulia::IsGoodQualityTrack (AliAODTrack *track)  {
    
    //Filterbit
    if(!track->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) return false;

    
    //Kinematic Cuts & Acceptance
    if ( track->Pt()<fPtMin || track->Pt()>fPtMax ) return false;
    if ( TMath::Abs(track->Eta()) > fEtaMax )       return false;
    
    //Track Selection Cuts
    if ( track->GetITSNcls() < fNumberOfClustersITSMin ) return false;
    if ( track->GetTPCNcls() < fNumberOfClustersTPCMin ) return false;
    if ( track->GetTRDncls() < fNumberOfClustersTRDMin ) return false;
    if ( track->GetTPCNCrossedRows() < fNumberOfCrossedRowsMin ) return false;
    if ( static_cast<Double_t>(track->GetTPCNCrossedRows())/static_cast<Double_t>(track->GetTPCNclsF()) < fCrossedRowsOverFindableClsMin) return false;
    if ( track->GetTPCsignalN() < fNumberOfClustersTPCdEdxMin ) return false;
    if ( track->Chi2perNDF() > fChiSquarePerNDFMax) return false;
    
    
    //DCA
    Double_t DCAz  = GetDCAz  (track);
    Double_t DCAxy = GetDCAxy (track);
    
    if (TMath::Abs(DCAz)  > fDCAzMax)  return false;
    if (TMath::Abs(DCAxy) > fDCAxyMax) return false;

    return true;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskJulia::GetDCAxy (AliAODTrack *track)  {
    
    Double_t impactParameter[2];
    Double_t covarianceMatrix[3];
    if (!track->PropagateToDCA (fAODevent->GetPrimaryVertex(),fAODevent->GetMagneticField(),10000,impactParameter,covarianceMatrix)) return -999;
    
    Double_t DCAxy = impactParameter[0];
    
    return DCAxy;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskJulia::GetDCAz (AliAODTrack *track)  {
    
    Double_t impactParameter[2];
    Double_t covarianceMatrix[3];
    if (!track->PropagateToDCA (fAODevent->GetPrimaryVertex(),fAODevent->GetMagneticField(),10000,impactParameter,covarianceMatrix)) return -999;
    
    Double_t DCAz = impactParameter[1];
    
    return DCAz;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskJulia::IsDeuteronCandidate (AliAODTrack *track)  {
    
    Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC (track,AliPID::kDeuteron);
    if (TMath::Abs(nsigmaTPC) > fTPCnsigmaDeuteronMax) return false;

    return true;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskJulia::Terminate(Option_t *)  {
    
    fOutputList = dynamic_cast<TList*> (GetOutputData(1));
    if (!fOutputList) return;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________

