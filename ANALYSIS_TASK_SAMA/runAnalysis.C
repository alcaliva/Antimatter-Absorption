#define ALIPHYSICS_VER  "vAN-20190301-1"
#define GRIDWorkingDir  "SAMA_DIR"
#define AnalysisMacro   "AnalysisSama"
#define AnalysisTask    "AliAnalysisTaskSama"
#define TTL             20000
#define nRunsPerMaster  1


//__________________________________________________________________________________________________________________________________________________________________________
void LoadAnalysisTask (AliAnalysisManager *mgr)  {
    
    
    //Analysis Cuts
    Double_t CentralityMin =  0.0;
    Double_t CentralityMax = 90.0;
    Double_t VertexZmin = -10.0;
    Double_t VertexZmax = +10.0;
    Int_t    NumberOfVertexContributorsMin = 1;
    Double_t MagFieldSign = 0.0;
    Double_t PtMin =  0.5;
    Double_t PtMax = 50.0;
    Double_t EtaMax = 0.8;
    Double_t DCAxyMax = 0.1;
    Double_t DCAzMax  = 1.0;
    Int_t    NumberOfClustersITSMin = 2;
    Int_t    NumberOfClustersTPCMin = 80;
    Int_t    NumberOfClustersTRDMin = 2;
    Int_t    NumberOfClustersTPCdEdxMin = 50;
    Int_t    NumberOfCrossedRowsMin = 70;
    Double_t CrossedRowsOverFindableClsMin = 0.8;
    Double_t ChiSquarePerNDFMax = 4.0;
    Double_t ChiSquareTRDMax = 20.0;
    Int_t    NumberOfTRDtrackletsPIDMin = 0;
    Bool_t   IsTRDrefit = (false);
    Double_t TPCnsigmaDeuteronMax = 3.0;
    Double_t TOFnsigmaDeuteronMax = 3.0;

    
    //Load Analysis Tasks
    gROOT->LoadMacro("AliAnalysisTaskSama.cxx+g");
    
    //Input container
    AliAnalysisDataContainer *input = mgr -> GetCommonInputContainer();

    //Analysis Task
    AliAnalysisTaskSama *task = new AliAnalysisTaskSama  ("AliAnalysisTaskSama");
    task -> SelectCollisionCandidates (AliVEvent::kINT7);
    task -> AliAnalysisTaskSama::SetAnalysisCuts (CentralityMin,CentralityMax,VertexZmin,VertexZmax,NumberOfVertexContributorsMin,MagFieldSign,DCAxyMax,DCAzMax,PtMin,PtMax,EtaMax,NumberOfClustersITSMin,NumberOfClustersTPCMin,NumberOfClustersTRDMin,NumberOfClustersTPCdEdxMin,NumberOfCrossedRowsMin,CrossedRowsOverFindableClsMin,ChiSquarePerNDFMax,ChiSquareTRDMax,NumberOfTRDtrackletsPIDMin,IsTRDrefit,TPCnsigmaDeuteronMax,TOFnsigmaDeuteronMax);
    mgr -> AddTask(task);
    AliAnalysisDataContainer *output   = mgr -> CreateContainer("Input",TList::Class(),AliAnalysisManager::kOutputContainer,"InputFile.root");
    AliAnalysisDataContainer *outputQA = mgr -> CreateContainer("InputQA",TList::Class(),AliAnalysisManager::kOutputContainer,"Input_QA.root");
    mgr -> ConnectInput (task,0,input);
    mgr -> ConnectOutput(task,1,output);
    mgr -> ConnectOutput(task,2,outputQA);
    
}
//__________________________________________________________________________________________________________________________________________________________________________
void runAnalysis (const char *mode="full", Bool_t merge=true)  {
    
    //Load Libraries
    LoadLibraries();
    
    //Grid Connection
    TGrid::Connect("alien://");

    //Alien Handler
    AliAnalysisGrid *alienHandler = CreateAlienHandler (mode,merge);
    if (!alienHandler) return;
    
    //Analysis Manager
    AliAnalysisManager *mgr = new AliAnalysisManager("AnalysisManager");
    mgr->SetGridHandler(alienHandler);

    //Event Handler, PID & Centrality
    EventHandler(mgr);
    LoadPhysicsSelection ();
    LoadCentrality ();
    LoadEPSelection();
    LoadPIDResponse();
    
    //Analysis Task
    LoadAnalysisTask(mgr);
    
    //Start analysis
    if (!mgr->InitAnalysis())  return;
    mgr->PrintStatus();
    mgr->StartAnalysis("grid");
};
//__________________________________________________________________________________________________________________________________________________________________________
AliAnalysisGrid* CreateAlienHandler (const char *mode, Bool_t merge )  {
    
    //Alien Handler
    AliAnalysisAlien *alien = new AliAnalysisAlien();
    SetAdditionalLibraries (alien);
    alien->SetOverwriteMode();
    alien->SetCheckCopy(false);
    alien->SetRunMode(mode);
    alien->SetNtestFiles(10);
    alien->SetAPIVersion("V1.1x");
    alien->SetAliPhysicsVersion(ALIPHYSICS_VER);
    alien->AddIncludePath("$ALICE_PHYSICS/include");
    alien->AddIncludePath("$ALICE_ROOT/include");
    alien->SetGridDataDir("/alice/data/2015/LHC15o");
    alien->SetDataPattern("*/pass1_pidfix/AOD/*/AliAOD.root");
    alien->SetRunPrefix("000");
    SetInputRuns (alien,mode);
    alien->SetNrunsPerMaster(nRunsPerMaster);
    alien->SetGridWorkingDir(Form("%s",GRIDWorkingDir));
    alien->SetGridOutputDir("OUTPUT");
    alien->SetAnalysisSource(Form("%s.cxx",AnalysisTask));
    alien->SetAdditionalLibs(Form("%s.cxx  %s.h",AnalysisTask,AnalysisTask));
    alien->SetMergeViaJDL(merge);
    alien->SetMaxMergeStages(2);
    alien->SetAnalysisMacro(Form("%s.C",AnalysisMacro));
    alien->SetSplitMaxInputFileNumber(50);
    alien->SetMasterResubmitThreshold(90);
    alien->SetTTL(TTL);
    alien->SetExecutable(Form("%s.sh",AnalysisMacro));
    alien->SetInputFormat("xml-single");
    alien->SetJDLName(Form("%s.jdl",AnalysisMacro));
    alien->SetPrice(1);
    alien->SetSplitMode("se");
    return alien;
}
//__________________________________________________________________________________________________________________________________________________________________________
void EventHandler (AliAnalysisManager *mgr)  {
    
    AliAODInputHandler* aodH = new AliAODInputHandler();
    mgr->SetInputEventHandler(aodH);
}
//__________________________________________________________________________________________________________________________________________________________________________
void LoadPhysicsSelection()  {
    
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
    AliPhysicsSelectionTask *physSel = AddTaskPhysicsSelection(false, true);
}
//__________________________________________________________________________________________________________________________________________________________________________
void LoadPIDResponse ()  {
    
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
    AliAnalysisTaskPIDResponse *taskPID = AddTaskPIDResponse(false);
}
//__________________________________________________________________________________________________________________________________________________________________________
void LoadEPSelection ()  {

    AliVZEROEPSelectionTask *vzeroEPSelectionTask;
    gROOT->ProcessLine(".L $ALICE_ROOT/ANALYSIS/macros/AddTaskVZEROEPSelection.C");
    vzeroEPSelectionTask = AddTaskVZEROEPSelection();
}
//__________________________________________________________________________________________________________________________________________________________________________
void LoadCentrality()  {
    
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
    AliMultSelectionTask *taskMult = AddTaskMultSelection(false);
    taskMult->SetSelectedTriggerClass(AliVEvent::kINT7);
}
//__________________________________________________________________________________________________________________________________________________________________________
void SetInputRuns (AliAnalysisAlien *alien, const char *mode)  {
    
    //Test Mode Option
    if (mode=="test") {alien->AddRunNumber(245145); return;}
    
    //Run List
    Int_t run[] = { 245145, 245146, 245151, 245152, 245231, 245232, 245259, 245343, 245345, 245346,
                    245347, 245349, 245353, 245396, 245397, 245401, 245407, 245409, 245411, 245439,
                    245441, 245446, 245450, 245452, 245454, 245496, 245497, 245501, 245504, 245505,
                    245507, 245535, 245540, 245542, 245543, 245544, 245545, 245554 };
    
    Int_t nRuns = sizeof(run)/sizeof(Int_t);
    nRuns = 1;
    for ( Int_t iRun=0 ; iRun<nRuns ; iRun++ )
        alien->AddRunNumber(run[iRun]);
        
}
//__________________________________________________________________________________________________________________________________________________________________________
void LoadLibraries()
{
    gSystem->Load("libCore.so");
    gSystem->Load("libTree.so");
    gSystem->Load("libGeom.so");
    gSystem->Load("libVMC.so");
    gSystem->Load("libPhysics.so");
    gSystem->Load("libMinuit.so");
    gSystem->Load("libGui.so");
    gSystem->Load("libXMLParser.so");
    gSystem->Load("libSTEERBase.so");
    gSystem->Load("libESD.so");
    gSystem->Load("libCDB.so");
    gSystem->Load("libAOD.so");
    gSystem->Load("libANALYSIS.so");
    gSystem->Load("libANALYSISalice.so");
    
    gROOT->ProcessLine(".include $ALICE_PHYSICS/include");
    gROOT->ProcessLine(".include $ALICE_ROOT/ANALYSIS/");
    gROOT->ProcessLine(".include $ALICE_ROOT/include");
    gROOT->ProcessLine(".include $ROOTSYS/include");
    
    gSystem->SetIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS -I$ALICE_PHYSICS/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/STEER/STEER -I$ALICE_ROOT/STEER/STEERBase -I$ALICE_ROOT/STEER/ESD -I$ALICE_ROOT/STEER/AOD -I$ALICE_ROOT/TRD -I$ALICE_ROOT/macros -I$ALICE_ROOT/ANALYSIS  -I$ALICE_PHYSICS/OADB -I$ALICE_PHYSICS/PWGHF -I$ALICE_PHYSICS/PWGHF/base -I$ALICE_PHYSICS/PWGHF/vertexingHF -I$ALICE_PHYSICS/PWG/FLOW/Base -I$ALICE_PHYSICS/PWG/FLOW/Tasks -I$ALICE_PHYSICS/JETAN -I$ALICE_PHYSICS/PWG/Tools -g");
}
//__________________________________________________________________________________________________________________________________________________________________________
void SetAdditionalLibraries (AliAnalysisAlien *alien)  {
    
    alien->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS -I$ALICE_PHYSICS/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/STEER/STEER -I$ALICE_ROOT/STEER/STEERBase -I$ALICE_ROOT/STEER/ESD -I$ALICE_ROOT/STEER/AOD -I$ALICE_ROOT/TRD -I$ALICE_ROOT/macros -I$ALICE_ROOT/ANALYSIS  -I$ALICE_PHYSICS/OADB -I$ALICE_PHYSICS/PWGHF -I$ALICE_PHYSICS/PWGHF/base -I$ALICE_PHYSICS/PWGHF/vertexingHF -I$ALICE_PHYSICS/PWG/FLOW/Base -I$ALICE_PHYSICS/PWG/FLOW/Tasks -I$ALICE_PHYSICS/JETAN -I$ALICE_PHYSICS/PWG/Tools -g");
    
    alien->SetAdditionalLibs("AliAnalysisTaskFlowITSTPCTOFQCSP_Template.cxx AliAnalysisTaskFlowITSTPCTOFQCSP_Template.h libPWGHFhfe.so libCDB.so libSTEER.so  libCORRFW.so libPWGflowBase.so libPWGflowTasks.so libGui.so libProof.so libMinuit.so libXMLParser.so libRAWDatabase.so libRAWDatarec.so libCDB.so libSTEERBase.so libSTEER.so libTPCbase.so libTOFbase.so libTOFrec.so libTRDbase.so libVZERObase.so libVZEROrec.so libT0base.so libT0rec.so ");
}
//__________________________________________________________________________________________________________________________________________________________________________
