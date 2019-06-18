#include "TFile.h"
#include "TList.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"




//_______________________________________________________________________________________________________________________________________________________________________
void data_analysis ()  {
    
    
    
    //Get Input
    TList *input = GetInput ();
    
    
    //Get Input Histograms
    TH2F *hdEdx_vs_momentum_positive = (TH2F*) input -> FindObject ("hdEdx_vs_momentum_positive");
    hdEdx_vs_momentum_positive -> SetStats(false);
    hdEdx_vs_momentum_positive -> GetXaxis() -> SetTitle("#it{p} (GeV/#it{c})");
    hdEdx_vs_momentum_positive -> GetXaxis() -> SetTitleOffset(1.6);
    hdEdx_vs_momentum_positive -> SetTitleSize(0.045,"x");
    hdEdx_vs_momentum_positive -> SetLabelSize(0.045,"x");
    hdEdx_vs_momentum_positive -> GetXaxis() -> CenterTitle();
    hdEdx_vs_momentum_positive -> GetYaxis() -> SetTitle("TPC d#it{E} / d#it{x} (a.u.)");
    hdEdx_vs_momentum_positive -> GetYaxis() -> SetTitleOffset(1.8);
    hdEdx_vs_momentum_positive -> SetTitleSize(0.045,"y");
    hdEdx_vs_momentum_positive -> SetLabelSize(0.045,"y");
    hdEdx_vs_momentum_positive -> GetYaxis() -> CenterTitle();

    
    //Plot dE/dx for positive tracks
    TCanvas *cdEdx_vs_momentum_positive = new TCanvas ("cdEdx_vs_momentum_positive","",900,700);
    cdEdx_vs_momentum_positive -> cd();
    cdEdx_vs_momentum_positive -> SetTickx(1);
    cdEdx_vs_momentum_positive -> SetTicky(1);
    cdEdx_vs_momentum_positive -> SetLeftMargin(0.18);
    cdEdx_vs_momentum_positive -> SetRightMargin(0.1);
    cdEdx_vs_momentum_positive -> SetBottomMargin(0.15);
    cdEdx_vs_momentum_positive -> SetTopMargin(0.1);
    hdEdx_vs_momentum_positive -> Draw("colz");
    Text() -> Draw();
}
//_______________________________________________________________________________________________________________________________________________________________________
TList *GetInput ()  {
    
    TFile *inputfile = TFile::Open ("input/InputFile.root");
    TList *inputlist = (TList*) inputfile -> Get ("Input");
    return inputlist;
}
//_______________________________________________________________________________________________________________________________________________________________________
TPaveText *Text ()  {
    
    TPaveText *text = new TPaveText(0.35,0.75,0.7,0.88, "NDC");
    text -> SetTextSize(0.045);
    text -> SetTextFont(42);
    text -> SetTextColor(1);
    text -> SetFillColor(0);
    text -> SetBorderSize(0);
    text -> AddText ("ALICE work in progress");
    text -> AddText ("Pb#font[122]{-}Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV");

    return text;
}
//_______________________________________________________________________________________________________________________________________________________________________




