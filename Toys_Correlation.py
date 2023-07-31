import ROOT 
import ROOT as root
import atlasplots as aplt 
from ROOT import gROOT, TCanvas, TFile, THStack, TH1F, TPad, TLine, TAttFill, TF1, TGraph, gROOT, gRandom
from numpy import c_
import numpy as np
from array import array
import matplotlib.pyplot as plt
from  math import *
import copy

################################################################################################
 
def get_bootstrap_hist(file_name, hist_name_format, n_toys, nominal_hist):

    bootstrap_hists = []
    for i in range(0, n_toys):
        hist = file_name.Get(hist_name_format.format(i)).Clone()
        hist.Add(nominal_hist.Clone())
        bootstrap_hists.append(hist.Clone())

    return bootstrap_hists

##############################################################################################

def calculate_chi2_for_bootstrap_hists(bootstrap_hists, filename, hist_names,nominal_hist):
    chi2_list = [[] for i in range(len(bootstrap_hists))]
    for i in range(len(bootstrap_hists)):
        for j in range(len(hist_names)):
            MC_hist = filename.Get(hist_names[j]).Clone()
            MC_hist.Add(nominal_hist.Clone()) 
            chi2_list[i].append(bootstrap_hists[i].Chi2Test(MC_hist, "CHI2WW"))
    return chi2_list
    
#################################################################################################

def create_graphs(chi2_list):

    MassSteps=[-0.2,-0.15,-0.1,-0.05,-0.025,0,0.025,0.05,0.1,0.15,0.2]
    graphs = []
    for i in range(0, len(chi2_list)):
        grph = ROOT.TGraph(len(MassSteps))
        for j in range(len(MassSteps)):
            grph.SetPoint(j, 80.399+ MassSteps[j], chi2_list[i][j])
        graphs.append(grph)
    return graphs

###################################################################################################

# Define the fitting function
def create_fitchi2():
    fitchi2 = ROOT.TF1("fitchi2","(x-[0])*(x-[0])/[1]/[1] + [2] +[3]*x*x*x",80.150,80.650) #(x-[0])*(x-[0])/[1]/[1] + [2]
    fitchi2.SetParameter(0,80.399)
    fitchi2.SetParameter(1,10)
    fitchi2.SetParameter(2,1.)
    fitchi2.SetParameter(3, 1.)
    return fitchi2

##################################################################################################

# Fit the function to each TGraph and return a list of fitted parameter 0
def fit_graphs(graphs, fitchi2):
    param0_list = []
    for i in range(len(graphs)):
        graphs[i].Fit(fitchi2)
        param0_list.append(fitchi2.GetParameter(0))
    return param0_list          

##########################################################################################"""

def GraphPloting(listelPt, listemT):

        n = len(listelPt)
        x, y = array( 'd' ), array( 'd' )
        for i in range( n ):
             x.append( listelPt[i] )
             y.append( listemT[i]  )

        gr = TGraph(n, x, y)
        axis = gr.GetXaxis()
        axis.SetLimits(80.370,80.430)
        gr.GetHistogram().GetYaxis().SetRangeUser(80.360, 80.430)  # along X

        c1 = TCanvas("c1", "Graph", 900, 700)
        gr.Draw("APL")
        gr.SetLineWidth(0)
        gr.SetMarkerStyle(20)
        gr.GetXaxis().SetTitle("p_{T}^{e}")
        gr.GetYaxis().SetTitle("m_{T}^{W}")
        gr.GetXaxis().CenterTitle()
        gr.GetYaxis().CenterTitle()
        gr.GetYaxis().SetTitleOffset(2.1)

        # Define left margin
        c1.SetLeftMargin(0.15)

        # Add titles
        title1 = ROOT.TLatex()
        title1.SetTextSize(0.04)
        title1.SetTextAlign(11)
        title1.SetTextFont(42)
        title1.DrawLatexNDC(0.22, 0.89, "ATLAS Internal")

        c1.SaveAs("corr.pdf")
        return gr

###################################################################################################################

def CorrelationCoefficient(pTlMwToys, mTMwToys):
 
    MwpTlmoy = 0
    MwmTmoy  = 0
    for i in range(0, len(pTlMwToys)):
            MwpTlmoy = MwpTlmoy + pTlMwToys[i]
            MwmTmoy  = MwmTmoy  + mTMwToys[i]

    mTErrorStat   = 0
    elPtErrorStat = 0
    for i in range(0, len(mTMwToys)):
            mTErrorStat   = mTErrorStat   + (mTMwToys[i]  - (MwmTmoy/len(mTMwToys)))*(mTMwToys[i] - (MwmTmoy/len(mTMwToys)))
            elPtErrorStat = elPtErrorStat + (pTlMwToys[i] - (MwpTlmoy/len(mTMwToys)))*(pTlMwToys[i] - (MwpTlmoy/len(mTMwToys)))

    Cor = 0
    for i in range(0, len(mTMwToys)):
            Cor = Cor + (mTMwToys[i]  - (MwmTmoy/len(mTMwToys)))*(pTlMwToys[i] - (MwpTlmoy/len(mTMwToys)))

    print("cor :", Cor/sqrt(mTErrorStat*elPtErrorStat))

    return Cor        

#################################################################################################################################

def main():
    # Set the ATLAS Style
    aplt.set_atlas_style()
    # Create a figure and axis
    fig, ax = aplt.subplots(name="fig1", figsize=(800, 800))
    ax.set_ylim(0, 800)
    
    # Open the MC files
    MC_files = ROOT.TFile.Open("mc16_5TeV.361100.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusenu.e4916_s3238_r10243_r10210_p3665.root")
    # Open the data file
    DATA_files = ROOT.TFile.Open("data17_WZ_lowMu_repro_5TeV.root")
    eta=1    #2,3,4
    # get elPt nominale MC
    MC_Nominal1 = MC_files.Get("WplusenuSelection/elPt_eta{}_cut8".format(eta))
    # get mT nominale MC
    MC_Nominal2 = MC_files.Get("WplusenuSelection/mT_eta{}_cut8".format(eta))  #_eta2
    # get elPt nominale DATA
    DATA_Nominal1 = DATA_files.Get("WplusenuSelection/elPt_eta{}_cut8".format(eta))
    # get mT nominale DATA
    DATA_Nominal2 = DATA_files.Get("WplusenuSelection/mT_eta{}_cut8".format(eta))  #_eta2                                                                                                                                                                                                                                                                                                      elPt_cut8_dM000_dG000 "WplusenuSelection/elPt_eta{}_cut8_dM000_dG000".format(eta)  "WplusenuSelection/mT_eta{}_cut8_dM000_dG000".format(eta)
    #hist_names elPt
    hist_names1 = ["WplusenuSelection_WeightVariations/elPt_eta{}_cut8_dM-200_dG000".format(eta), "WplusenuSelection_WeightVariations/elPt_eta{}_cut8_dM-150_dG000".format(eta), "WplusenuSelection_WeightVariations/elPt_eta{}_cut8_dM-100_dG000".format(eta),"WplusenuSelection_WeightVariations/elPt_eta{}_cut8_dM-50_dG000".format(eta),"WplusenuSelection_WeightVariations/elPt_eta{}_cut8_dM-25_dG000".format(eta),"WplusenuSelection_WeightVariations/elPt_eta{}_cut8_dM000_dG000".format(eta),"WplusenuSelection_WeightVariations/elPt_eta{}_cut8_dM025_dG000".format(eta),"WplusenuSelection_WeightVariations/elPt_eta{}_cut8_dM050_dG000".format(eta),"WplusenuSelection_WeightVariations/elPt_eta{}_cut8_dM100_dG000".format(eta),"WplusenuSelection_WeightVariations/elPt_eta{}_cut8_dM150_dG000".format(eta),"WplusenuSelection_WeightVariations/elPt_eta{}_cut8_dM200_dG000".format(eta)]

    #hist_names mT
    hist_names2 = ["WplusenuSelection_WeightVariations/mT_eta{}_cut8_dM-200_dG000".format(eta), "WplusenuSelection_WeightVariations/mT_eta{}_cut8_dM-150_dG000".format(eta), "WplusenuSelection_WeightVariations/mT_eta{}_cut8_dM-100_dG000".format(eta),"WplusenuSelection_WeightVariations/mT_eta{}_cut8_dM-50_dG000".format(eta),"WplusenuSelection_WeightVariations/mT_eta{}_cut8_dM-25_dG000".format(eta),"WplusenuSelection_WeightVariations/mT_eta{}_cut8_dM000_dG000".format(eta),"WplusenuSelection_WeightVariations/mT_eta{}_cut8_dM025_dG000".format(eta),"WplusenuSelection_WeightVariations/mT_eta{}_cut8_dM050_dG000".format(eta),"WplusenuSelection_WeightVariations/mT_eta{}_cut8_dM100_dG000".format(eta),"WplusenuSelection_WeightVariations/mT_eta{}_cut8_dM150_dG000".format(eta),"WplusenuSelection_WeightVariations/mT_eta{}_cut8_dM200_dG000".format(eta)]
     
    #Open the bootstrap files
    bootstrap_file = ROOT.TFile.Open("DataForM2/MC_Bootstrap/mc16_5TeV.361100.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusenu.e4916_s3238_r10243_r10210_p3665.root")

    #appl func1 get toys and add mC nominale to toys
    bootstrap_hist1 = get_bootstrap_hist( bootstrap_file, "WplusenuSelection_WeightVariations/elPt_eta{}_cut8_toy{{0}}".format(eta), 400, MC_Nominal1) 
    bootstrap_hist2 = get_bootstrap_hist( bootstrap_file, "WplusenuSelection_WeightVariations/mT_eta{}_cut8_toy{{0}}".format(eta), 400,   MC_Nominal2 )
   
    #appl fun2 chi2 calcul   
    chi2_list1 = calculate_chi2_for_bootstrap_hists( bootstrap_hist1, MC_files,  hist_names1, MC_Nominal1) 
    chi2_list2 = calculate_chi2_for_bootstrap_hists( bootstrap_hist2, MC_files,  hist_names2, MC_Nominal2) 
    
    for i in range(0, 400):             
         print("Chi2 list for Bootstrap elPt {0}: {1}".format(i, chi2_list1[i]))
         print("Chi2 list for Bootstrap mT {0}: {1}".format(i, chi2_list2[i]))

    #func 3
    graphs1 =create_graphs(chi2_list1)
    graphs2 =create_graphs(chi2_list2)
    fitchi2 =create_fitchi2()
    # Fit the function to each TGraph and get a list of fitted parameter 0 for elPt
    param0_listelPt = fit_graphs(graphs1, fitchi2)
    print("mW from elPt:")
    print(param0_listelPt)
    print(len(param0_listelPt))

    # Fit the function to each TGraph and get a list of fitted parameter 0 for mT
    param0_listemT = fit_graphs(graphs2, fitchi2)
    
    print("mW from mT :")
    print(param0_listemT)
    print(len(param0_listemT))

    #ghraph plot
    GraphPloting(param0_listelPt, param0_listemT)

    # correlation calcul
    CorrelationCoefficient(param0_listelPt, param0_listemT)      
   
if __name__ == '__main__':
    ROOT.gROOT.SetBatch()
    main()

