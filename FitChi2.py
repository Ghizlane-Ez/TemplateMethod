import ROOT 
import ROOT as root
import atlasplots as aplt 
from ROOT import gROOT, TCanvas, TFile, THStack, TH1F, TPad, TLine, TAttFill, TF1, TGraph, gROOT, gRandom

def calChi2(data, MC_tem):
    print("chi2 Calculation")
    print(data.Chi2Test(MC_tem,"CHI2WW"))

def fitChi2(var,labels,limy1,limy2):
    # Set the ATLAS Style
    aplt.set_atlas_style()

    # Create a figure and axis
    fig, ax = aplt.subplots(name="fig1", figsize=(800, 800))
    ax.set_ylim(0, 800)

    # Open the data file
    data_file = ROOT.TFile.Open("data17_WZ_lowMu_repro_5TeV.root")
    data_hist = data_file.Get("WplusenuSelection/"+var+"_cut8")
    
    # Open the MC files
    MC_files = ROOT.TFile.Open("mc16_5TeV.361100.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusenu.e4916_s3238_r10243_r10210_p3665.root")
    MC_Nominal = MC_files.Get("WplusenuSelection/"+var+"_cut8")

    # Define a list of histogram names
    hist_names = ["WplusenuSelection_WeightVariations/"+var+"_cut8_dM-200_dG000", "WplusenuSelection_WeightVariations/"+var+"_cut8_dM-150_dG000", "WplusenuSelection_WeightVariations/"+var+"_cut8_dM-100_dG000","WplusenuSelection_WeightVariations/"+var+"_cut8_dM-50_dG000","WplusenuSelection_WeightVariations/"+var+"_cut8_dM-25_dG000","WplusenuSelection/"+var+"_cut8","WplusenuSelection_WeightVariations/"+var+"_cut8_dM025_dG000","WplusenuSelection_WeightVariations/"+var+"_cut8_dM050_dG000","WplusenuSelection_WeightVariations/"+var+"_cut8_dM100_dG000","WplusenuSelection_WeightVariations/"+var+"_cut8_dM150_dG000","WplusenuSelection_WeightVariations/"+var+"_cut8_dM200_dG000"]
    
    #define a list of chi2 values
    chi2_list1 = []
    chi2_list2 = []

    # Loop over the histograms
    for i in range(len(hist_names)):
        # get the histogram
        MC_hist = MC_files.Get(hist_names[i])
        if i != 5:
            MC_hist.Add(MC_Nominal)
        # Calculate the chi2
        print(i, MC_hist)
        #calChi2(MC_Nominal, MC_hist)
        chi2_list1.append(MC_Nominal.Chi2Test(MC_hist,"CHI2WW"))
        chi2_list2.append(data_hist.Chi2Test(MC_hist,"CHI2WW"))

    # Close the files
    data_file.Close()
    
    #define mass step
    MassSteps=[-0.2,-0.15,-0.1,-0.05,-0.025,0,0.025,0.05,0.1,0.15,0.2]

    #define the graphs
    graph1 = TGraph(len(chi2_list1)) 
    for i in range(0, len(chi2_list1)):
         for i in range(0, len(MassSteps)):
            graph1.SetPoint(i, 80.399+MassSteps[i], chi2_list1[i])
    graph2 = TGraph(len(chi2_list2)) 
    for i in range(0, len(chi2_list2)):
         for i in range(0, len(MassSteps)):
            graph2.SetPoint(i, 80.399+MassSteps[i], chi2_list2[i])        

     # Fit function
    fitchi21  = ROOT.TF1("fitchi21","(x-[0])*(x-[0])/[1]/[1] + [2]",80.150,80.650)
    fitchi21.SetParameter(0,80.399)
    fitchi21.SetParameter(1,10)
    fitchi21.SetParameter(2,1.)
    #fit graph 1
    graph1.Fit("fitchi21","rQ")
    graph1.Fit("fitchi21","rQ")
    graph1.Fit("fitchi21","rQ")
    # print the results.
    print("Nominal : ",   fitchi21.GetParameter(0))
    print("Stat error : ",fitchi21.GetParameter(1))

    # Fit function
    fitchi22  = ROOT.TF1("fitchi22","(x-[0])*(x-[0])/[1]/[1] + [2]",80.150,80.650)
    fitchi22.SetParameter(0,80.399)
    fitchi22.SetParameter(1,10)
    fitchi22.SetParameter(2,1.)
    #fit graph 2
    graph2.Fit("fitchi22","rQ")
    graph2.Fit("fitchi22","rQ")
    graph2.Fit("fitchi22","rQ")
    # print the results.
    print("Nominal : ",   fitchi22.GetParameter(0))
    print("Stat error : ",fitchi22.GetParameter(1))

    results1 = []
    results1.append(fitchi21.GetParameter(0))
    results1.append(fitchi21.GetParameter(1))
    print(results1)
    results2 = []
    results2.append(fitchi22.GetParameter(0))
    results2.append(fitchi22.GetParameter(1))
    print(results2)
    c1 = TCanvas("c1","c1", 1000, 800)
   
    # Define x an y titles
    graph1.GetXaxis().SetTitle(labels) # p_{T}^{e}[GeV]   m_{T}^{W} [GeV]
    graph1.GetXaxis().SetTitleOffset(1.6)
    graph1.GetYaxis().SetTitle("\chi^{2}")
    graph1.GetYaxis().SetTitleOffset(1.2)
    graph1.GetYaxis().SetTitleSize(0.05)
    graph1.GetYaxis().SetTitleFont(42)

    #color of fit 
    fitchi21.SetLineColor(ROOT.kBlue)
    fitchi22.SetLineColor(ROOT.kRed)

    yaxis = graph1.GetYaxis()

    yaxis.SetRangeUser(limy1,limy2)

    #draw graphs and fit functions in the same canvas
    graph1.Draw("AL*")
    fitchi21.Draw("same")
    graph2.Draw("L*same")
    fitchi22.Draw("same")

    # Set left margin
    c1.SetLeftMargin(0.15)

    # Add titles
    title1 = ROOT.TLatex()
    title1.SetTextSize(0.04)
    title1.SetTextAlign(11)
    title1.SetTextFont(42)
    title1.DrawLatexNDC(0.22, 0.89, "ATLAS Internal")
    title2 = ROOT.TLatex()
    title2.SetTextSize(0.036)
    title2.SetTextAlign(11)
    title2.SetTextFont(42)
    title2.DrawLatexNDC(0.23, 0.84, "#sqrt{s} = 5 TeV, 256.827 pb^{-1}")

    # draw two lines at top right
    line_mc = ROOT.TLine(0.9, 0.8, 0.9, 0.8)
    line_mc.SetLineColor(ROOT.kBlue)
    line_mc.Draw()

    line_data = ROOT.TLine(0.8, 0.75, 0.9, 0.75)
    line_data.SetLineColor(ROOT.kRed)
    line_data.Draw()

    # add a legend
    legend = ROOT.TLegend(0.72, 0.8, 0.8, 0.9)
    legend.AddEntry(line_mc, "MC", "l")
    legend.AddEntry(line_data, "DATA", "l")
    legend.Draw()

    #print c1
    c1.Print("Outputs/chi2FitMC-DATA-"+var+".png")

