import ROOT 
import ROOT as root
import atlasplots as aplt

def truthreco(var1,var,mid,xlim1, xlim2,ylim1, ylim2,labele):
    # Set the ATLAS Style
    aplt.set_atlas_style()

    # Create a figure and axis
    fig, ax = aplt.subplots(name="fig1", figsize=(800, 800))

    # Define distributions
    MC_Input1   = ROOT.TFile.Open("mc16_5TeV.361100.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusenu.e4916_s3238_r10243_r10210_p3665.root")
    hist1 = MC_Input1.Get("TruthSelection/truth"+var1+"_cut4")
  
    MC_Input2   = root.TFile.Open("mc16_5TeV.361100.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusenu.e4916_s3238_r10243_r10210_p3665.root")
    hist2 = MC_Input2.Get("WplusenuSelection/"+var+"_cut8")

    # Find the peak height of the truth histogram
    maxbin = hist1.GetMaximum()
    binnum = mid 
    
    # Draw a vertical line to the peak of the truth histogram
    line = root.TLine(binnum, 0, binnum, maxbin)
    line.SetLineWidth(2)
    line.SetLineStyle(2)
    line.SetLineColor(2)
    ax.plot(line)

    # Draw the histograms on the same axis
    ax.plot(hist1, linecolor='green', linewidth=2, label="Truth", labelfmt="L")
    ax.plot(hist2, linecolor=root.kBlue+1, label="Reconstructed", labelfmt="L")
    ax.set_xlim(xlim1, xlim2)
    ax.set_ylim(ylim1, ylim2)

    # Add extra space at top of plot to make room for labels
    ax.add_margins(top=0.16)

    # Set axis titles
    ax.set_ylabel("Events ", titleoffset=2.2)
    ax.set_xlabel(labele, titleoffset=1.4)  #m_{T}^{W}[GEV]  p_{T}^{l}

    # Go back to top axes to add labels
    ax.cd()

    # Add the ATLAS Label
    aplt.atlas_label(text="Internal", loc="upper left")
    ax.text(0.2, 0.84, "#sqrt{s} = 5 TeV, 256.827 pb^{-1}", size=22, align=13)

    # Add legend
    ax.legend(loc=(0.78, 0.78, 0.6, 0.90))

    # Save the plot as a Png
    fig.savefig("Outputs/truth-recon-MC-"+var+".png")


