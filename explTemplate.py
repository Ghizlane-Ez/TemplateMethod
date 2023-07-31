import ROOT 
import ROOT as root
import atlasplots as aplt 

def varnom(var,xmin, xmax,labele):

    # Define a distribution
    MC_Input1   = ROOT.TFile.Open("mc16_5TeV.361100.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Wplusenu.e4916_s3238_r10243_r10210_p3665.root")
    hist_Nomin  = MC_Input1.Get("WplusenuSelection/"+var+"_cut8")
    hist_delta1 = MC_Input1.Get("WplusenuSelection_WeightVariations/"+var+"_cut8_dM-200_dG000")
    hist_delta2 = MC_Input1.Get("WplusenuSelection_WeightVariations/"+var+"_cut8_dM200_dG000")    

    # Set the ATLAS Style
    aplt.set_atlas_style()

    # Create a figure and axes
    fig, (ax1, ax2) = aplt.ratio_plot(name="fig1", figsize=(800, 800), hspace=0.05)

    # Draw the histograms on these axes
    ax1.plot(hist_delta1, linecolor=2, label=""+var+" - 200 MeV", labelfmt="L")
    ax1.plot(hist_Nomin , linecolor=1, label=""+var+" ", labelfmt="L")
    ax1.plot(hist_delta2, linecolor=4, label=""+var+" + 200 MeV", labelfmt="L")

    hist_delta1.Add(hist_Nomin)
    hist_delta2.Add(hist_Nomin)

    ax1.set_xlim(xmin, xmax)

    # Draw line at y=1 in ratio panel
    line = root.TLine(ax1.get_xlim()[0], 1, ax1.get_xlim()[1], 1)
    ax2.plot(line)

    # Calculate and draw the ratio
    ratio_hist = hist_delta1.Clone("ratio_hist")
    ratio_hist_delta1 = hist_delta2.Clone("ratio_hist")

    ratio_hist.Divide(hist_Nomin)
    ratio_hist_delta1.Divide(hist_Nomin)

    ax2.plot(ratio_hist, linecolor=2)
    ax2.plot(ratio_hist_delta1, linecolor=4)

    ax2.set_xlim(xmin, xmax)
    ax2.set_ylim(0.94, 1.06)

    # Add extra space at top of plot to make room for labels
    ax1.add_margins(top=0.16)

    # Set axis titles
    ax1.set_ylabel("Events ", titleoffset=2.2)
    ax2.set_ylabel("dM/"+var+"", loc="centre")
    ax2.set_xlabel(labele, titleoffset=1)

    # Add extra space at top and bottom of ratio panel
    ax2.add_margins(top=0.1, bottom=0.1)

    # Go back to top axes to add labels
    ax1.cd()

    # Add the ATLAS Label
    aplt.atlas_label(text="Internal", loc="upper left")
    ax1.text(0.2, 0.84, "#sqrt{s} = 5 TeV, 256.827 pb^{-1}", size=22, align=13)

    # Add legend
    ax1.legend(loc=(0.78, 0.74, 0.6, 0.90))

    # Save the plot as a Png 
    fig.savefig("Outputs/hist-nom-200-MC-"+var+".png")







   