import ROOT as rt

def process_pave(pave, font=42, size=0.06):

    pave.SetFillStyle(0)
    pave.SetBorderSize(0)
    #pave.SetTextColor(kBlack)
    pave.SetTextSize(size)
    pave.SetTextFont(font)
    
    return pave


def process_tf1(tf1, col=1, wid=2, sty=1):

    tf1.SetLineColor(col)
    tf1.SetLineStyle(sty)
    tf1.SetLineWidth(wid)


def process_hist(hist, size=1.4, col=1, style=20):

    rt.gPad.SetTickx()
    rt.gPad.SetTicky()
    hist.SetMarkerSize(size)
    hist.SetMarkerColor(col)
    hist.SetLineColor(col)
    hist.SetMarkerStyle(style) 
    hist.GetYaxis().SetTitleOffset(1.1) 
    hist.GetYaxis().SetTitleSize(0.055) 
    hist.GetYaxis().SetLabelSize(0.055) 
    hist.GetYaxis().SetLabelFont(42) 
    hist.GetXaxis().SetLabelFont(42) 
    hist.GetYaxis().SetTitleFont(42) 
    hist.GetXaxis().SetTitleFont(42) 
    hist.GetXaxis().SetTitleOffset(1.0) 
    hist.GetXaxis().SetTitleSize(0.055) 
    hist.GetXaxis().SetLabelSize(0.055) 

def process_hist2D(hist):

    #  hist.SetLogz() 
    hist.GetYaxis().SetTitleOffset(1.1) 
    hist.GetYaxis().SetTitleSize(0.055) 
    hist.GetYaxis().SetLabelSize(0.05) 
    hist.GetYaxis().SetLabelFont(42) 
    hist.GetXaxis().SetLabelFont(42) 
    hist.GetYaxis().SetTitleFont(42) 
    hist.GetXaxis().SetTitleFont(42) 
    hist.GetZaxis().SetLabelFont(42) 
    hist.GetZaxis().SetLabelSize(0.04) 
    hist.GetXaxis().SetTitleOffset(1.1) 
    hist.GetXaxis().SetTitleSize(0.055) 
    hist.GetXaxis().SetLabelSize(0.05) 
 
def process_tgraph(tgraph, font_size=2):

    tgraph.GetXaxis().CenterTitle(True)
    tgraph.GetYaxis().CenterTitle(True)
    tgraph.GetXaxis().SetTitleOffset(1.4)
    tgraph.GetYaxis().SetTitleOffset(1.6)
    #tgraph.GetXaxis().SetLabelSize(font_size)
    #tgraph.GetYaxis().SetLabelSize(font_size)

    tgraph.GetYaxis().SetTitleOffset(1.1) 
    tgraph.GetYaxis().SetTitleSize(0.055) 
    tgraph.GetYaxis().SetLabelSize(0.055) 
    tgraph.GetYaxis().SetLabelFont(42) 
    tgraph.GetXaxis().SetLabelFont(42) 
    tgraph.GetYaxis().SetTitleFont(42) 
    tgraph.GetXaxis().SetTitleFont(42) 
    tgraph.GetXaxis().SetTitleOffset(1.0) 
    tgraph.GetXaxis().SetTitleSize(0.055) 
    tgraph.GetXaxis().SetLabelSize(0.055) 