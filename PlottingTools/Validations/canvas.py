from ROOT import TCanvas, SetOwnership

def newCanvas():
    H_ref = 600;
    W_ref = 800;
    H  = H_ref
    W = W_ref

    # references for T, B, L, R
    T = 0.08*H_ref
    B = 0.12*H_ref
    L = 0.12*W_ref
    R = 0.04*W_ref

    c = TCanvas("c","c",50,50,W,H)
    c.Clear()
    c.SetLeftMargin( L/W )
    c.SetRightMargin( R/W )
    c.SetTopMargin( T/H )
    c.SetBottomMargin( B/H )

    SetOwnership(c, False)
    return c

def newCanvas2D():
    H_ref = 800;
    W_ref = 1000;
    H  = H_ref
    W = W_ref

    # references for T, B, L, R
    T = 0.08*H_ref
    B = 0.12*H_ref
    L = 0.12*W_ref
    R = 0.12*W_ref

    c = TCanvas("c","c",50,50,W,H)
    c.Clear()
    c.SetLeftMargin( L/W )
    c.SetRightMargin( R/W )
    c.SetTopMargin( T/H )
    c.SetBottomMargin( B/H )

    SetOwnership(c, False)
    return c
