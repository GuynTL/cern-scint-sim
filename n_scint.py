#import necessary libraries
import numpy as np
import matplotlib.pyplot as plt
import ROOT
import math as mt
import time
##create constants dimesions for the scintillators etc !global variables!##
start = time.time()
# give scintillator dimensions(width, length, area)
s_w = 18.5
s_l = 340
s_a = s_l * s_w
sc_sep = 25

def main():
    ##get the relevant values from the user##
    
    # get the number of muons
    response = input("How many cosmic muons are falling per square centimeter per minute?\n")
    mins = input("How many minutes are we running the simulation for?\n")
    # mpscm is the muons per square centimeter per minute
    mpscm = float(response) * float(mins)
    # create a rectangle one meter above the center of the scintilator for point generation
    # get desired elevation for point generation (z)
    n_r = input("At what elevation are we generating muons? (centimeters)\n")
    # z coordinate for top and bottom planes
    z_t = sc_sep/2
    z_b = -z_t
    z_el = float(n_r)
    
    ##calculate the rectangle of generation edges based on the height above the top scintillator##
    
    # calculate slopes of lines allowing lowest angle of entry
    xslope = sc_sep/s_w
    yslope = sc_sep/s_l
    # calculate furthest y
    dy = (z_el/yslope)
    # multiply by 2 to get rectangle length
    r_l = dy * 2
    print("--- %s centimeters length ---" % (r_l))
    # calculate furthest x
    dx = (z_el/xslope)
    # multiply by two to get the rectangle length
    r_w = dx * 2
    print("--- %s centimeters width ---" % (r_w))
    #calculate area
    r_a = r_w * r_l
    print("--- %s centimeters squared ---" % (r_a))
    # create a raw number of descending muons for the rectangular area we are looking at
    muon_num = r_a * mpscm
    print("--- We are generating %s muons ---" % (muon_num))
    # create cosine distribution
    cos_d = ROOT.TF1( 'cos_d', '(cos(x))**2', -np.pi/2, np.pi/2 )
    
    ##set up counters##
    
    tot = 0
    triggers = 0
    out_s = 0
    
    ##create histograms to be filled.##
    
    h1 = ROOT.TH1D("h1", "Theta Distribution", 100, -np.pi, np.pi)
    h2 = ROOT.TH1D("h2", "Phi Distribution", 100, -np.pi, np.pi)
    denominator = ROOT.TH2D("h3", "X and Y Distribution", 250, -r_l/4, r_l/4, 250, -r_l/4, r_l/4)
    numerator   = ROOT.TH2D("h4", "Xin Yin Distribution", 250, -r_l/4, r_l/4, 250, -r_l/4, r_l/4)
    h5 = ROOT.TH2D("h5", "XY in/out Distribution", 100, -s_w, s_w, 100, -r_l, r_l)
    
    ##iterate to make relevant values for each muon##
    
    for i in range (int(round_up(muon_num))):
        # create random x and y coordinates
        (x_rand, y_rand) = np.random.uniform(-r_w/2,r_w/2), np.random.uniform(-r_l/2,r_l/2)
        # fill the histogram with values
        denominator.Fill(x_rand,y_rand)
        # tally that a muon has been generated
        tot += 1
        # create angles for each point
        theta = cos_d.GetRandom()
        phi = cos_d.GetRandom()
        # fill the histograms with values
        h1.Fill(theta)
        h2.Fill(phi)
        # get the slope from the angle
        mx = mt.tan(theta)
        my = mt.tan(phi)
        # get the x and y coordinates on the top and bottom layers of the scintillators for each trigger
        top_x = get_val(x_rand, mx, z_t, z_el)
        bot_x = get_val(x_rand, mx, z_b, z_el)
        top_y = get_val(y_rand, my, z_t, z_el)
        bot_y = get_val(y_rand, my, z_b, z_el)
        # see if the points pass through both
        t_pass = intersection(top_x,top_y)
        b_pass = intersection(bot_x,bot_y)
        # count the number of triggers
        if t_pass == True and b_pass == True:
            triggers += 1
            # fill the histogram with values
            numerator.Fill(x_rand,y_rand)
        else:
            out_s += 1
    percentage = (triggers/tot) * 100
    # convert the raw muon pass through to hertz
    seconds = float(mins) * 60
    Hertz = triggers/seconds
    
    ##Draw the histograms##
    
    # theta distribution
    thcan = ROOT.TCanvas("c1","c1",1000, 1000)
    h1.GetXaxis().SetTitle("Theta")
    h1.Draw()
    thcan.Print("theta_dist.pdf")
    # phi distribution
    phcan = ROOT.TCanvas("c2","c2",1000, 1000)
    h2.GetXaxis().SetTitle("Phi")
    h2.Draw()
    phcan.Print("phi_dist.pdf")
    # overall points distribution
    xycan = ROOT.TCanvas("c3","c3",1000, 1000)
    denominator.GetXaxis().SetTitle("X-coordinate")
    denominator.GetYaxis().SetTitle("Y-coorcinate")
    denominator.Draw("colz")
    xycan.Print("XY.pdf")
    # points in distribution
    incan = ROOT.TCanvas("c4","c4",1000, 1000)
    numerator.GetXaxis().SetTitle("X-coordinate")
    numerator.GetYaxis().SetTitle("Y-coorcinate")
    numerator.Draw("colz")
    incan.Print("XYin.pdf")
    # divide the xy histgram by the xy in histogram, make efficiency plot
    # clone the histogram
    h6 = numerator.Clone("h6")
    h6.GetXaxis().SetTitle("X-coordinate")
    h6.GetYaxis().SetTitle("Y-coorcinate")
    h5 = h6.Divide(denominator)
    inoutcan = ROOT.TCanvas("c5","c5",1000, 1000)
    h6.Draw("colz")
    
    ##put a box on the graphs for the scintillators##

    inoutcan.Update()
    # create lines to denote the scintillator

    # bottom line
    b_line = ROOT.TLine(-s_w/2,-s_l/2,s_w/2,-s_l/2)
    b_line.SetLineColor(ROOT.kBlack)
    b_line.SetLineWidth(2)
    b_line.Draw()
    # top line
    t_line = ROOT.TLine(-s_w/2,s_l/2,s_w/2,s_l/2)
    t_line.SetLineColor(ROOT.kBlack)
    t_line.SetLineWidth(2)
    t_line.Draw()
    # left line
    l_line = ROOT.TLine(-s_w/2,-s_l/2,-s_w/2,s_l/2)
    l_line.SetLineColor(ROOT.kBlack)
    l_line.SetLineWidth(2)
    l_line.Draw()
    # right line
    r_line = ROOT.TLine(s_w/2,-s_l/2,s_w/2,s_l/2)
    r_line.SetLineColor(ROOT.kBlack)
    r_line.SetLineWidth(2)
    r_line.Draw()
    # print the histogram
    inoutcan.Print("Efficiency.pdf")

    # display relevant values
    # display relevant values
    print("--- %s muons total ---" % (tot))
    print("--- %s triggers ---" % (triggers))
    print("--- %s misses ---" % (out_s))
    print("--- %s percent of muons passed through ---" % (percentage))
    print("--- %s hertz ---" % (Hertz))
    print("--- The runtime was %s seconds ---" % (time.time() - start))

# function to create the values on the line
def get_val(x, m, z, z_i):
    new_x = x + m * (z - z_i)
    return new_x

# function that sees if the point intersected
def intersection(x,y):
    #set up maximum and minimum xs for the representative plane of the scintillator
    global s_w
    global s_l
    xmax = s_w/2
    xmin = -xmax
    ymax = s_l/2
    ymin = -ymax
    # compare the values to see if it went through the plane
    if xmin <= x <= xmax and ymin <= y <= ymax:
        return True
    else:
        return False
    
# function that rounds up
def round_up(n, decimals=0):
    multiplier = 10 ** decimals
    return mt.ceil(n * multiplier) / multiplier

if __name__ == "__main__":
    main()
        
