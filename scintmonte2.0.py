#import necessary libraries
import numpy as np
import matplotlib.pyplot as plt
import ROOT
ROOT.gROOT.SetBatch(True)
import math as mt
import time
##create constants dimesions for the scintillators !global variables!##
start = time.time()

def main():
    
## get the relevant values from the user ##
    # give scintillator dimensions(width, length, area)
    s_w = 15
    s_l = 15
    s_a = s_l * s_w
    # get desired elevation for point generation (z)
    # z_el = float(input("At what elevation are we generating muons? (centimeters)\n"))
    z_el = 500
    # get the center points from the user
    # y_center = float(input("What is the center y coordinate of our scintillator setup? (centimeters)\n"))
    y_center = 0
    # x_c_1 = float(input("What is the center x coordinate of our first scintillator? (centimeters)\n"))
    x_c_1 = 0
    # z_c_1 = float(input("What is the center z coordinate of our first scintillator? (centimeters)\n"))
    z_c_1 = 300
    # x_c_2 = float(input("What is the center x coordinate of our second scintillator? (centimeters)\n"))
    x_c_2 = 20
    # z_c_2 = float(input("What is the center z coordinate of our second scintillator? (centimeters)\n"))
    z_c_2 = 320
    # mpscm is the muons per square centimeter
    mins = float(input("How many minutes are we running the simulation for?\n"))
    mpscm = float(input("How many cosmic muons are falling per square centimeter per minute?\n")) * mins
    # calculate scintillator seperation
    sc_sep = x_c_2 - x_c_1
    print("--- %s seperation ---" % (sc_sep))
    print("--- %s area ---" % (s_a))
    # figure out how the scintillators are positioned.
    x_or = orientation_check(x_c_2, x_c_1)
    if x_or == False:
        print("Please place the second scintillator to the right of the first!")
        return 1
    # set up the upper and lower bounds for our scintillators
    z_edgetop1 = z_c_1 + s_w/2
    print("--- %s top edge 1 ---" % (z_edgetop1))
    z_edgebot1 = z_c_1 - s_w/2
    print("--- %s bot edge 1  ---" % (z_edgebot1))
    z_edgetop2 = z_c_2 + s_w/2
    print("--- %s top edge 2 ---" % (z_edgetop2))
    z_edgebot2 = z_c_2 - s_w/2
    print("--- %s bot edge 2  ---" % (z_edgebot2))
    z_or = orientation_check(z_edgebot2, z_edgetop1)
    if z_or == False:
        print("Please place the second scintillator so that the edges of the two dont overlap in the z plane!")
        return 1
    y_edgel = y_center - s_l/2
    print("--- %s left edge ---" % (y_edgel))
    y_edger = y_center + s_l/2
    print("--- %s right edge ---" % (y_edger))
    #find centerpoint of the scintillators
    x_center = (x_c_1 + (sc_sep/2))
    print("--- %s x coordinate center ---" % (x_center))
    # z center point
    z_center = z_edgetop1 + ((z_edgebot2-z_edgetop1)/2)
    print("--- %s z coordinate center ---" % (z_center))
    
## calculate the rectangle of generation edges based on the height above the top scintillator ##
    
    # calculate slopes of lines allowing lowest angle of entry according to the orientation
    xzslope = (z_edgebot2 - z_edgetop1)/sc_sep
    print("--- %s xz slope ---" % (xzslope))
    yzslope = (z_edgebot2 - z_edgetop1)/s_l
    print("--- %s yz slope ---" % (yzslope))
    #find the intercept of the lowest slope line
    yz_int = find_int(yzslope, y_center, z_center)
    print("--- %s yz intercept ---" % (yz_int))
    # calculate furthest y
    dy = ((z_el - yz_int)/yzslope)
    print("--- %s furthest y ---" % (dy))
    # multiply by two to get rectangle length
    r_l = dy * 2
    print("--- %s centimeters length ---" % (r_l))
    #find the intercept of the lowest slope line
    xz_int = find_int(xzslope, x_c_1, z_edgetop1)
    print("--- %s xz intercept ---" % (xz_int))
    # calculate furthest x
    dx = ((z_el - xz_int)/xzslope)
    print("--- %s furthest x ---" % (dx))
    # multiply by two to get the rectangle width
    r_w = dx * 2
    print("--- %s centimeters width ---" % (r_w))
    # calculate area
    r_a = r_w * r_l
    print("--- %s centimeters squared ---" % (r_a))
    # create a raw number of descending muons for the rectangular area we are looking at
    muon_num = r_a * mpscm
    print("--- We are generating %s muons ---" % (muon_num))
    # create cosine distribution
    cos_d = ROOT.TF1( 'cos_d', '(cos(x))**2', -np.pi/2, np.pi/2 )
    
## set up counters ##
    
    tot = 0
    triggers = 0
    out_s = 0
    
## create histograms to be filled ##
    
    h1 = ROOT.TH1D("h1", "Theta Distribution", 100, -np.pi, np.pi)
    h2 = ROOT.TH1D("h2", "Phi Distribution", 100, -np.pi, np.pi)
    denominator = ROOT.TH2D("h3", "X and Y Distribution", 250, -r_l/4, r_l/4, 250, -r_l/4, r_l/4)
    numerator   = ROOT.TH2D("h4", "Xin Yin Distribution", 250, -r_l/4, r_l/4, 250, -r_l/4, r_l/4)
    h5 = ROOT.TH2D("h5", "XY in/out Distribution", 100, -s_w, s_w, 100, -r_l, r_l)
    
## iterate to make relevant values for each muon ##
    
    for i in range (int(round_up(muon_num))):
        # create random x and y coordinates
        (x_rand, y_rand) = np.random.uniform( -r_w + x_center, r_w + x_center), np.random.uniform(-r_l/2 + y_center,r_l/2 + y_center)
        # fill the histogram with values
        denominator.Fill(x_rand,y_rand)
        # create angles for each point
        theta = cos_d.GetRandom()
        phi = cos_d.GetRandom()
        # generate a slope for the muon
        # get the slope from the angle
        mx = mt.tan(theta)
        my = mt.tan(phi)
        
    ## figure out if the line that is the scintillator and the line that is the muon path intersect ##
        # tally that a muon has been generated
        tot += 1
        # fill the histograms with values
        h1.Fill(theta)
        h2.Fill(phi)
        # get the x and y and z coordinates on the top and bottom layers of the scintillators for each trigger
        z_1 = get_newcoord(z_el, mx, x_c_1, x_rand)
        z_2 = get_newcoord(z_el, mx, x_c_2, x_rand)
        y_t2 = get_newcoord(y_rand, my, z_edgetop2, z_el)
        y_b2 = get_newcoord(y_rand, my, z_edgebot2, z_el)
        y_t1 = get_newcoord(y_rand, my, z_edgetop1, z_el)
        y_b1 = get_newcoord(y_rand, my, z_edgebot1, z_el)
        # see if the points pass through both
        f_pass = intersection(z_1, y_t1, y_b1, z_edgetop1, z_edgebot1, y_center + s_l/2, y_center - s_l/2)
        s_pass = intersection(z_2, y_t2, y_b2, z_edgetop2, z_edgebot2, y_center + s_l/2, y_center - s_l/2)
        # count the number of triggers
        if f_pass == True and s_pass == True:
            triggers += 1
            # fill the histogram with values
            numerator.Fill(x_rand,y_rand)
        else:
            out_s += 1
    percentage = (triggers/tot) * 100
    # convert the raw muon pass through to hertz
    seconds = float(mins) * 60
    Hertz = triggers/seconds
    
## Draw the histograms ##
    
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
    denominator.GetYaxis().SetTitle("Y-coordinate")
    denominator.Draw("colz")
    xycan.Print("XY.pdf")
    # points in distribution
    incan = ROOT.TCanvas("c4","c4",1000, 1000)
    numerator.GetXaxis().SetTitle("X-coordinate")
    numerator.GetYaxis().SetTitle("Y-coordinate")
    numerator.Draw("colz")
    incan.Print("XYin.pdf")
    # divide the xy histgram by the xy in histogram, make efficiency plot
    # clone the histogram
    h6 = numerator.Clone("h6")
    h6.GetXaxis().SetTitle("X-coordinate")
    h6.GetYaxis().SetTitle("Y-coordinate")
    h5 = h6.Divide(denominator)
    inoutcan = ROOT.TCanvas("c5","c5",1000, 1000)
    inoutcan.Draw()
    #inoutcan.cd()
    h6.Draw("colz")
    inoutcan.Print("Efficiency.pdf")

    #display relevant values
    print("--- %s muons total ---" % (tot))
    print("--- %s triggers ---" % (triggers))
    print("--- %s misses ---" % (out_s))
    print("--- %s percent of muons passed through ---" % (percentage))
    print("--- %s hertz ---" % (Hertz))
    print("--- The runtime was %s seconds ---" % (time.time() - start))

# function to create the values on the line
def get_newcoord(x, m, z, z_i):
    new_x = x + m * (z - z_i)
    return new_x

#function that finds intercept
def find_int(m,x,z):
    b = z - (m*x)
    return b

# function that sees if the point intersected
def intersection(z, y, y2, z_max, z_min, y_max, y_min):
    # compare the values to see if it went through the plane
    if z_min <= z <= z_max and y_min <= y <= y_max and y_min <= y2 <= y_max:
        return True
    else:
        return False
    
# function that rounds up
def round_up(n, decimals=0):
    multiplier = 10 ** decimals
    return mt.ceil(n * multiplier) / multiplier
                             
# function that checks to see what orientation is of the scintillators is
def orientation_check(x,y):
    if x > y:
        return True
    else:
        return False
if __name__ == "__main__":
    main()
