import math
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import scipy.stats as stats
import collections
import matplotlib as mpl
from matplotlib import rc
from itertools import groupby
rc('font', **{'family': 'serif', 'serif': ['Palatino']})
rc('text', usetex=True)


def drawconfig(txtname, filename):
    '''use data from the textfile to draw the configuration'''
    X = []  #Xpos
    Y = []  #Ypos
    Sxx = []  #sig_xx
    Syy = []  #sig_yy
    Sxy = []  #sig_xy
    P = [] #hydrostatic pressure
    with open(txtname, "r") as file:
        for line in file:
            nums = line.split()
            x = float(nums[0])  
            y = float(nums[1])  
            sxx = float(nums[2])  
            syy = float(nums[3]) 
            sxy = float(nums[4]) 
            p = -0.5*(syy+sxx)

            X.append(x)  
            Y.append(y)  
            Sxx.append(sxx)
            Syy.append(syy)
            Sxy.append(sxy)
            P.append(p)

    fig = plt.figure()
    plt.plot(X, Y, 'ro', linewidth=0.8)
    plt.ylim([-4, 80])
    plt.xlim([-4, 90])
    xlabel = r'$x_{coordinate}\ unit =\AA $'
    ylabel = r'$y_{coordinate}\ unit =\AA $'
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title('Configuration of Atoms')
    fig.savefig(filename, dpi=180, bbox_inches='tight')

    #Plot the stresses as heatmap
    #---------------simgaxx--------------------------------
    fig, ax = plt.subplots(1)
    plt.scatter(X, Y, s=150, marker='8', c=Sxx, alpha=0.5)
    cbar = plt.colorbar()
    cbar.set_label(r"$eV/{\AA}^2$")
    xlabel = r'X_{pos}'
    ylabel = r'Y_{pos}'
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.locator_params(axis='y', nticks=3)
    plt.locator_params(axis='x', nticks=8)
    fname = "sig_xx" + filename
    titlename = r"$\sig_{xx}$"
    plt.title(titlename)
    fig.savefig(fname, dpi=180)

    #---------------simgayy--------------------------------
    fig, ax = plt.subplots(1)
    plt.scatter(X, Y, s=150, marker='8', c=Syy, alpha=0.5)
    cbar = plt.colorbar()
    cbar.set_label(r"$eV/{\AA}^2$")
    xlabel = r'X_{pos}'
    ylabel = r'Y_{pos}'
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.locator_params(axis='y', nticks=3)
    plt.locator_params(axis='x', nticks=8)
    fname = "sig_yy" + filename
    titlename = r"$\sig_{yy}$"
    plt.title(titlename)
    fig.savefig(fname, dpi=180)

    #---------------simgaxy--------------------------------
    fig, ax = plt.subplots(1)
    plt.scatter(X, Y, s=150, marker='8', c=Sxy, alpha=0.5)
    cbar = plt.colorbar()
    cbar.set_label(r"$eV/{\AA}^2$")
    xlabel = r'X_{pos}'
    ylabel = r'Y_{pos}'
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.locator_params(axis='y', nticks=3)
    plt.locator_params(axis='x', nticks=8)
    fname = "sig_xy" + filename
    titlename = r"$\sig_{xy}$"
    plt.title(titlename)
    fig.savefig(fname, dpi=180)

    #---------------  P -----------------------------------
    fig, ax = plt.subplots(1)
    plt.scatter(X, Y, s=150, marker='8', c=P, alpha=0.5)
    cbar = plt.colorbar()
    cbar.set_label(r"$eV/{\AA}^2$")
    xlabel = r'X_{pos}'
    ylabel = r'Y_{pos}'
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.locator_params(axis='y', nticks=3)
    plt.locator_params(axis='x', nticks=8)
    fname = "hydrop" + filename
    titlename = r"$\sig_{xy}$"
    plt.title(titlename)
    fig.savefig(fname, dpi=180)


def graphrdf(txtname, filename):
    '''take the rdf.txt to plot the rdf function'''

    R = []  # a list of distance distribution
    G = []  # a list of distance distribution
    with open(txtname, "r") as file:
        for line in file:
            nums = line.split()
            r = float(nums[0])  # take the value
            g = float(nums[1])  # take the value
            R.append(r)  # append r value into r
            G.append(g)  # append g(r) value into G
    fig = plt.figure()

    if txtname == "rdf.txt":
        plt.stem(R, G, linefmt='b-', linewidth=0.4,
                 markerfmt='rx', basefmt='r-', label=r'RDF')
    else:
        plt.plot(R, G, 'rx-', linewidth=0.4, markersize=1, label=r'RDF')
    xlabel = r'$r\ unit =\AA $'
    ylabel = r'$g(r)\ $'
    plt.title('Radial distribution function g(r)')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    fig.savefig(filename, dpi=180, bbox_inches='tight')


def plot_vsN(tx1, tx2, tx3, tx4, tx5, EN, SxxN, SyyN, SxyN, PN, FN, E10, F10):
    '''make plots for stress vs N; potential energy vs N
    maximum Force vs N and hydrostatic pressure vs N'''

    #-------------------------------------------------------------
    # step 1 load files
    #-------------------------------------------------------------
    #--------------------lambda = 0.5 ----------------------------
    N = []  
    E = []  
    SXX = [] 
    SYY = [] 
    SXY = []  
    P = []  
    F = [] 
    with open(tx1, "r") as file:
        for line in file:
            nums = line.split()
            n = float(nums[0])  
            e = float(nums[1])  
            sxx = float(nums[2])  
            syy = float(nums[3])  
            sxy = float(nums[4])  
            p = float(nums[5])  
            f = float(nums[6])  

            N.append(n)  
            E.append(e)  
            SXX.append(sxx)  
            SYY.append(syy)  
            SXY.append(sxy) 
            P.append(p)
            F.append(f)

    #--------------------lambda = 1 ------------------------------
    N1 = []  # a list of interation
    E1 = []  # a list of PE
    SXX1 = []  # sig_xx
    SYY1 = []  # sig_yy
    SXY1 = []  # sig_xy
    P1 = []  # hydrostatic pressure
    F1 = []  # force
    with open(tx2, "r") as file:
        for line in file:
            nums = line.split()
            n = float(nums[0])  
            e = float(nums[1])  
            sxx = float(nums[2])  
            syy = float(nums[3])  
            sxy = float(nums[4]) 
            p = float(nums[5])  
            f = float(nums[6]) 

            N1.append(n)  
            E1.append(e)  
            SXX1.append(sxx)  
            SYY1.append(syy)  
            SXY1.append(sxy)  
            P1.append(p)
            F1.append(f)

    #--------------------lambda = 1.5 ----------------------------
    N2 = []  
    E2 = []  
    SXX2 = []  
    SYY2 = [] 
    SXY2 = []  
    P2 = []  
    F2 = []  
    with open(tx3, "r") as file:
        for line in file:
            nums = line.split()
            n = float(nums[0])  
            e = float(nums[1])  
            sxx = float(nums[2]) 
            syy = float(nums[3])  
            sxy = float(nums[4])  
            p = float(nums[5])  
            f = float(nums[6])  

            N2.append(n)  
            E2.append(e)  
            SXX2.append(sxx)  
            SYY2.append(syy)  
            SXY2.append(sxy)  
            P2.append(p)
            F2.append(f)

    #--------------------lambda = 5 ------------------------------
    N3 = []  
    E3 = []  
    SXX3 = [] 
    SYY3 = []  
    SXY3 = []  
    P3 = []  
    F3 = []  
    with open(tx4, "r") as file:
        for line in file:
            nums = line.split()
            n = float(nums[0])  
            e = float(nums[1])  
            sxx = float(nums[2]) 
            syy = float(nums[3])  
            sxy = float(nums[4])  
            p = float(nums[5])  
            f = float(nums[6])  

            N3.append(n)  
            E3.append(e)  
            SXX3.append(sxx)  
            SYY3.append(syy) 
            SXY3.append(sxy)
            P3.append(p)
            F3.append(f)

    #--------------------lambda = 10 ------------------------------
    N4 = []  # a list of interation
    E4 = []  # a list of PE
    F4 = []  # force
    with open(tx5, "r") as file:
        for line in file:
            nums = line.split()
            n = float(nums[0])  
            e = float(nums[1])  
            p = float(nums[5])  
            f = float(nums[6])  

            N4.append(n)  
            E4.append(e)  
            F4.append(f)

    #-------------------------------------------------------------
    #step 2
    #-------------------------------------------------------------

    #plot the plots

    #--------------------plot Potential energy vs N---------------

    fig = plt.figure()

    plt.plot(N, E, 'rx', linewidth=0.3, markersize=1, label=r'$\lambda = 0.5$')
    plt.plot(N, E1, 'gx', linewidth=0.3, markersize=1, label=r'$\lambda = 1$')
    plt.plot(N, E2, 'bx', linewidth=0.3,markersize=1, label=r'$\lambda = 1.5$')
    plt.plot(N, E3, 'cx', linewidth=0.3, markersize=1, label=r'$\lambda = 5$')

    leg = plt.legend(prop={'size': 10})
    leg.get_frame().set_alpha(0.5)
    xlabel = r'Number of Iterations'
    ylabel = r'Potential Energy (eV)'
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(r'Potential Energy vs N')
    fig.savefig(EN, dpi=180, bbox_inches='tight')

    #--------------------plot stresses vs N ----------------------
    # Sxx
    fig = plt.figure()
    plt.plot(N, SXX, 'rx', linewidth=0.3,
             markersize=1, label=r'$\lambda = 0.5$')
    plt.plot(N, SXX1, 'gx', linewidth=0.3,
             markersize=1, label=r'$\lambda = 1$')
    plt.plot(N, SXX2, 'bx', linewidth=0.3,
             markersize=1, label=r'$\lambda = 1.5$')
    plt.plot(N, SXX3, 'cx', linewidth=0.3,
             markersize=1, label=r'$\lambda = 5$')

    leg = plt.legend(prop={'size': 10})
    leg.get_frame().set_alpha(0.5)
    xlabel = r'Number of Iterations'
    ylabel = r'$\sig_{xx}$ ($eV/\AA^{2}$)'
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(r'$\sig_{xx}$ vs N')
    fig.savefig(SxxN, dpi=180, bbox_inches='tight')

    # Syy
    fig = plt.figure()
    plt.plot(N, SYY, 'rx', linewidth=0.3,markersize=1, label=r'$\lambda = 0.5$')
    plt.plot(N, SYY1, 'gx', linewidth=0.3,markersize=1, label=r'$\lambda = 1$')
    plt.plot(N, SYY2, 'bx', linewidth=0.3,markersize=1, label=r'$\lambda = 1.5$')
    plt.plot(N, SYY3, 'cx', linewidth=0.3,markersize=1, label=r'$\lambda = 5$')

    leg = plt.legend(prop={'size': 10})
    leg.get_frame().set_alpha(0.5)
    xlabel = r'Number of Iterations'
    ylabel = r'$\sig_{yy}$ ($eV/\AA^{2}$)'
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(r'$\sig_{yy}$ vs N')
    fig.savefig(SyyN, dpi=180, bbox_inches='tight')

    # Sxy
    fig = plt.figure()
    plt.plot(N, SXY, 'rx', linewidth=0.3, markersize=1, label=r'$\lambda = 0.5$')
    plt.plot(N, SXY1, 'gx', linewidth=0.3, markersize=1, label=r'$\lambda = 1$')
    plt.plot(N, SXY2, 'bx', linewidth=0.3, markersize=1, label=r'$\lambda = 1.5$')
    plt.plot(N, SXY3, 'cx', linewidth=0.3, markersize=1, label=r'$\lambda = 5$')

    leg = plt.legend(prop={'size': 10})
    leg.get_frame().set_alpha(0.5)
    xlabel = r'Number of Iterations'
    ylabel = r'$\sig_{xy}$ ($eV/\AA^{2}$)'
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(r'$\sig_{xy}$ vs N')
    fig.savefig(SxyN, dpi=180, bbox_inches='tight')

    # P
    fig = plt.figure()
    plt.plot(N, P, 'rx', linewidth=0.3, markersize=1, label=r'$\lambda = 0.5$')
    plt.plot(N, P1, 'gx', linewidth=0.3, markersize=1, label=r'$\lambda = 1$')
    plt.plot(N, P2, 'bx', linewidth=0.3,markersize=1, label=r'$\lambda = 1.5$')
    plt.plot(N, P3, 'cx', linewidth=0.3, markersize=1, label=r'$\lambda = 5$')

    leg = plt.legend(prop={'size': 10})
    leg.get_frame().set_alpha(0.5)
    xlabel = r'Number of Iterations'
    ylabel = r'Hydrostatic pressure ($eV/\AA^{2}$)'
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(r'Hydrostatic pressure vs N')
    fig.savefig(PN, dpi=180, bbox_inches='tight')

    #--------------------plot F_max vs N -------------------------
    fig = plt.figure()
    plt.plot(N, F, 'r-', linewidth=0.6, markersize=1, label=r'$\lambda = 0.5$')
    plt.plot(N, F1, 'g-', linewidth=0.6, markersize=1, label=r'$\lambda = 1$')
    plt.plot(N, F2, 'b-', linewidth=0.6, markersize=1, label=r'$\lambda = 1.5$')
    plt.plot(N, F3, 'c-', linewidth=0.6, markersize=1, label=r'$\lambda = 5$')

    leg = plt.legend(prop={'size': 10})
    leg.get_frame().set_alpha(0.5)
    xlabel = r'Number of Iterations'
    ylabel = r'Maximum Force ($eV/\AA$)'
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(r'Maximum Force vs N')
    fig.savefig(FN, dpi=180, bbox_inches='tight')

    #--------------------plot F_10 vs N ---------------------------
    fig = plt.figure()
    plt.plot(N, F4, 'r-', linewidth=0.6, markersize=1, label=r'$\lambda = 10$')
    leg = plt.legend(prop={'size': 10})
    leg.get_frame().set_alpha(0.5)
    xlabel = r'Number of Iterations'
    ylabel = r'Maximum Force ($eV/\AA$)'
    plt.ylim([0, 0.14])
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(r'Maximum Force vs N')
    fig.savefig(F10, dpi=180, bbox_inches='tight')

    #--------------------plot E_10 vs N ---------------------------
    fig = plt.figure()
    plt.plot(N, E4, 'r-', linewidth=0.6, markersize=1, label=r'$\lambda = 10$')
    leg = plt.legend(prop={'size': 10})
    leg.get_frame().set_alpha(0.5)
    xlabel = r'Number of Iterations'
    ylabel = r'Potential Energy (eV)'
    plt.ylim([-12, 12])
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(r'Potential Energy vs N')
    fig.savefig(E10, dpi=180, bbox_inches='tight')


def main():
    #-------------------------------------------------------------
    #plot the initial configuration and rdf
    #-------------------------------------------------------------
    drawconfig("config.txt", 'initial_config.png')
    graphrdf("rdf.txt", "initial_rdf.png")

    #-------------------------------------------------------------
    #plot the final configuration and rdf of the 4 experiments
    #with different lambdas
    #-------------------------------------------------------------

    #----------------lambda = 0.5--------------------
    drawconfig("l=0p5_config.txt", "l=0p5_config.png")
    graphrdf("l=0p5_rdf.txt", "l=0p5_rdf.png")

    #----------------lambda = 1 ---------------------
    drawconfig("l=1_config.txt", "l=1_config.png")
    graphrdf("l=1_rdf.txt", "l=1_rdf.png")

    #----------------lambda = 1.5--------------------
    drawconfig("l=1p5_config.txt", "l=1p5_config.png")
    graphrdf("l=1p5_rdf.txt", "l=1p5_rdf.png")

    #----------------lambda = 5 ---------------------
    drawconfig("l=5_config.txt", "l=5_config.png")
    graphrdf("l=5_rdf.txt", "l=5_rdf.png")

    #----------------lambda = 10 ---------------------
    drawconfig("l=10_config.txt", "l=10_config.png")
    graphrdf("l=10_rdf.txt", "l=10_rdf.png")

    #------------------------------------------------
    #study the evolution lambda = 1.5
    #------------------------------------------------

    #iteration 100
    drawconfig("100config.txt", "100config.png")
    graphrdf("100rdf.txt", "100rdf.png")

    #iteration 1000
    drawconfig("1000config.txt", "1000config.png")
    graphrdf("1000rdf.txt", "1000rdf.png")

    #iteration 5000
    drawconfig("5000config.txt", "5000config.png")
    graphrdf("5000rdf.txt", "5000rdf.png")

    #iteration 10000
    drawconfig("10000config.txt", "10000config.png")
    graphrdf("10000rdf.txt", "10000rdf.png")

    #iteration 20000
    drawconfig("20000config.txt", "20000config.png")
    graphrdf("20000rdf.txt", "20000rdf.png")

    #iteration 50000
    drawconfig("50000config.txt", "50000config.png")
    graphrdf("50000rdf.txt", "50000rdf.png")

    #------------------------------------------------
    #plot variables vs N
    #------------------------------------------------
    plot_vsN("l=0p5_vsN.txt", "l=1_vsN.txt", "l=1p5_vsN.txt",
             "l=5_vsN.txt", "l=10_vsN.txt", "E_vs_N.png",
             "Sxx_vs_N.png", "Syy_vs_N.png", "Sxy_vs_N.png",
             "P_vs_N.png", "F_vs_N.png", "P10_vs_N.png", "F10_vs_N.png")


main()
