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
    X = []  # Xpos
    Y = []  # Ypos
    Sxx = []  # sig_xx
    Syy = []  # sig_yy
    Sxy = []  # sig_xy
    P = []  # hydrostatic pressure
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
    plt.plot(N, E2, 'bx', linewidth=0.3,
             markersize=1, label=r'$\lambda = 1.5$')
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
    plt.plot(N, SYY, 'rx', linewidth=0.3,
             markersize=1, label=r'$\lambda = 0.5$')
    plt.plot(N, SYY1, 'gx', linewidth=0.3,
             markersize=1, label=r'$\lambda = 1$')
    plt.plot(N, SYY2, 'bx', linewidth=0.3,
             markersize=1, label=r'$\lambda = 1.5$')
    plt.plot(N, SYY3, 'cx', linewidth=0.3,
             markersize=1, label=r'$\lambda = 5$')

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
    plt.plot(N, SXY, 'rx', linewidth=0.3,
             markersize=1, label=r'$\lambda = 0.5$')
    plt.plot(N, SXY1, 'gx', linewidth=0.3,
             markersize=1, label=r'$\lambda = 1$')
    plt.plot(N, SXY2, 'bx', linewidth=0.3,
             markersize=1, label=r'$\lambda = 1.5$')
    plt.plot(N, SXY3, 'cx', linewidth=0.3,
             markersize=1, label=r'$\lambda = 5$')

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
    plt.plot(N, P2, 'bx', linewidth=0.3,
             markersize=1, label=r'$\lambda = 1.5$')
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
    plt.plot(N, F2, 'b-', linewidth=0.6,
             markersize=1, label=r'$\lambda = 1.5$')
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


def totE_VS_T():
    #data collect from each Monte Carlo simulations when reached equlibirum states
    KE_n = [0.17234, 0.34468, 0.51702, 0.68936, 0.8617, 1.03404, 1.20638,1.37872, 1.55106, 1.7234, 1.89574, 2.06808]
    #this is obtained by a MD simulator:http://physics.weber.edu/schroeder/md/ since the potential energy generated appaers way too big
    PE_n = [-13.746, -12.631, -11.654, -11.832, -11.246, -11.134, -10.824, -10.231, -9.604, -8.876, -7.341, -6.743]

    E = []  # a list of E
    for i in range(len(KE_n)):
        e = KE_n[i] + PE_n[i]
        E.append(e)

    T = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60]

    fig = plt.figure()
    plt.plot(T, PE_n, 'bo-', linewidth=1, markersize=5, label=r'PE')
    leg = plt.legend(prop={'size': 10})
    leg.get_frame().set_alpha(0.5)
    xlabel = r'Temperature (K)'
    ylabel = r'Potential Energy (eV)'
    plt.ylim([-13, 0])
    plt.xlim([0, 90])
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(r'Potential Energy vs T')
    fig.savefig("PEvsT", dpi=300, bbox_inches='tight')

    fig = plt.figure()
    plt.plot(T, E, 'go-', linewidth=1, markersize=5, label=r'PE')
    leg = plt.legend(prop={'size': 10})
    leg.get_frame().set_alpha(0.5)
    xlabel = r'Temperature (K)'
    ylabel = r'Total Energy (eV)'
    #  plt.ylim ([-13,0])
    plt.xlim([0, 90])
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(r'Total Energy vs T')
    fig.savefig("EvsT", dpi=300, bbox_inches='tight')



def main():
    #-------------------------------------------------------------
    #plot the initial configuration and rdf
    #-------------------------------------------------------------
    drawconfig("config.txt", 'initial_config.png')
    graphrdf("rdf.txt", "initial_rdf.png")

    #iteration 5K
    drawconfig("5K_config.txt", "5K_config.png")
    graphrdf("5K_rdf.txt", "5K_rdf.png")

    #iteration 10K
    drawconfig("10K_config.txt", "10K_config.png")
    graphrdf("10K_rdf.txt", "10K_rdf.png")

    #iteration 15K
    drawconfig("15K_config.txt", "15K_config.png")
    graphrdf("15K_rdf.txt", "15K_rdf.png")

    #iteration 20K
    drawconfig("20K_config.txt", "20K_config.png")
    graphrdf("20K_rdf.txt", "20K_rdf.png")

    #iteration 25K
    drawconfig("25K_config.txt", "25K_config.png")
    graphrdf("25K_rdf.txt", "25K_rdf.png")

    #iteration 30K
    drawconfig("30K_config.txt", "30K_config.png")
    graphrdf("30K_rdf.txt", "30K_rdf.png")

    #iteration 35K
    drawconfig("35K_config.txt", "35K_config.png")
    graphrdf("35K_rdf.txt", "35K_rdf.png")

    #iteration 40K
    drawconfig("40K_config.txt", "40K_config.png")
    graphrdf("40K_rdf.txt", "40K_rdf.png")

    #iteration 45K
    drawconfig("45K_config.txt", "45K_config.png")
    graphrdf("45K_rdf.txt", "45K_rdf.png")

    #iteration 50k
    drawconfig("50K_config.txt", "50K_config.png")
    graphrdf("50K_rdf.txt", "50K_rdf.png")

    #iteration 55K
    drawconfig("55K_config.txt", "55K_config.png")
    graphrdf("55K_rdf.txt", "55K_rdf.png")

    #iteration 60K
    drawconfig("60K_config.txt", "60K_config.png")
    graphrdf("60K_rdf.txt", "60K_rdf.png")

    totE_VS_T()



main()
