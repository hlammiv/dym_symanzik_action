import numpy as np
import os
import matplotlib.pyplot as plt
import scipy
from scipy import signal
# opens file (the path will need to be changed according to workstation)
# the path to out_b. . . files in the workstation: /home/guest/dym_par_adj/hotstartdata
# on my computer: C:\\Users\\annac\\OneDrive\\Desktop\\SIST Project\\Computer Code\\data\\data
themeans = []
thestandarddevs = []
beta1 = []
beta2 = []
data = []
dataformatrix = []
# files = os.listdir("/home/guest/dym_par_adj/hotstartdata")
files = os.listdir("/home/guest/dym_par_adj/datatofindcriticalval")
# print(files)
for file in files:
    # f = open(os.path.join("/home/guest/dym_par_adj/hotstartdata", file), 'r')
    f = open(os.path.join("/home/guest/dym_par_adj/datatofindcriticalval", file), 'r')
    file = file.split('b1')

    file = file[1].split('b2')

    file1 = file[0].split('_')
    file2 = file[1].split('_')

    betaone = 3 * float(file1[0])
    betatwo = 3 * float(file2[0])
    # print(betaone)
    # print(betatwo)
    beta1.append(3 * float(file1[0]))
    beta2.append(3 * float(file2[0]))
    measurementlines = []
    for line in f:
        if 'GMES: ' in line:
            measurementlines.append(line)
# splits each line at a space and creates an array of arrays of strings
    for i in range(len(measurementlines)):
        measurementlines[i] = measurementlines[i].split(' ')

        simpplaqs = []
    for i in range(len(measurementlines)):
        simpplaqs.append((float(measurementlines[i][4]))/3 + 1)

    mean = sum(simpplaqs)/len(simpplaqs)
    sd = np.std(simpplaqs)

    themeans.append(mean)
    thestandarddevs.append(sd)
    data.append([betaone, betatwo, mean, sd])
    dataformatrix.append([int(float(file1[0]) * 10), int(float(file2[0]) * 10) + 17, mean, sd])
        # print('yes')
# print(beta2)
# print(len(themeans))
# print(len(dataformatrix))

# files = os.listdir("/home/guest/dym_par_adj/forthefirstcontourplotnegative")
# print(files)
# for file in files:
#     # f = open(os.path.join("/home/guest/dym_par_adj/hotstartdata", file), 'r')
#     f = open(os.path.join("/home/guest/dym_par_adj/forthefirstcontourplotnegative", file), 'r')
#     file = file.split('b1')
#
#     file = file[1].split('b2')
#
#     file1 = file[0].split('_')
#     file2 = file[1].split('_')
#
#     betaone = 3 * float(file1[0])
#     betatwo = 3 * float(file2[0])
#     # print(betaone)
#     # print(betatwo)
#     beta1.append(3 * float(file1[0]))
#     beta2.append(3 * float(file2[0]))
#     measurementlines = []
#     for line in f:
#         if 'GMES: ' in line:
#             measurementlines.append(line)
# # splits each line at a space and creates an array of arrays of strings
#     for i in range(len(measurementlines)):
#         measurementlines[i] = measurementlines[i].split(' ')
#
#         simpplaqs = []
#     for i in range(len(measurementlines)):
#         simpplaqs.append((float(measurementlines[i][4]))/3 + 1)
#
#     mean = sum(simpplaqs)/len(simpplaqs)
#     sd = np.std(simpplaqs)
#
#     themeans.append(mean)
#     thestandarddevs.append(sd)
#     data.append([betaone, betatwo, mean])
#     dataformatrix.append([int(float(file1[0]) * 10), int(float(file2[0]) * 10) + 17, mean])

a = np.empty([35, 61])
for x in dataformatrix:
    column = x[0]
    row = x[1]
    # print(row)
    # print(column)
    # print(row)
    # print(column)
    value = x[2]
    # print(row)
    # print(column)
    a[row, column] = value

# print(a)
# print(len(a))
y = 3 * 0.1 * np.arange(-17.0, 18.0)
x = 3 * 0.1 * np.arange(0.0, 61.0)

def getpeaklocations(pine, b2list):
    peaklocations = []
    for b2 in b2list:
        a = np.empty([2, 61])
        i = 0
        for x in pine:
            bb = x[1]
            # print(bb)
            betatwoforarray = (bb - 17) * 0.1
            # print(betatwoforarray)
            if betatwoforarray == b2:
                column = x[0]
                value = x[2]
                # print('yes')
            # print(row)
            # print(column)
                a[0, column] = value
                a[1, column] = x[3]
                i = i + 1
            # print(i)
        delEs = []
        betaones = []
        errors = []
        E = []
        Eerr = []
        for i in range(i - 1):
            if i > 0:
                delE = abs(a[0, i + 1] - a[0, i - 1])
                betaone = i * 0.1
                error = ((a[1, i + 1])**2 + (a[1, i - 1])**2)**0.5
                delEs.append(delE)
                betaones.append(betaone)
                errors.append(error)
                    # E.append(a[0, i])
                    # Eerr.append(a[1, i])
        # print(delEs)
        peaks = scipy.signal.find_peaks(delEs, height = 0.3)
            # print(peaks)
        for x in peaks[0]:
            peaklocations.append([betaones[x], b2, delEs[x], errors[x]])
    return peaklocations

b2 = 0.1 * np.arange(-34, 4)
peaklocations = getpeaklocations(dataformatrix, b2)
# print(peaklocations)
betaonespeak = []
betatwospeak = []
peakerrors = []
betaonespeakforlin = []
betatwospeakforlin = []
peakerrorsforlin = []
for x in peaklocations:
    betaonespeak.append(3 * x[0])
    betatwospeak.append(3 * x[1])
    peakerrors.append(x[3])
    if (3 * x[0]) >= 5:
        betaonespeakforlin.append(3 * x[0])
        betatwospeakforlin.append(3 * x[1])
        peakerrorsforlin.append(x[3])
x = 3 * 0.1 * np.arange(0.0, 61.0)
a = np.polyfit(betaonespeakforlin, betatwospeakforlin, 1, w = peakerrorsforlin)
y = a[0] * x
theslope = a[0]
# print(theslope)
theslopemult = a[0]*18
# print(theslopemult)
doodle = 11.0 * a[0]
# print(a)

abc = np.polyfit(betaonespeak, betatwospeak, 4, w = peakerrors)
poodle = abc[0]*(11**4) + abc[1]*(11**3) + abc[2]*(11**2) + abc[3]*(11) + abc[4]
yack = poodle - doodle
xyz = abc[0]*(x**4) + abc[1]*(x**3) + abc[2]*(x**2) + abc[3]*(x) + abc[4] - yack
moodle = abc[0]*(18**4) + abc[1]*(18**3) + abc[2]*(18**2) + abc[3]*(18) + abc[4]
# print(moodle)

betatwosfordiststuff = 0.1 * np.arange(-17.0, 0.0)
# tea = 3 * 0.1 * np.arange(0.0, 61.0)
Eult = []
Dist = []
Eulterr = []
for a in betatwosfordiststuff:
    E = []
    Eerr = []
    bwon = []

    for stuff in data:
        if round(stuff[1], 1) == round(3 * a, 1):
            E.append(stuff[2])
            Eerr.append(stuff[3])
            bwon.append(stuff[0])
            Eult.append([stuff[0], stuff[1], stuff[2], stuff[3]])
            # beta 1, beta 2, E, Eerr
            # print(a)
    if len(E) > 0:
        # print(len(tea))
        # print(len(E))
        plt.errorbar(bwon, E, Eerr, marker = 'o', markersize = '3', label = r'$\beta_2$ = ' + str(round(3 * a, 1)), linestyle = 'None')
        plt.show()
# E = []
# Eerr = []
# bwon = []
# for stuff in data:
#     if stuff[1] < -1.7 and stuff[1] > -2.0:
#         E.append(stuff[2])
#         Eerr.append(stuff[3])
#         bwon.append(stuff[0])
#         Eult.append([stuff[0], stuff[1], stuff[2], stuff[3]])
#         # Eulterr.append(stuff[3])
#         # Dist.append((stuff[0]**2 + stuff[1]**2)**0.5)
#         # print(stuff[1])
# plt.errorbar(bwon, E, Eerr, marker = 'o', markersize = 3, label = r'$\beta_2$ = -1.8', linestyle = 'None')
# plt.show()
#
# E = []
# Eerr = []
# bwon = []
# for stuff in data:
#     if stuff[1] < -0.8 and stuff[1] > -1.0:
#         E.append(stuff[2])
#         Eerr.append(stuff[3])
#         bwon.append(stuff[0])
#         Eult.append([stuff[0], stuff[1], stuff[2], stuff[3]])
#         # Eulterr.append(stuff[3])
#         # Dist.append((stuff[0]**2 + stuff[1]**2)**0.5)
#         # print(stuff[1])
# plt.figure()
# plt.errorbar(bwon, E, Eerr, marker = 'o', markersize = 3, label = r'$\beta_2$ = -0.9', linestyle = 'None')
# plt.show()
#
# E = []
# Eerr = []
# bwon = []
# for stuff in data:
#     if stuff[1] < -2.0 and stuff[1] > -2.2:
#         E.append(stuff[2])
#         Eerr.append(stuff[3])
#         bwon.append(stuff[0])
#         Eult.append([stuff[0], stuff[1], stuff[2], stuff[3]])
#         # Eulterr.append(stuff[3])
#         # Dist.append((stuff[0]**2 + stuff[1]**2)**0.5)
#         # print(stuff[1])
plt.errorbar(bwon, E, Eerr, marker = 'o', markersize = 3, label = r'$\beta_2$ = -2.1', linestyle = 'None')
plt.show()
plt.xlabel(r'$\beta_1$')
plt.ylabel('<E>')
plt.legend(loc = 'upper right')
plt.savefig('EvsB1forB2online.pdf')
plt.figure().clear()
plt.clf()
Eultultcurve = []
Distcurve = []
Eultulterrcurve = []
Eultult = []
Dist = []
Eultulterr = []
EultultWil = []
DistWil = []
EultulterrWil = []
lineB1 = []
curveB1 = []
WilB1 = []
for stuff in data:
    if round(stuff[0]*theslope, 1) == round(stuff[1], 1):
        if stuff[1] < 0 and stuff[1] > 18 * theslope:
            Eultult.append(stuff[2])
            Dist.append((stuff[0]**2 + stuff[1]**2)**0.5)
            Eultulterr.append(stuff[3])
            lineB1.append(stuff[0])

    if round(abc[0]*(stuff[0]**4) + abc[1]*(stuff[0]**3) + abc[2]*(stuff[0]**2) + abc[3]*(stuff[0]) + abc[4] - yack, 1) == round(stuff[1], 1):
        if stuff[1] < abc[4] - yack and stuff[1] > abc[0]*(18**4) + abc[1]*(18**3) + abc[2]*(18**2) + abc[3]*(18) + abc[4] - yack:
            Eultultcurve.append(stuff[2])
            Distcurve.append((stuff[0]**2 + stuff[1]**2)**0.5)
            curveB1.append(stuff[0])
            Eultulterrcurve.append(stuff[3])
    if stuff[1] == 0:
        EultultWil.append(stuff[2])
        DistWil.append((stuff[0]**2 + stuff[1]**2)**0.5)
        EultulterrWil.append(stuff[3])
        WilB1.append(stuff[0])

# bwonders = 3 * betatwosfordiststuff/theslope
# Eultult = []
# Dist = []
# Eultulterr = []
# for yoda in bwonders:
#     for luke in Eult:
#         if round(luke[0], 1) == round(yoda, 1):
#             Eultult.append(luke[2])
#             Eultulterr.append(luke[3])
#             Dist.append((luke[0]**2 + luke[1]**2)**0.5)
#             print('yes')
themeans = []
thestandarddevs = []
beta1 = []
beta2 = []
data = []
dataformatrix = []
# files = os.listdir("/home/guest/dym_par_adj/hotstartdata")
files = os.listdir("/home/guest/dym_par_adj/su3data")
# print(files)
su3beta = []
su3energy = []
for file in files:
    f = open(os.path.join("/home/guest/dym_par_adj/su3data", file), 'r')
    for x in f:
        meaning = x.split(' ')
        # print(meaning)
        su3beta.append(float(meaning[0]))
        su3energy.append(float(meaning[1]))
        if float(meaning[0]) == 6:
            yoyo = float(meaning[1])
# print(su3beta)


print(yoyo)
plt.figure()
plt.axvspan(6, 20, alpha = 0.2, color = '#FED141')
plt.show()
plt.errorbar(np.arange(0, 21, 1), yoyo*np.ones(21), color = '#AF272F')
# plt.errorbar(curveB1, Eultultcurve, Eultulterrcurve, linestyle = 'None', marker = 'o', markersize = 2, label = r'$S(1080)$ and LW along curve', color = '#F68D2E')
# plt.errorbar(lineB1, Eultult, Eultulterr, linestyle = 'None', marker = 'o', markersize = 2, label = r'$S(1080)$ and LW along line', color = '#4C8C2B')
plt.errorbar(WilB1, EultultWil, EultulterrWil, linestyle = 'None', marker = 'o', markersize = 2, label = r'$S(1080)$ and Wilson', color = '#41B6E6')
plt.errorbar(su3beta, su3energy,  linestyle = 'None', marker = 'o', markersize = 0.7, label = r'$SU(3)$ and Wilson', color = '#000000')
plt.xlabel(r"$\beta_1$", fontsize = 14)
plt.xlim((0, 20))
plt.ylim((-0.2, 1.1))
plt.show()

# plt.title(r'Average Plaquette Energy with Different $\beta$ Combinations')
plt.ylabel(r'$a^4<E>$', fontsize = 14)
plt.legend(loc = 'upper right', fontsize = 12)
plt.savefig('EvsBdistalongcurveSU3andS1080.pdf', bbox_inches = 'tight')
