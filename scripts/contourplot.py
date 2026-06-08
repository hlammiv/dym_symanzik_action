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
    data.append([betaone, betatwo, mean])
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
# print(tool)
# for x in range(len(data)):
#     for m in tool:
#         for n in tool:
#             if data[x][0] == m and data[x][2] == n:


# X, Y = np.meshgrid(beta1, beta2)
# fig, ax = plt.subplots()
# CS = ax.contour(X, Y, themeans)
# fig = plt.figure()
# fig.set_size_inches(14,8)

CS = plt.contourf(x, y, a, 20)
plt.contourf(x, y, a, 20)
plt.xlabel(r'$\beta_1$', fontsize = 14)
plt.ylabel(r'$\beta_2$', fontsize = 14)
proxy = [plt.Rectangle((1, 1), 2, 2, fc=pc.get_facecolor()[0]) for pc in CS.collections]
lgd = plt.legend(proxy, ['0.0 - 0.06', '0.06 - 0.12', '0.12 - 0.18', '0.18 - 0.24', '0.24 - 0.30', '0.30 - 0.36', '0.36 - 0.42', '0.42 - 0.48', '0.48 - 0.54', '0.54 - 0.60', '0.60 - 0.66', '0.66 - 0.72', '0.72 - 0.78', '0.78 - 0.84', '0.84 - 0.90', '0.90 - 0.96', '0.96 +'],  bbox_to_anchor = (1.3, 1), title = r'$a^4 <E>$')
plt.show()

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
plt.errorbar(x, y, color = 'w', linewidth = 4)
print(a[0])
plt.show()
doodle = 11.0 * a[0]
# print(a)

abc = np.polyfit(betaonespeak, betatwospeak, 4, w = peakerrors)
slopeslope = (betatwospeak[-1] - betatwospeak[-2])/(betaonespeak[-1] - betaonespeak[-2])
poodle = abc[0]*(11**4) + abc[1]*(11**3) + abc[2]*(11**2) + abc[3]*(11) + abc[4]
yack = poodle - doodle
xyz = abc[0]*(x**4) + abc[1]*(x**3) + abc[2]*(x**2) + abc[3]*(x) + abc[4] - yack
print(yack)

point = abc[0]*(18**4) + abc[1]*(18**3) + abc[2]*(18**2) + abc[3]*(18) + abc[4] - yack

lmnop = slopeslope * (x - 18) + point

slopeslopeslope = ((a[0] * 17.5) - (-0.6))/17.5

cookie = slopeslopeslope * x - 0.6

plt.errorbar(x, xyz, color = '#F68D2E', linewidth = 4)
print(slopeslopeslope)
# plt.errorbar(x, lmnop, color = 'r', linewidth = 4)
plt.errorbar(x, cookie, color = 'r', linewidth = 4)
plt.show()

# plt.errorbar(betaonespeak, betatwospeak, peakerrors, linestyle = 'None', marker = 'o', markersize = 4, mfc = 'white', mec = 'None', ecolor = 'white')
# # plt.xlim([0, 17.6])
# plt.show()
# plt.savefig('contourplotwithpeaks.pdf')

# CS = plt.contour(x, y, a, 20)
# plt.clabel(CS)
plt.ylim(top = 2)
plt.text(6, 0.8, r'$\beta_2$ = '+ str("{:.2e}".format(abc[0])) + r'$\beta_1^4$ + ' + str("{:.2e}".format(abc[1])) + r'$\beta_1^3$ + ' + str('\n') + str("{:.2e}".format(abc[2])) + r'$\beta_1^2$ + ' + str("{:.2e}".format(abc[3])) + r'$\beta_1$ + 1.23', c = '#F68D2E', fontsize = 12)
plt.text(7.5, 0, r'$\beta_2$ = ' + str("{:.2e}".format(a[0])) + r'$\beta_1$', c = 'w', fontsize = 12)
plt.savefig('contourplotforK1000N10000withdelEpeak.pdf', bbox_extra_artists=(lgd,), bbox_inches='tight')
