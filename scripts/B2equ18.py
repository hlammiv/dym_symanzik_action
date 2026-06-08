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
bwon = []
# files = os.listdir("/home/guest/dym_par_adj/hotstartdata")
files = os.listdir("/home/guest/dym_par_adj/KNtendata")
# print(files)
for file in files:
    # f = open(os.path.join("/home/guest/dym_par_adj/hotstartdata", file), 'r')
    f = open(os.path.join("/home/guest/dym_par_adj/KNtendata", file), 'r')
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
    bwon.append(betatwo)
    dataformatrix.append([int(float(file1[0]) * 10), int(float(file2[0]) * 10) + 17, mean, sd])
plt.errorbar(bwon, themeans, thestandarddevs, linestyle = 'None', marker = 'o', markersize = 3)
plt.title(r'<E> vs. $\beta_2$ for $\beta_1$ = 18')
plt.ylabel('<E>')
plt.xlabel(r'$\beta_2$')

themeans = []
thestandarddevs = []
beta1 = []
beta2 = []
data = []
dataformatrix = []
bwon = []
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
    if round(betaone, 1) == 18:

        themeans.append(mean)
        thestandarddevs.append(sd)
        data.append([betaone, betatwo, mean, sd])
        bwon.append(betatwo)
        dataformatrix.append([int(float(file1[0]) * 10), int(float(file2[0]) * 10) + 17, mean, sd])
plt.errorbar(bwon, themeans, thestandarddevs, linestyle = 'None', marker = 'o', markersize = 3)
plt.savefig('B1is18graph.pdf')
