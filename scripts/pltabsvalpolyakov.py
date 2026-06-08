import numpy as np
import os
import matplotlib.pyplot as plt
# opens file (the path will need to be changed according to workstation)
# the path to out_b. . . files in the workstation: /home/guest/dym_par_adj/hotstartdata
# on my computer: C:\\Users\\annac\\OneDrive\\Desktop\\SIST Project\\Computer Code\\data\\data
betas = []
themeans = []
thestandarddevs = []
betaprimes = []
# files = os.listdir("/home/guest/dym_par_adj/hotN10000K1000")
files = os.listdir("/home/guest/dym_par_adj/Nt4Ns12datathird")
# print(files)
for file in files:
    absvalpolyakov = []
    # f = open(os.path.join("/home/guest/dym_par_adj/hotN10000K1000", file), 'r')
    f = open(os.path.join("/home/guest/dym_par_adj/Nt4Ns12datathird", file), 'r')
    file = file.split('b1')

    file = file[1].split('b2')

    file1 = file[0].split('_')
    file2 = file[1].split('_')

    betaone = 3 * float(file1[0])
    betatwo = 3 * float(file2[0])
    # print(betaone)
    # print(betatwo)
    # beta1.append(3 * float(file1[0]))
    # beta2.append(3 * float(file2[0]))
    # print(file)
    # adds only the lines of the log file that contain measurements
    measurementlines = []
    for line in f:
        if 'GMES: ' in line:
            measurementlines.append(line)

    # splits each line at a space and creates an array of arrays of strings
    for i in range(len(measurementlines)):
        measurementlines[i] = measurementlines[i].split(' ')

    for i in range(len(measurementlines)):
        # absval = (((float(measurementlines[i][2]))**2) + ((float(measurementlines[i][3]))**2))**0.5
        absval = (float(measurementlines[i][2]))
        absvalpolyakov.append(absval)

    mean = abs(sum(absvalpolyakov)/len(absvalpolyakov))
    print(len(absvalpolyakov))
    sd = np.std(absvalpolyakov)
    themeans.append(mean)
    thestandarddevs.append(sd)
    betaprimes.append(((3 * float(file1[0]))**2 + (3 * float(file2[0]))**2)**0.5)
        # print(mean)
        # print(sd)
        # print(file)

plt.errorbar(betaprimes, themeans, thestandarddevs, linestyle = 'None', marker = 'o', markersize = 4, label = 'Nt = 4; Ns = 12', color = '#F68D2E')
plt.show()

# betas = []
# themeans = []
# thestandarddevs = []
# betaprimes = []
# # files = os.listdir("/home/guest/dym_par_adj/hotN10000K1000")
# files = os.listdir("/home/guest/dym_par_adj/Nt6Ns18data")
# # print(files)
# for file in files:
#     absvalpolyakov = []
#     # f = open(os.path.join("/home/guest/dym_par_adj/hotN10000K1000", file), 'r')
#     f = open(os.path.join("/home/guest/dym_par_adj/Nt6Ns18data", file), 'r')
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
#     # beta1.append(3 * float(file1[0]))
#     # beta2.append(3 * float(file2[0]))
#     # print(file)
#     # adds only the lines of the log file that contain measurements
#     measurementlines = []
#     for line in f:
#         if 'GMES: ' in line:
#             measurementlines.append(line)
#
#     # splits each line at a space and creates an array of arrays of strings
#     for i in range(len(measurementlines)):
#         measurementlines[i] = measurementlines[i].split(' ')
#
#     for i in range(len(measurementlines)):
#         absval = (((float(measurementlines[i][2]))**2) + ((float(measurementlines[i][3]))**2))**0.5
#         absvalpolyakov.append(absval)
#
#     mean = sum(absvalpolyakov)/len(absvalpolyakov)
#     print(len(absvalpolyakov))
#     sd = np.std(absvalpolyakov)
#     themeans.append(mean)
#     thestandarddevs.append(sd)
#     betaprimes.append(((3 * float(file1[0]))**2 + (3 * float(file2[0]))**2)**0.5)
#         # print(mean)
#         # print(sd)
#         # print(file)
#
# plt.errorbar(betaprimes, themeans, thestandarddevs, linestyle = 'None', marker = 'o', markersize = 4, label = 'Nt = 6; Ns = 18', color = '#4C8C2B')
# plt.show()
#
# betas = []
# themeans = []
# thestandarddevs = []
# betaprimes = []
# # files = os.listdir("/home/guest/dym_par_adj/hotN10000K1000")
# files = os.listdir("/home/guest/dym_par_adj/Nt8Ns24data")
# # print(files)
# for file in files:
#     absvalpolyakov = []
#     # f = open(os.path.join("/home/guest/dym_par_adj/hotN10000K1000", file), 'r')
#     f = open(os.path.join("/home/guest/dym_par_adj/Nt8Ns24data", file), 'r')
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
#     # beta1.append(3 * float(file1[0]))
#     # beta2.append(3 * float(file2[0]))
#     # print(file)
#     # adds only the lines of the log file that contain measurements
#     measurementlines = []
#     for line in f:
#         if 'GMES: ' in line:
#             measurementlines.append(line)
#
#     # splits each line at a space and creates an array of arrays of strings
#     for i in range(len(measurementlines)):
#         measurementlines[i] = measurementlines[i].split(' ')
#
#     for i in range(len(measurementlines)):
#         absval = (((float(measurementlines[i][2]))**2) + ((float(measurementlines[i][3]))**2))**0.5
#         absvalpolyakov.append(absval)
#
#     mean = sum(absvalpolyakov)/len(absvalpolyakov)
#     print(len(absvalpolyakov))
#     sd = np.std(absvalpolyakov)
#     themeans.append(mean)
#     thestandarddevs.append(sd)
#     betaprimes.append(((3 * float(file1[0]))**2 + (3 * float(file2[0]))**2)**0.5)
#         # print(mean)
#         # print(sd)
#         # print(file)
#
# plt.errorbar(betaprimes, themeans, thestandarddevs, linestyle = 'None', marker = 'o', markersize = 4, label = 'Nt = 8; Ns = 24', color = '#41B6E6')
# plt.show()


plt.xlabel(r"$\beta$'", fontsize = 14)
plt.ylabel('Abs Val Real Part Polyakov', fontsize = 14)
plt.legend(loc = 'upper left', fontsize = 12)
plt.savefig('PolyakovvBetaDist.pdf', bbox_inches = 'tight')
