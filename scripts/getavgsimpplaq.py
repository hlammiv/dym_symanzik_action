import numpy as np
import os
import matplotlib.pyplot as plt
# opens file (the path will need to be changed according to workstation)
# the path to out_b. . . files in the workstation: /home/guest/dym_par_adj/hotstartdata
# on my computer: C:\\Users\\annac\\OneDrive\\Desktop\\SIST Project\\Computer Code\\data\\data
themeans = []
thestandarddevs = []
betas = []
# files = os.listdir("/home/guest/dym_par_adj/hotstartdata")
files = os.listdir("/home/guest/dym_par_adj/comparisons/beta2is-0.1")
# print(files)
for file in files:
    # f = open(os.path.join("/home/guest/dym_par_adj/hotstartdata", file), 'r')
    f = open(os.path.join("/home/guest/dym_par_adj/comparisons/beta2is-0.1", file), 'r')
    file = file.split('b')
    file = file[1].split('_')

    print(file)
    # adds only the lines of the log file that contain measurements
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
    betas.append(3 * float(file[0]))
        # print(mean)
        # print(sd)
        # print(file)

plt.errorbar(betas, themeans, thestandarddevs, linestyle = 'None', marker = 'o', markersize = 4, label = 'Beta2 = -0.1')

# plt.title('Freeze Transition for S1080, the Wilson Action, and a Cold Start')


# themeans = []
# thestandarddevs = []
# betas = []
# # files = os.listdir("/home/guest/dym_par_adj/hotstartdata")
# files = os.listdir("/home/guest/dym_par_adj/coldsmallbetaincr")
# # print(files)
# for file in files:
#     # f = open(os.path.join("/home/guest/dym_par_adj/hotstartdata", file), 'r')
#     f = open(os.path.join("/home/guest/dym_par_adj/coldsmallbetaincr", file), 'r')
#     file = file.split('b')
#     file = file[1].split('_')
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
#     simpplaqs = []
#     for i in range(len(measurementlines)):
#         simpplaqs.append(float(measurementlines[i][5]))
#
#     mean = sum(simpplaqs)/len(simpplaqs)
#     sd = np.std(simpplaqs)
#
#     themeans.append(mean)
#     thestandarddevs.append(sd)
#     betas.append(3 * float(file[0]))
#         # print(mean)
#         # print(sd)
#         # print(file)
#
# plt.errorbar(betas, themeans, thestandarddevs, linestyle = 'None', marker = 'o', markersize = 4, label = 'Cold Start')
# plt.xlabel('Beta1')
# plt.ylabel('Average Plaquette Value')
# plt.title('Freeze Transition for S1080, the Wilson Action, and a Hot and Cold Start')
# # plt.title('Freeze Transition for S1080, the Wilson Action, and a Cold Start')
# plt.legend(loc = 'upper right')
# plt.show()
# # plt.savefig("hotstartdata.png")
# plt.savefig("hotandcoldsmallbetaincr.png")

themeans = []
thestandarddevs = []
betas = []
# files = os.listdir("/home/guest/dym_par_adj/hotstartdata")
files = os.listdir("/home/guest/dym_par_adj/comparisons/beta2is0.1")
# print(files)
for file in files:
    # f = open(os.path.join("/home/guest/dym_par_adj/hotstartdata", file), 'r')
    f = open(os.path.join("/home/guest/dym_par_adj/comparisons/beta2is0.1", file), 'r')
    file = file.split('b')
    file = file[1].split('_')

    print(file)
    # adds only the lines of the log file that contain measurements
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
    betas.append(3 * float(file[0]))
        # print(mean)
        # print(sd)
        # print(file)

plt.errorbar(betas, themeans, thestandarddevs, linestyle = 'None', marker = 'o', markersize = 4, label = 'Beta2 = 0.1')

themeans = []
thestandarddevs = []
betas = []
# files = os.listdir("/home/guest/dym_par_adj/hotstartdata")
files = os.listdir("/home/guest/dym_par_adj/comparisons/beta2is0")
# print(files)
for file in files:
    # f = open(os.path.join("/home/guest/dym_par_adj/hotstartdata", file), 'r')
    f = open(os.path.join("/home/guest/dym_par_adj/comparisons/beta2is0", file), 'r')
    file = file.split('b')
    file = file[1].split('_')

    print(file)
    # adds only the lines of the log file that contain measurements
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
    betas.append(3 * float(file[0]))
        # print(mean)
        # print(sd)
        # print(file)

plt.errorbar(betas, themeans, thestandarddevs, linestyle = 'None', marker = 'o', markersize = 4, label = 'Beta2 = 0.0')


themeans = []
thestandarddevs = []
betas = []
# files = os.listdir("/home/guest/dym_par_adj/hotstartdata")
files = os.listdir("/home/guest/dym_par_adj/comparisons/beta2is-0.2")
# print(files)
for file in files:
    # f = open(os.path.join("/home/guest/dym_par_adj/hotstartdata", file), 'r')
    f = open(os.path.join("/home/guest/dym_par_adj/comparisons/beta2is-0.2", file), 'r')
    file = file.split('b')
    file = file[1].split('_')

    print(file)
    # adds only the lines of the log file that contain measurements
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
    betas.append(3 * float(file[0]))
        # print(mean)
        # print(sd)
        # print(file)

plt.errorbar(betas, themeans, thestandarddevs, linestyle = 'None', marker = 'o', markersize = 4, label = 'Beta2 = -0.2')

themeans = []
thestandarddevs = []
betas = []
# files = os.listdir("/home/guest/dym_par_adj/hotstartdata")
files = os.listdir("/home/guest/dym_par_adj/comparisons/beta2is0.2")
# print(files)
for file in files:
    # f = open(os.path.join("/home/guest/dym_par_adj/hotstartdata", file), 'r')
    f = open(os.path.join("/home/guest/dym_par_adj/comparisons/beta2is0.2", file), 'r')
    file = file.split('b')
    file = file[1].split('_')

    print(file)
    # adds only the lines of the log file that contain measurements
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
    betas.append(3 * float(file[0]))
        # print(mean)
        # print(sd)
        # print(file)

plt.errorbar(betas, themeans, thestandarddevs, linestyle = 'None', marker = 'o', markersize = 4, label = 'Beta2 = 0.2')

plt.xlabel('Beta1')
plt.ylabel('Average Plaquette Value')
plt.legend(loc = 'upper right')
plt.show()
plt.savefig("Betasnegative0.2to0.2.png")
