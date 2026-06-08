import numpy as np
import os
import matplotlib.pyplot as plt
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
files = os.listdir("/home/guest/dym_par_adj/forthefirstcontourplot")
# print(files)
for file in files:
    # f = open(os.path.join("/home/guest/dym_par_adj/hotstartdata", file), 'r')
    f = open(os.path.join("/home/guest/dym_par_adj/forthefirstcontourplot", file), 'r')
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
    if float(file2[0]) >= 0:
        dataformatrix.append([int(float(file1[0]) * 20), int(float(file2[0]) * 20) + 120, mean, sd])
# print(beta2)
# print(len(themeans))
# print(len(dataformatrix))

files = os.listdir("/home/guest/dym_par_adj/forthefirstcontourplotnegative")
# print(files)
for file in files:
    # f = open(os.path.join("/home/guest/dym_par_adj/hotstartdata", file), 'r')
    f = open(os.path.join("/home/guest/dym_par_adj/forthefirstcontourplotnegative", file), 'r')
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
    dataformatrix.append([int(float(file1[0]) * 20), int(float(file2[0]) * 20) + 120, mean, sd])

a = np.empty([241, 121])
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

# to create an ordered list of E's for a fixed beta2
def energyderivative(pine, b2):
    a = np.empty([2, 121])
    i = 0
    for x in pine:
        bb = x[1]
        # print(bb)
        betatwoforarray = (bb - 120) * 0.05
        # print(betatwoforarray)
        if betatwoforarray == b2:
            column = x[0]
            value = x[2]
            print('yes')
    # print(row)
    # print(column)
            a[0, column] = value
            a[1, column] = x[3]
            i = i + 1
    print(i)
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
            E.append(a[0, i])
            Eerr.append(a[1, i])

    plt.errorbar(betaones, E, Eerr, linestyle = 'None', marker = 'o', markersize = 4)
    plt.xlabel('beta1')
    plt.ylabel('<E>')
    plt.title('<E> with Beta2 = 0.3')
    plt.show()
    plt.savefig('b2_' + str(b2) + '_Evb1_test.png')

    # plt.errorbar(betaones, delEs, errors, linestyle = 'None', marker = 'o', markersize = 4)
    # plt.xlabel('beta1')
    # plt.ylabel('|delat <E>|')
    # plt.title('Absolute Value of Change in <E> with Beta2 = 0.3')
    # plt.show()
    # plt.savefig('b2_' + str(b2) + '_Evb1_test.png')

energyderivative(dataformatrix, 0.1)
