import numpy as np
import os
import matplotlib.pyplot as plt
# opens file (the path will need to be changed according to workstation)
# the path to out_b. . . files in the workstation: /home/guest/dym_par_adj/hotstartdata
# on my computer: C:\\Users\\annac\\OneDrive\\Desktop\\SIST Project\\Computer Code\\data\\data
therealpt = []
theimgpt = []
betas = []
# files = os.listdir("/home/guest/dym_par_adj/hotN10000K1000")
files = os.listdir("/home/guest/dym_par_adj/coldK10004x4")
# print(files)
for file in files:
    # f = open(os.path.join("/home/guest/dym_par_adj/hotN10000K1000", file), 'r')
    f = open(os.path.join("/home/guest/dym_par_adj/coldK10004x4", file), 'r')
    file = file.split('b')
    file = file[1].split('_')
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
        therealpt.append(float(measurementlines[i][2]))
        theimgpt.append(float(measurementlines[i][3]))

    beta = 3 * float(file[0])
        # print(mean)
        # print(sd)
        # print(file)


    plt.errorbar(therealpt, theimgpt, linestyle = 'None', marker = 'o', markersize = 4, label = 'Beta = ' + str(beta))
    plt.xlabel('Real Pt Polyakov')
    plt.ylabel('Imaginary Pt Polyakov')
    plt.ylim((-0.6, 0.6))
    plt.xlim((-0.6, 0.6))
    plt.title('Freeze Transition for S1080, the Wilson Action, and a Cold Start')
    # plt.title('Freeze Transition for S1080, the Wilson Action, and a Hot Start')
    plt.legend(loc = 'upper right')
    plt.show()

    figtitle = 'realvimpolyakovCOLD4^4' + str(beta) + '.png'
    # figtitle = 'realvimpolyakovHOT' + str(beta) + '.png'
    plt.savefig('coldK10004^4graphs/' + figtitle)
    plt.figure().clear()
