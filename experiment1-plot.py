#! -*- coding: utf-8 -*-
import csv
import matplotlib.pyplot as plt

x = []
y1 = []
y2 = []
with open('experiment1.csv', 'r') as f:
    reader = csv.reader(f)
    # header = next(reader)

    for row in reader:
        x.append(int(row[0]))
        y1.append(float(row[1]))
        y2.append(float(row[2]))
plt.plot(x, y1, '+r', label="random")
plt.plot(x, y2, 'ob', label="BO")
plt.xlabel("experiment")
plt.ylabel("best score so far")
plt.legend(loc="lower right")
plt.show()
