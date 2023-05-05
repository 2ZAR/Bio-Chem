#Practice#
#ELISA#

%matplotlib inline
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd

df = pd.read_csv('./data/4PL_standard.csv')
df = df.iloc[::-1] 
df

def logistic4(x, A, B, C, D):
    return ((A-D)/(1.0+((x/C)**B))) + D

xdata = df['Conc']
ydata = df['Value']
popt, pcov = curve_fit(logistic4, xdata, ydata)
print(popt)

import numpy as np
x_fit = np.linspace(0.01,50000,50000)
y_fit = logistic4(x_fit, *popt)
y_fit

plt.plot(xdata, ydata, 'o', label = 'data')
plt.plot(x_fit, y_fit, label = 'fit')

plt.plot(xdata, ydata, 'o', label = 'data')
plt.plot(x_fit, y_fit, label = 'fit')
plt.xscale('log')

def solvex(y, A, B, C, D):
    return C*(((A-D)/(y-D)-1.0)**(1.0/B))

solvex(1.0,*popt)