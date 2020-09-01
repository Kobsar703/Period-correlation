#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 20 19:00:12 2019

@author: alex
"""

import numpy as np
import matplotlib.pyplot as plt

import tkinter as tk
import os
import math
from functools import partial
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from tkinter import filedialog
plt.rcParams.update({'font.size': 20})

def data(path):

    DATA = os.path.join(path)
    Data = np.loadtxt(DATA)
    x = np.array(Data[:, 0])
    y = np.array(Data[:, 1])
    y_err = np.array(Data[:, 2])
    return x, y, y_err


file_photo = ''


def callback():
    global file_photo
    name = filedialog.askopenfilename()
    file_photo = str(name)
    return file_photo


file_magnetic = ''


def callback1():
    global file_magnetic
    name = filedialog.askopenfilename()
    file_magnetic = str(name)
    return file_magnetic


def find_nearest(array, value):
    array = np.asarray(array)
    idex = (np.abs(array - value)).argmin()
    return array[idex]


def coef_cor(ph, ma):
    Pearson = np.sum((ph - np.sum(ph) / len(ph)) * (ma - np.sum(ma) / len(ma))) / math.pow(
        np.sum((ph - np.sum(ph) / len(ph))**2) * np.sum((ma - np.sum(ma) / len(ma))**2), 0.5)
    return Pearson


def amplitude(period):

    global file_photo
    global file_magnetic

    x, y, y_err = data(file_photo)
    x1, y1, y_err1 = data(file_magnetic)

    x_per = [((i - x[np.argmin(y)]) / period) % 1 for i in x]
    x_per1 = [((i - (2450000 + x[np.argmin(y)])) / period) % 1 for i in x1]
    y_per = [y[np.where(x_per == find_nearest(x_per, i))][0] for i in x_per1]

    Pearson = coef_cor(y_per, y1)

    a = 2.0 * np.sum(np.double(np.cos(2.0 * x * np.pi / period) * y)) / len(y)
    b = 2.0 * np.sum(np.double(np.sin(2.0 * x * np.pi / period) * y)) / len(y)

    c = 2.0 * \
        np.sum(np.double(np.cos(2.0 * x1 * np.pi / period) * y1)) / len(y1)
    d = 2.0 * \
        np.sum(np.double(np.sin(2.0 * x1 * np.pi / period) * y1)) / len(y1)

    a_omg = np.double(math.pow(a**2 + b**2, 0.5))
    ampl_mag = np.double(math.pow(c**2 + d**2, 0.5))
    phase = np.double(math.atan(b / a))
    phase_mag = np.double(math.atan(d / c))

    return a_omg, phase, Pearson, ampl_mag, phase_mag


def create_image(period):

    global file_photo
    global file_magnetic

    x, y, y_err = data(file_photo)
    x1, y1, y_err1 = data(file_magnetic)

    x_per = [((i - x[np.argmax(y)]) / period) % 1 for i in x]
    x_per1 = [((i - (2450000 + x[np.argmax(y)])) / period) % 1 for i in x1]
    y_per = [y[np.where(x_per == find_nearest(x_per, i))][0] for i in x_per1]

    Pearson = coef_cor(y_per, y1)

    fig = plt.figure()
    fig.set_dpi(100)

    plt.subplot(2, 1, 1)

    plt.plot(x_per, y, '.', color='green', label='Period')
    plt.ylim([max(y), min(y)])
    plt.xlim([0, 1])
    plt.ylabel('Aplitude, mmag')
    plt.title(
        f'Period = {round(period, 5)}, Pearson coefficient = {round(Pearson, 3)}')
    plt.subplot(2, 1, 2)
    plt.plot(x_per1, y1, '.', color='blue')
    plt.grid(True)
    plt.ylim([max(y1), min(y1)])
    plt.xlim([0, 1])
    plt.xlabel('time')
    plt.ylabel('$B_{z}$, G')
    plt.savefig(f'Best_correlation{round(period, 5)}.png')

    return x_per, y, x_per1, y1, Pearson


def period_new(l, r, p):
    period1 = (l.get())
    period2 = (r.get())
    points = (p.get())

    period_try = np.linspace(period1, period2, points)
    Amplitude, phase, Pearson, Amplitude_mag, Phase_mag = [], [], [], [], []
    for i in period_try:
        Amp, ph, Pear, ampl_m, ph_m = amplitude(i)
        Amplitude.append(Amp)
        phase.append(ph)
        Pearson.append(Pear)
        Amplitude_mag.append(ampl_m)
        Phase_mag.append(ph_m)
    answer = {Amplitude[i]: period_try[i] for i in range(len(Amplitude))}
    answer_mag = {Amplitude_mag[i]: period_try[i]
                  for i in range(len(Amplitude_mag))}
    Searched_P = answer.get(max(Amplitude))
    Searched_P_mag = answer_mag.get(max(Amplitude_mag))

    cof = {Pearson[i]: period_try[i] for i in range(len(Pearson))}
    Min_coef = cof.get(min(Pearson))
    Max_coef = cof.get(max(Pearson))

    x_per4, y4, x_per5, y5, Pearson4 = create_image(Min_coef)
    x_per2, y2, x_per3, y3, Pearson1 = create_image(Max_coef)

    window = tk.Tk()

    window.title('Correlation')
    figure = plt.Figure(figsize=(6, 5), dpi=100)
    ax1 = figure.add_subplot(111)
    bar1 = FigureCanvasTkAgg(figure, window)
    bar1.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH)

    ax1.plot(period_try, Pearson, '.', color='blue')
    ax1.grid(True)
    ax1.set_xlabel('Period, days')
    ax1.set_ylabel('Pearson coef')

    toolbar = NavigationToolbar2TkAgg(bar1, window)
    toolbar.update()
    bar1._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)

    window1 = tk.Tk()

    window1.title('Find maximum Amplitude (Fourier transformation)')
    figure1 = plt.Figure(figsize=(12, 7), dpi=100)
    ax = figure1.add_subplot(111)
    ax1 = figure1.add_subplot(211)
    ax2 = figure1.add_subplot(212)
    bar1 = FigureCanvasTkAgg(figure1, window1)
    bar1.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH)

    ax1.plot(period_try, Amplitude, '.', color='green', label='Photometry')
    ax2.plot(period_try, Amplitude_mag, '.',
             color='blue', label='Magnetic field')
    ax1.grid(True)
    ax2.grid(True)
    ax.set_xticks([])
    ax.set_yticks([])
    ax2.set_xlabel('Period, days')
    ax2.set_ylabel('$B_{z}$, G')
    ax1.set_ylabel('Amplitude, mmag')
    ax.set_title('Photometry - TOP, Magnetic field - BOTTOM')

    toolbar = NavigationToolbar2TkAgg(bar1, window1)
    toolbar.update()
    bar1._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)

    window2 = tk.Tk()

    window2.title('Best anti-correlation')
    figure2 = plt.Figure(figsize=(12, 7), dpi=100)
    ax = figure2.add_subplot(111)
    ax1 = figure2.add_subplot(211)
    ax2 = figure2.add_subplot(212)
    bar1 = FigureCanvasTkAgg(figure2, window2)
    bar1.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH)

    ax1.plot(x_per4, y4, '.', color='green')
    ax2.plot(x_per5, y5, '.', color='blue')
    ax1.grid(True)
    ax2.grid(True)
    ax.set_xticks([])
    ax.set_yticks([])
    ax2.set_xlabel('Phase')
    ax2.set_ylabel('$B_{z}$, G')
    ax1.set_ylabel('Amplitude, mmag')
    ax.set_title('Photometry - TOP, Magnetic field - BOTTOM \n' +
                 f'Period = {round(Min_coef, 5)}, Pearson cof = {round(Pearson4, 3)}')

    toolbar = NavigationToolbar2TkAgg(bar1, window2)
    toolbar.update()
    bar1._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)

    window3 = tk.Tk()

    window3.title('Best correlation')
    figure3 = plt.Figure(figsize=(12, 7), dpi=100)
    ax = figure3.add_subplot(111)
    ax1 = figure3.add_subplot(211)
    ax2 = figure3.add_subplot(212)
    bar1 = FigureCanvasTkAgg(figure3, window3)
    bar1.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH)

    ax1.plot(x_per2, y2, '.', color='green')
    ax2.plot(x_per3, y3, '.', color='blue')
    ax1.grid(True)
    ax2.grid(True)
    ax.set_xticks([])
    ax.set_yticks([])
    ax2.set_xlabel('Phase')
    ax2.set_ylabel('$B_{z}$, G')
    ax1.set_ylabel('Amplitude, mmag')
    ax.set_title('Photometry - TOP, Magnetic field - BOTTOM \n' +
                 f'Period = {round(Max_coef, 5)}, Pearson cof = {round(Pearson1, 3)}')

    toolbar = NavigationToolbar2TkAgg(bar1, window3)
    toolbar.update()
    bar1._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)

    Label7 = tk.Label(
        root, text='Searched Period from photometry', bg='darkcyan')
    Label7.grid(row=6, column=0)
    Label8 = tk.Label(root, text=str(
        round(Searched_P, 5)) + ' days', bg='darkcyan')
    Label8.grid(row=6, column=1)
    Label9 = tk.Label(
        root, text='Maximum Amplitude from photometry', bg='darkcyan')
    Label9.grid(row=7, column=0)
    Label10 = tk.Label(root, text=str(
        round(max(Amplitude), 2)) + ' mmag', bg='darkcyan')
    Label10.grid(row=7, column=1)

    Label11 = tk.Label(
        root, text='Searched Period of magnetic field', bg='darkcyan')
    Label11.grid(row=8, column=0)
    Label12 = tk.Label(root, text=str(
        round(Searched_P_mag, 5)) + ' days', bg='darkcyan')
    Label12.grid(row=8, column=1)
    Label13 = tk.Label(
        root, text='Maximum Amplitude of magnetic field', bg='darkcyan')
    Label13.grid(row=9, column=0)
    Label14 = tk.Label(root, text=str(
        round(max(Amplitude_mag), 2)) + ' G', bg='darkcyan')
    Label14.grid(row=9, column=1)

    window3.mainloop()
    window2.mainloop()
    window1.mainloop()
    window.mainloop()

    return


root = tk.Tk()
root.title("Try correlation")
root.configure(background='darkcyan')

photo = 'Open file with photometry'
magnetic = 'Open file with magnetic field'
left_limit = 'Left limit of period'
right_limit = 'Right limit of periof'
try_period = 'Number of points between'


Label2 = tk.Label(root, justify=tk.LEFT, text=photo, padx=50, bg='darkcyan')
Label2.grid(row=0, column=0)
Label3 = tk.Label(root, justify=tk.LEFT, text=magnetic, padx=50, bg='darkcyan')
Label3.grid(row=1, column=0)
Label4 = tk.Label(root, justify=tk.LEFT, text=left_limit,
                  padx=50, bg='darkcyan')
Label4.grid(row=2, column=0)
Label5 = tk.Label(root, justify=tk.LEFT, text=right_limit,
                  padx=50, bg='darkcyan')
Label5.grid(row=3, column=0)
Label6 = tk.Label(root, justify=tk.LEFT, text=try_period,
                  padx=50, bg='darkcyan')
Label6.grid(row=4, column=0)


left_limit1 = tk.DoubleVar()
Left_limit = tk.Entry(root, textvariable=left_limit1)
Left_limit.grid(row=2, column=1, columnspan=2)

right_limit1 = tk.DoubleVar()
Right_limit = tk.Entry(root, textvariable=right_limit1)
Right_limit.grid(row=3, column=1, columnspan=2)

try_period1 = tk.IntVar()
Try_period = tk.Entry(root, textvariable=try_period1)
Try_period.grid(row=4, column=1, columnspan=2)

period_new = partial(period_new, left_limit1, right_limit1, try_period1)

button1 = tk.Button(root, text='Open File', bg='darkred',
                    fg='white', padx=50, command=callback)
button1.grid(row=0, column=1, columnspan=2)

button2 = tk.Button(root, text='Open File', bg='darkred',
                    fg='white', padx=50, command=callback1)
button2.grid(row=1, column=1, columnspan=2)

button = tk.Button(root, text='Give me Period, please',
                   bg='darkred', fg='white', padx=170, command=period_new)
button.grid(row=5, column=0, columnspan=2)

root.mainloop()
