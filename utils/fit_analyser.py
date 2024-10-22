#!/usr/bin/python
###########################################################################
# Author: Rogerio Riffel                                                  #
# This is an adaptation of fit_scrutinizer written by Daniel Ruschel Duta #
# and available at: https://bitbucket.org/danielrd6/ifscube/             #
# IA - ChatGTP has been used to help in the development of the code      #
###########################################################################

import argparse
import tkinter as tk

# third party
import matplotlib
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk
from matplotlib.figure import Figure
import numpy as np
from numpy import ma
from mpl_toolkits.axes_grid1 import make_axes_locatable
import re
import matplotlib.pyplot as plt
# local
import sys
from astropy.io import fits as pf

class mclass:

    def __init__(self, window, megacube='None', PopBins='PoPBins', flxNorm='FLXOBS', flxSyn='FLXSYN', PoPVecsL='PoPVecsL', PoPVecsM='PoPVecsM'):
        self.cube_header = pf.getheader(megacube, PopBins)
        self.cube_data = pf.getdata(megacube, PopBins)
        self.l0 = pf.getheader(megacube, flxNorm)['CRVAL3']
        self.dl = pf.getheader(megacube, flxNorm)['CD3_3']

        self.data_cube_Lfrac = pf.getdata(megacube, PoPVecsL)
        self.data_cube_Mfrac = pf.getdata(megacube, PoPVecsM)
        self.data_cube_flux = pf.getdata(megacube, flxNorm)
        self.data_cube_synt = pf.getdata(megacube, flxSyn)
        (self.size, self.ypix, self.xpix) = np.shape(self.data_cube_flux)
        self.lamb = np.arange(self.l0, (self.dl * self.size + self.l0), self.dl)
        self.z, self.y, self.x = np.shape(self.cube_data)
        
        self.aux = pf.getdata(megacube, 'BaseAgeMetal')
        self.popage = self.aux[:, 0]
        self.popZ = self.aux[:, 1]
        self.Zs = np.unique(self.popZ)[np.unique(self.popZ) != 0]

        self.vecs = []
        for i in range(0, np.shape(self.cube_data)[0], 1):
            if 'FC' in self.cube_header['DATA' + str(i)]:
                self.vecs.append(i)
            elif 'BB' in self.cube_header['DATA' + str(i)]:
                self.vecs.append(i)
            elif 'light' in self.cube_header['DATA' + str(i)]:
                self.vecs.append(i)

        self.window = window
        self.buttonPlot = tk.Button(window, text='Image plot', command=self.plot)
        self.buttonContourPlot = tk.Button(window, text='Overplot contours', command=self.plot_contours)
        self.buttonFollow = tk.Button(window, text='Follow mouse', command=self.follow)
        self.buttonSinglePlot = tk.Button(window, text='Plot on click', command=self.singleplot)

        self.text = tk.Text(window, bg='white', height=10, width=80, font=('Fixedsys', 12))

        s_map = tk.Listbox(window, selectmode='single', exportselection=0, bg='yellow')
        for i in range(0, self.z, 1):
            s_map.insert('end', self.cube_header['DATA' + str(i)])

        self.fig = Figure(figsize=(6, 6))
        self.fitplot = Figure(figsize=(12, 3.5))
        self.hists = Figure(figsize=(12, 3.5))
        self.ax1 = self.fig.add_subplot(111)
        self.ax2 = self.fitplot.add_subplot(111)
        self.ax3 = self.hists.add_subplot(121)
        self.ax4 = self.hists.add_subplot(122)

        div = make_axes_locatable(self.ax1)
        self.cax = div.append_axes('right', size='5%', pad=0)
        self.cax.set_xticks([])
        self.cax.set_yticks([])

        canvas = FigureCanvasTkAgg(self.fig, master=self.window)
        canvas2 = FigureCanvasTkAgg(self.fitplot, master=self.window)
        canvas3 = FigureCanvasTkAgg(self.hists, master=self.window)
        canvas.get_tk_widget().grid(row=0, column=0, rowspan=6, columnspan=6)
        canvas2.get_tk_widget().grid(row=0, column=7, rowspan=3, columnspan=6)
        canvas3.get_tk_widget().grid(row=3, column=7, rowspan=3, columnspan=6)

        tb_frame = tk.Frame(window)
        tb_frame.grid(row=7, column=0, columnspan=6, sticky='W')
        toolbar = NavigationToolbar2Tk(canvas, tb_frame)
        toolbar.update()

        other_tb_frame = tk.Frame(window)
        other_tb_frame.grid(row=8, column=7, columnspan=12, sticky='W')
        other_toolbar = NavigationToolbar2Tk(canvas2, other_tb_frame)
        other_toolbar.update()

        self.buttonPlot.grid(row=8, column=0, sticky='W')
        self.buttonContourPlot.grid(row=9, column=0, sticky='W')
        self.buttonFollow.grid(row=10, column=0, sticky='W')
        self.buttonSinglePlot.grid(row=11, column=0, sticky='W')
        s_map.grid(row=8, column=1, rowspan=4)
        self.text.grid(row=9, column=7, rowspan=4)

        self.s_map = s_map

        # Initialize zmin and zmax to None for auto-scaling
        self.zmin = None
        self.zmax = None

        # Scale for rescaling
        self.rescale_var = tk.BooleanVar()
        self.rescale_checkbox = tk.Checkbutton(window, text="Rescale", variable=self.rescale_var, command=self.toggle_scales)
        self.rescale_checkbox.grid(row=3, column=6)

        self.scale_zmin = tk.Scale(window, from_=100, to=0, orient='vertical', command=self.update_zmin)
        self.scale_zmin.grid(row=2, column=6, sticky='NS')
        self.scale_zmin.config(state=tk.DISABLED)

        self.scale_zmax = tk.Scale(window, from_=100, to=0, orient='vertical', command=self.update_zmax)
        self.scale_zmax.grid(row=1, column=6, sticky='NS')
        self.scale_zmax.config(state=tk.DISABLED)

        # Scale for controlling number of contour levels
        self.levels_scale = tk.Scale(window, from_=3, to=20, orient='horizontal', label='Contour Levels',
                                     command=self.update_contour_levels)
        self.levels_scale.grid(row=5, column=6, sticky='WE')
        self.levels_scale.set(10)  # Default is 10 levels

        self.contour_levels = 5  # Initial number of contour levels

        self.canvas = canvas
        self.canvas2 = canvas2
        self.canvas3 = canvas3

    def toggle_scales(self):
        if self.rescale_var.get():
            # Enable scales if checkbox is checked
            self.scale_zmin.config(state=tk.NORMAL)
            self.scale_zmax.config(state=tk.NORMAL)
            self.scale_zmin.set(0)
            self.scale_zmax.set(100)
        else:
            # Disable scales if checkbox is unchecked and reset to auto-scale
            self.scale_zmin.config(state=tk.DISABLED)
            self.scale_zmax.config(state=tk.DISABLED)
            self.zmin = None
            self.zmax = None
            self.plot()  # Redraw with auto-scaling

    def update_contour_levels(self, value):
        self.contour_levels = int(value)
        self.plot_contours()

    def singleplot(self):
        try:
            self.canvas.mpl_disconnect(self.connect_id)
        except:
            pass

        self.connect_id = self.canvas.mpl_connect('button_press_event', self.onclick)

    def follow(self):
        try:
            self.canvas.mpl_disconnect(self.connect_id)
        except:
            pass

        self.connect_id = self.canvas.mpl_connect('motion_notify_event', self.onclick)

    def getim(self):
        # Check if any selection is made in the Listbox
        if not self.s_map.curselection():
            # If no selection, return None to indicate invalid state
            self.text.delete('1.0', 'end')
            self.text.insert('insert', 'No data selection made in the Listbox!\n')
            return None, None, None

        cube_data = self.cube_data
        cm = 'viridis'
        par = self.s_map.curselection()[0]
        d = self.cube_data[par, :, :]
        t = self.cube_header['DATA' + str(par)]

        return d, cm, t
        
    def plot(self, contrast=1):
        a = self.ax1
        a.cla()
        self.cax.cla()
        d, cm, t = self.getim()
        a.set_title(t)

        # Use zmin and zmax for scaling, or auto-scaling if they are None
        im = a.pcolormesh(d, cmap=cm, vmin=self.zmin if self.zmin is not None else d.min(),
                          vmax=self.zmax if self.zmax is not None else d.max())

        matplotlib.colorbar.Colorbar(self.cax, im)

        a.set_aspect('equal', 'datalim')

        self.canvas.draw()

    def plot_contours(self):
        # Get the data, colormap, and title
        d, _, _ = self.getim()

        # If getim() returned None, exit the function
        if d is None:
            return

        a = self.ax1

        # Generate contour data based on selected number of levels
        contour_levels = np.linspace(d.min(), d.max(), self.contour_levels)
        a.contour(d, levels=contour_levels, colors='red', linewidths=2)  # Contour over existing plot

        self.canvas.draw()

    def update_zmin(self, value):
        self.zmin = float(value)
        self.plot()

    def update_zmax(self, value):
        self.zmax = float(value)
        self.plot()

    def onclick(self, event):
        try:
            if event.xdata is not None and event.ydata is not None:
                i, j = [int(np.floor(x) + 0.5) for x in (event.xdata, event.ydata)]
                if np.any(np.array([i, j]) < 0):
                    self.text.delete('1.0', 'end')
                    self.text.insert('insert', 'Index Error!')
                    return
            else:
                self.text.delete('1.0', 'end')
                self.text.insert('insert', 'You clicked outside the plot!')
                return
        except Exception as e:
            self.text.delete('1.0', 'end')
            self.text.insert('insert', f'Error: {str(e)}')
            return

        self.text.delete('1.0', 'end')
        a = self.ax2
        a.cla()

        self.text.insert('insert', ('Best fit for Spaxel x= {:2} and y= {:2}\n').format(i, j))
        self.text.insert('insert', '{:10} {:10} {:10} {:10}\n'.format('  Age', '  Z', '     Lfrac', '   Mfrac'))

        try:
            self.flx = self.data_cube_flux[:, j, i]
            self.synt = self.data_cube_synt[:, j, i]
            self.lfrac = self.data_cube_Lfrac[:, j, i]
            self.mfrac = self.data_cube_Mfrac[:, j, i]

            a.plot(self.lamb, self.flx, color='blue',label='Obs')
            a.plot(self.lamb, self.synt, color='red',label='Synt')
            a.legend(loc=0, frameon=False, ncol=1, prop={'size': 12})
            a.set_xlabel('Wavelength', labelpad=1)
            a.set_ylabel('Flux')

            s = []
            for u in range(0, len(self.popage)):
                try:
                    s.append((self.popage[u], self.popZ[u]))
                except:
                    ' '
            for l in range(0, len(s)):
                self.text.insert('insert', '{:.2E}    '.format(s[l][0]))
                self.text.insert('insert', '{:10}'.format(s[l][1]))
                try:
                    self.text.insert('insert', '{:8.2f}'.format(self.lfrac[l]))
                    self.text.insert('insert', '{:8.2f}'.format(self.mfrac[l]))
                except:
                    ' '
                self.text.insert('insert', '\n')

        except IndexError:
            self.text.insert('insert', 'Index Error!')

        self.canvas2.draw()

        b = self.ax3
        b.cla()

        xlim = np.log10([1e5, 20e9])  # do not change)
        ylim = [0, 100]  # do not change

        b.set_xlim(xlim)
        b.set_ylim(ylim)
        b.set_xlabel('log(Age)', labelpad=1)
        b.set_ylabel('%')

        self.popx = self.data_cube_Lfrac[:, j, i]
        self.popm = self.data_cube_Mfrac[:, j, i]

        self.unique_ages = np.unique(self.popage[self.popage > 1.])
        self.summed_ages_light = np.zeros(len(self.unique_ages))
        self.summed_ages_mass = np.zeros(len(self.unique_ages))

        for z in self.Zs:
            for self.i in np.arange(0, len(self.unique_ages)):
                for self.k in np.arange(0, len(self.popage)):
                    if (self.popage[self.k] == self.unique_ages[self.i]) & (self.popZ[self.i] == z):
                        self.summed_ages_light[self.i] += self.popx[self.k]
                        self.summed_ages_mass[self.i] += self.popm[self.k]

        b.bar(np.log10(self.unique_ages), self.summed_ages_light, width=0.8, align='center', color='None', edgecolor='blue', label=r'$\Sigma x_j$')
        b.bar(np.log10(self.unique_ages), self.summed_ages_mass, width=0.2, align='center', color='None', edgecolor='red', label=r'$\Sigma\mu_j$', ls='dotted')
        b.legend(loc=0, frameon=False, ncol=1, prop={'size': 12})

        self.canvas3.draw()

        c = self.ax4
        c.cla()
        c.set_xlabel('Vecs',labelpad=1) 
        c.set_ylabel('%')
        c.set_ylim(0, 100)
        for self.k in self.vecs:
            c.bar(self.k, self.cube_data[self.k, j, i], alpha=0.5)
            c.text(self.k, 50, self.cube_header['DATA' + str(self.k)], rotation=90, fontsize=9, color='black')

        self.canvas3.draw()

if __name__ == '__main__':
    try:
        megacube = sys.argv[1]
        window = tk.Tk()
        window.title('Urutau Fit Analyser:' + megacube)
        start = mclass(window, megacube=megacube)
        window.mainloop()
    except:
        print("\n\
                  ***********************\n\
                  To run fit_analyser do: \n\
                  python fit_analyser urutau_output_cube_name.fits \n\
                   ***********************\n\
                  ")
