# coding=gbk

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
import Tkinter as tk
from tkFileDialog import askopenfilename
from mpl_toolkits.basemap import Basemap
from matplotlib import animation
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkMessageBox
from subprocess import call
import os


depthName = ''
boundryName = ''
numX = 0
numY = 0
numSteps = 120
resolution = 10
totalT = 44712
enableAdvection = 1
enableViscosity = 1
enableFriction = 1
enableCoriolisForce = 2
autoXStep = 1
viscosityCoefficient = 1e3
bottomFriction = 1.8e-3
alpha = 0.5
latitude = 34.
longitude = 117.5

eps = 1e-3
invalid = -1.e5

ratio = 7. / 6.

class GUI(tk.Tk):
    def __init__(self, *args, **kwargs):
        tk.Tk.__init__(self, *args, **kwargs)

        self.width = tk.IntVar()
        self.width.set(6)
        self.depthname = tk.StringVar()
        self.depthname.set(depthName)
        self.boundryname = tk.StringVar()
        self.boundryname.set(boundryName)

        self.numX = tk.IntVar()
        self.numX.set(numX)
        self.numY = tk.IntVar()
        self.numY.set(numY)
        self.numSteps = tk.IntVar()
        self.numSteps.set(numSteps)

        self.resolution = tk.DoubleVar()
        self.resolution.set(resolution)
        self.totalT = tk.DoubleVar()
        self.totalT.set(totalT)

        self.enableAdvection = tk.IntVar()
        self.enableAdvection.set(enableAdvection)
        self.enableViscosity = tk.IntVar()
        self.enableViscosity.set(enableViscosity)
        self.enableFriction = tk.IntVar()
        self.enableFriction.set(enableFriction)
        self.enableCoriolisForce = tk.IntVar()
        self.enableCoriolisForce.set(enableCoriolisForce)
        self.autoXStep = tk.IntVar()
        self.autoXStep.set(autoXStep)

        self.viscosityCoefficient = tk.DoubleVar()
        self.viscosityCoefficient.set(viscosityCoefficient)
        self.bottomFriction = tk.DoubleVar()
        self.bottomFriction.set(bottomFriction)
        self.alpha = tk.DoubleVar()
        self.alpha.set(alpha)
        self.latitude = tk.DoubleVar()
        self.latitude.set(latitude)
        self.longitude = tk.DoubleVar()
        self.longitude.set(longitude)

        self.index = 0
        self.anim = None
        self.title(u'二维潮波模拟')
        self.protocol("WM_DELETE_WINDOW", lambda: (self.quit(), self.destroy()))
        self.makeWidgets()
        return

    def makeWidgets(self):
        frame = tk.Frame(self)
        # tk.Button(frame, text='test', command=self.makeMap).pack(side=tk.TOP)
        self.makeNameEntry(frame)
        self.inputXYS(frame)
        self.makeSelect(frame)
        self.makeButton(frame)
        self.makeLogo(frame)
        # self.makeMap()
        frame.grid()
        return

    def makeLogo(self, frame):
        f = tk.Frame(frame)
        tk.Label(f, text=u"作者:黄龑 海洋科学2011级", foreground='blue', font=('YaHei', 11)).grid(row=0, column=0)
        f.grid()

    def makeNameEntry(self, frame):
        f = tk.Frame(frame)
        tk.Label(f, text=u'水深数据文件名').grid(row=0, column=0, sticky=tk.W)
        tk.Entry(f, width=9, textvariable=self.depthname).grid(row=0, column=1)
        tk.Button(f, text='...', command=self.getdepthname).grid(row=0, column=2)
        tk.Label(f, text=u'宽度', width=4).grid(row=0, column=3, sticky=tk.W)

        tk.Label(f, text=u'开边界数据文件名').grid(row=1, column=0, sticky=tk.W)
        tk.Entry(f, width=9, textvariable=self.boundryname).grid(row=1, column=1)
        tk.Button(f, text='...', command=self.getboundryname).grid(row=1, column=2)
        tk.Entry(f, width=3, textvariable=self.width).grid(row=1, column=3)
        f.grid()

    def inputXYS(self, frame):
        f = tk.Frame(frame)
        tk.Label(f, text=u'X方向网格数').grid(row=2, column=0)
        tk.Entry(f, width=5, textvariable=self.numX).grid(row=2, column=1)
        tk.Label(f, text=u'时间步数').grid(row=2, column=2)
        tk.Entry(f, width=5, textvariable=self.numSteps).grid(row=2, column=3)
        tk.Label(f, text=u'Y方向网格数').grid(row=3, column=0)
        tk.Entry(f, width=5, textvariable=self.numY).grid(row=3, column=1)
        tk.Label(f, text=u'周期(秒)').grid(row=3, column=2)
        tk.Entry(f, width=6, textvariable=self.totalT).grid(row=3, column=3)
        f.grid()

    def getdepthname(self):
        filename = askopenfilename(title=u'选择数据文件')
        self.depthname.set(filename.split('/')[-1])
        tmp = np.loadtxt(self.depthname.get())
        self.numX.set(tmp.shape[1])
        self.numY.set(tmp.shape[0])

    def getboundryname(self):
        filename = askopenfilename(title=u'选择数据文件')
        self.boundryname.set(filename.split('/')[-1])

    def makeSelect(self, frame):
        f = tk.Frame(frame)

        tk.Checkbutton(f, text=u'平流项', variable=self.enableAdvection).grid(row=4, column=0, sticky=tk.W)
        tk.Checkbutton(f, text=u'粘性项', variable=self.enableViscosity).grid(row=5, column=0, sticky=tk.W)
        tk.Checkbutton(f, text=u'底摩擦', variable=self.enableFriction).grid(row=6, column=0, sticky=tk.W)
        tk.Radiobutton(f, text=u'忽略科氏力', value=0, variable=self.enableCoriolisForce).grid(row=7, column=0, sticky=tk.W)
        tk.Radiobutton(f, text=u'常数科氏力', value=1, variable=self.enableCoriolisForce).grid(row=8, column=0, sticky=tk.W)
        tk.Radiobutton(f, text=u'可变科氏力', value=2, variable=self.enableCoriolisForce).grid(row=9, column=0, sticky=tk.W)

        tk.Label(f, text=u'分辨率(分)').grid(row=4, column=1, sticky=tk.E)
        tk.Label(f, text=u'涡动粘性系数').grid(row=5, column=1, sticky=tk.E)
        tk.Label(f, text=u'底摩擦系数').grid(row=6, column=1, sticky=tk.E)
        tk.Label(f, text=u'半隐半显系数').grid(row=7, column=1, sticky=tk.E)
        tk.Label(f, text=u'经度').grid(row=8, column=1, sticky=tk.E)
        tk.Label(f, text=u'纬度').grid(row=9, column=1, sticky=tk.E)

        tk.Entry(f, width=7, textvariable=self.resolution).grid(row=4, column=2, sticky=tk.W)
        tk.Entry(f, width=7, textvariable=self.viscosityCoefficient).grid(row=5, column=2, sticky=tk.W)
        tk.Entry(f, width=7, textvariable=self.bottomFriction).grid(row=6, column=2, sticky=tk.W)
        tk.Entry(f, width=7, textvariable=self.alpha).grid(row=7, column=2, sticky=tk.W)
        tk.Entry(f, width=7, textvariable=self.longitude).grid(row=8, column=2, sticky=tk.W)
        tk.Entry(f, width=7, textvariable=self.latitude).grid(row=9, column=2, sticky=tk.W)

        f.grid()

    def makeButton(self, frame):
        f = tk.Frame(frame)
        tk.Button(f, text=u'绘制水深', command=self.plotDepth, font=('YaHei', 15)).grid(row=10, column=0)
        tk.Button(f, text=u'开始模拟', command=self.work, font=('YaHei', 15)).grid(row=10, column=1)
        tk.Button(f, text=u'绘制动画', command=self.plotAnim, font=('YaHei', 15)).grid(row=11, column=0)
        tk.Button(f, text=u'保存动画', command=self.saveAnim, font=('YaHei', 15)).grid(row=11, column=1)
        tk.Button(f, text=u'绘同潮图', command=self.tongchaotu, font=('YaHei', 15)).grid(row=12, column=0)
        tk.Button(f, text=u'清空数据', command=self.clean, font=('YaHei', 15)).grid(row=12, column=1)
        f.grid()

    def clean(self):
        self.index = 0
        files = os.listdir('.')
        for i in files:
            if i.split('.')[-1] in ['txt', 'png', 'gif', 'mod']:
                os.remove(i)

    def work(self):
        try:
            self.realWork()
        except:
            tkMessageBox.showwarning(u"启动数值模式", u"无法运行数值模式，请检查参数！")

    def plotDepth(self):
        try:
            self.realPlotDepth()
        except:
            tkMessageBox.showwarning(u"绘制水深", u"无水深数据或数据有误！")

    def plotAnim(self):
        try:
            self.realPlotAnim()
        except:
            tkMessageBox.showwarning(u"播放动画", u"未运行模式或数据有误！")

    def tongchaotu(self):
        try:
            self.realTongchaotu()
        except:
            tkMessageBox.showwarning(u"同潮图", u"未运行模式或数据有误！")

    def realPlotDepth(self):
        if self.enableCoriolisForce.get() != 0:
            flag = True
        else:
            flag = False
        name_h = self.depthname.get()
        lon = self.longitude.get()
        lat = self.latitude.get()
        res = self.resolution.get() / 60.
        x = self.numX.get()
        y = self.numY.get()
        h = np.loadtxt(name_h)

        LON = np.linspace(lon, lon + res*x, x)
        LAT = np.linspace(lat + res*y, lat, y)
        lons, lats = np.meshgrid(LON, LAT)
        levels = np.linspace(eps, h.max(), 15)

        top = tk.Toplevel()
        top.title(u'水深分布图')

        my_cmap = matplotlib.cm.get_cmap('rainbow')
        fig = plt.figure(figsize=(self.width.get(), self.width.get()*ratio))
        ax = fig.add_axes([.02, .02, 1, 0.9])
        canvas = FigureCanvasTkAgg(fig, master=top)
        canvas.get_tk_widget().grid()

        plt.title('Water Depth', fontsize=16)
        if flag:
            m = Basemap(projection='merc', llcrnrlat=lat, urcrnrlat=lat+y*res, \
                        llcrnrlon=lon, urcrnrlon=lon+x*res, resolution='l', ax=ax)
            m.fillcontinents(color='black', lake_color='black')
            m.drawcoastlines(linewidth=2.0)
            m.drawmapboundary()
            m.drawcountries(linewidth=1.0)
            m.contour(lons, lats, h, levels, linewidths=0.5, colors='k', latlon=True)
            CS = m.contourf(lons, lats, h, levels, cmap=my_cmap, latlon=True)
        else:
            ax.patch.set_color('0.25')
            ax.contour(lons, lats, h, levels, linewidths=0.5, colors='k')
            CS = ax.contourf(lons, lats, h, levels, cmap=my_cmap)
        ticks = [float('%.0f' % i) for i in levels[::2]]
        cb = plt.colorbar(CS, drawedges=True, ticks=ticks)
        cb.ax.set_yticklabels(['%.0f' % i for i in ticks])

    def realWork(self):
        self.index += 1
        cmd = 'tide.exe'
        cmd += ' ' + str(self.depthname.get())
        cmd += ' ' + str(self.boundryname.get())
        cmd += ' ' + str(self.numX.get())
        cmd += ' ' + str(self.numY.get())
        cmd += ' ' + str(self.numSteps.get())
        cmd += ' ' + str(self.resolution.get())
        cmd += ' ' + str(self.totalT.get())
        cmd += ' ' + str(self.enableAdvection.get())
        cmd += ' ' + str(self.enableViscosity.get())
        cmd += ' ' + str(self.enableFriction.get())
        cmd += ' ' + str(self.enableCoriolisForce.get())
        cmd += ' ' + str(self.autoXStep.get())
        cmd += ' ' + str(self.viscosityCoefficient.get())
        cmd += ' ' + str(self.bottomFriction.get())
        cmd += ' ' + str(self.alpha.get())
        cmd += ' ' + str(self.latitude.get())
        cmd += ' ' + str(self.index)
        print cmd
        returnValue = call(cmd)
        if returnValue == 0:
            tkMessageBox.showinfo(u"恭喜你", u"模式运行成功！")
        else:
            tkMessageBox.showwarning(u"错误", u"模式运行失败，请检查参数！")

    def realPlotAnim(self):
        if self.enableCoriolisForce.get() != 0:
            flag = True
        else:
            flag = False
        my_cmap = matplotlib.cm.get_cmap('rainbow')
        name_zeta = 'output_Zeta_' + ('%03d' % self.index) + '.txt'
        name_u = 'output_U_' + ('%03d' % self.index) + '.txt'
        name_v = 'output_V_' + ('%03d' % self.index) + '.txt'
        x = self.numX.get()
        y = self.numY.get()
        s = self.numSteps.get()
        lon = self.longitude.get()
        lat = self.latitude.get()
        res = self.resolution.get() / 60.

        zeta = np.loadtxt(name_zeta)
        zeta = zeta.reshape((s, y, x))
        u = np.loadtxt(name_u)
        u = u.reshape((s, y, x))
        v = np.loadtxt(name_v)
        v = v.reshape((s, y, x))
        u = np.ma.array(u, mask=u < invalid)
        v = np.ma.array(v, mask=v < invalid)

        LON = np.linspace(lon, lon + res*x, x)
        LAT = np.linspace(lat + res*y, lat, y)

        lons, lats = np.meshgrid(LON, LAT)
        min_v = zeta[zeta > invalid].min()
        max_v = zeta[zeta > invalid].max()
        max_speed = np.sqrt((u*u + v*v).max())
        levels = np.linspace(min_v, max_v, 31)

        top = tk.Toplevel()
        top.title(u'水位&流场可视化')
        fig = plt.figure(figsize=(self.width.get(), self.width.get()*ratio))
        ax = fig.add_axes([.02, .02, 1, 0.85])
        canvas = FigureCanvasTkAgg(fig, master=top)
        canvas.get_tk_widget().grid()

        def init():
            global CS1, CS2, Q, m, txt
            if flag:
                m = Basemap(projection='merc', llcrnrlat=lat, urcrnrlat=lat+y*res, \
                            llcrnrlon=lon, urcrnrlon=lon+x*res, resolution='l', ax=ax)
                m.fillcontinents(color='black', lake_color='black')
                m.drawcoastlines(linewidth=1.0)
                m.drawmapboundary()
                m.drawcountries(linewidth=1.0)
                CS1 = m.contour(lons, lats, zeta[0, :, :], levels, linewidths=0.5, colors='k', latlon=True)
                CS2 = m.contourf(lons, lats, zeta[0, :, :], levels, cmap=my_cmap, latlon=True)
                Q = m.quiver(lons, lats, u[0, :, :], v[0, :, :], zorder=10, latlon=True)
            else:
                ax.patch.set_color('.25')
                CS1 = ax.contour(lons, lats, zeta[0, :, :], levels, linewidths=0.5, colors='k')
                CS2 = ax.contourf(lons, lats, zeta[0, :, :], levels, cmap=my_cmap)
                Q = ax.quiver(lons, lats, u[0, :, :], v[0, :, :])

            cb = plt.colorbar(CS2, drawedges=True, ticks=levels[::2])
            cb.set_clim(min_v, max_v)
            cb.ax.set_yticklabels(['%.1f' % i for i in levels[::2]])
            txt = ax.set_title('Water Level and Flow Velocity at Step ' + str(1), fontsize=15)
            plt.quiverkey(Q, 1.15, 1.05, max_speed, ('%.2f m/s' % max_speed), labelpos='W')

        def updatefig(nt):
            global CS1, CS2, Q, txt, m
            for c in CS1.collections:
                c.remove()
            for c in CS2.collections:
                c.remove()
            if flag:
                CS1 = m.contour(lons, lats, zeta[nt, :, :], levels, linewidths=0.5, colors='k', latlon=True)
                CS2 = m.contourf(lons, lats, zeta[nt, :, :], levels, cmap=my_cmap, latlon=True)
                Q.set_UVC(u[nt, :, :], v[nt, :, :])
            else:
                CS1 = ax.contour(lons, lats, zeta[nt, :, :], levels, linewidths=0.5, colors='k')
                CS2 = ax.contourf(lons, lats, zeta[nt, :, :], levels, cmap=my_cmap)
                Q.remove()
                Q = ax.quiver(lons, lats, u[nt, :, :], v[nt, :, :])
            txt.set_text('Water Level and Flow Velocity at Step ' + str(nt + 1))

        self.anim = animation.FuncAnimation(fig, updatefig, frames=s, init_func=init)

    def saveAnim(self):
        try:
            self.anim.save('animation' + ('%03d' % self.index) + '.gif', writer='imagemagick', fps=10)
            # self.anim.save('animation.mp4', fps=20, extra_args=['-vcodec', 'libx264'], writer=animation.FFMpegWriter())
        except:
            tkMessageBox.showwarning(u"保存动画", u"意外错误，无法保存动画！")

    def realTongchaotu(self):
        if self.enableCoriolisForce.get() != 0:
            flag = True
        else:
            flag = False
        amplitude = np.loadtxt('amplitude.txt')
        amplitude = np.ma.array(amplitude, mask=amplitude < invalid)
        arg = np.loadtxt('arg.txt')
        arg = np.ma.array(arg, mask=arg < invalid)
        lon = self.longitude.get()
        lat = self.latitude.get()
        res = self.resolution.get() / 60.
        x = self.numX.get()
        y = self.numY.get()

        LON = np.linspace(lon, lon + res*x, x)
        LAT = np.linspace(lat + res*y, lat, y)
        lons, lats = np.meshgrid(LON, LAT)
        levels1 = np.linspace(eps, amplitude.max(), 15)
        levels2 = np.linspace(eps, arg.max(), 15)

        top = tk.Toplevel()
        top.title(u'同潮图')

        fig = plt.figure(figsize=(self.width.get(), self.width.get()*ratio))
        ax = fig.add_axes([.05, .05, 0.9, 0.9])
        canvas = FigureCanvasTkAgg(fig, master=top)
        canvas.get_tk_widget().grid()

        def save():
            fig.savefig(u'同潮图' + str('%03d' % self.index) + u'.png', dpi=300)
        tk.Button(master=top, text=u'保存图片', font=('YaHei', 15), command=save).grid()

        if flag:
            m = Basemap(projection='merc', llcrnrlat=lat, urcrnrlat=lat+y*res, \
                        llcrnrlon=lon, urcrnrlon=lon+x*res, resolution='l', ax=ax)
            m.fillcontinents(color='aqua', lake_color='aqua')
            m.drawcoastlines(linewidth=2.0)
            m.drawmapboundary()
            m.drawcountries(linewidth=1.0)
            m.contour(lons, lats, amplitude, levels1, linewidths=1, colors='b', latlon=True, linestyles='dashed')
            m.contour(lons, lats, arg, levels2, linewidths=1, latlon=True)
        else:
            ax.patch.set_color('0.5')
            ax.contour(lons, lats, amplitude, levels1, linewidths=1, colors='green', linestyles='dashed')
            ax.contour(lons, lats, arg, levels2, linewidths=1)

if __name__ == '__main__':
    GUI().mainloop()
