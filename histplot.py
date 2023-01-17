
import matplotlib.pyplot as plt
import matplotlib as mpl

from matplotlib import colors, cm
import matplotlib.font_manager as fm


import numpy as np
import pandas as pd
import os

mpl.rcParams.update({'font.size': 18})

class custom_plot():
    def __init__(self):  # 类的初始化

        self.currentdir = os.getcwd()  # 得到当前文件夹路径
        self.base = "D:\\path\\to\\your\\file"
        self.base = "C:\\Users\\ZRZ\\Desktop\\python_scientific_plot"
        self.address = 'test.xlsx'
        self.fulladdress = os.path.join(self.base, self.address)
        self.waddresss = self.base


    def read_Information(self):
        self.df = pd.read_excel(self.fulladdress)

    def pandas_practice(self):

        self.blank_dic = {}

        """ uncomment above to see what is index"""
        for idx in range(len(self.df.index)):  # for loop


            if idx % 2 == 0:  # idx（行）是偶数
                self.blank_dic[idx] = self.df[idx]["column_name"]


        wdf = pd.DataFrame.from_dict(self.blank_dic, orient='index', columns=["column_name"])

        wdf.to_csv(self.waddresss + "output.csv", sep=',')


    def pandas_select(self):

        self.value1 = ""
        self.condition1 = ''
        self.column11 = ''
        self.column12 = ''

        self.x1 = self.df[self.df[self.condition1] == self.value1][self.column11].to_list()
        self.y1 = self.df[self.df[self.condition1] == self.value1][self.column12].to_list()



        self.value2 = ""
        self.condition2 = ''
        self.column21 = ''
        self.column22 = ''

        self.x2 = self.df[self.df[self.condition2] == self.value2][self.column21].to_list()
        self.y2 = self.df[self.df[self.condition2] == self.value2][self.column22].to_list()


    def plot_line_example(self):

        fig, axs = plt.subplots(1, 1, figsize=(3, 3))
        self.label1 = "Correct"
        self.label2 = "Session"
        self.title = "Sample of Plot Line"
        self.xlabel = r"$\mathregular{\lambda}$ (nm)"
        self.ylabel = 'y value'
        self.xlow = 0
        self.xhigh = 10
        self.ylow = 0
        self.yhigh = 10
        axs.plot(self.x1, self.y1, color="red", label=self.label1)
        axs.plot(self.x2, self.y2, color='green', label=self.label2)

        axs.legend()
        axs.set_title(self.title)
        axs.set_xlabel(self.xlabel, labelpad=10)
        axs.set_ylabel(self.ylabel, labelpad=10)


        font_names = [f.name for f in fm.fontManager.ttflist]
        print('font list：', font_names)

        mpl.rcParams['font.family'] = 'DejaVu Sans Mono'
        plt.rcParams['font.size'] = 18
        plt.rcParams['axes.linewidth'] = 2

        axs.spines['right'].set_visible(True)
        axs.spines['top'].set_visible(True)


        axs.xaxis.set_tick_params(which='major', size=10, width=2, direction='in', top='on')

        axs.xaxis.set_tick_params(which='minor', size=5, width=2, direction='in', top='on')

        axs.yaxis.set_tick_params(which='major', size=10, width=2, direction='in', right='on')

        axs.yaxis.set_tick_params(which='minor', size=5, width=2, direction='in', right='on')

        plt.show()

    def plot_line(self):

        fig, axs = plt.subplots(1, 1, figsize=(3, 3), sharex=True)
        self.label1 = ""

        self.label2 = " "
        self.title = ""
        self.xlabel = r"$\mathregular{\lambda}$ (nm)"

        self.ylabel = ''
        self.xlow = 0
        self.xhigh = 100
        self.ylow = 0
        self.yhigh = 100
        axs.plot(self.x1, self.y1, color="red", label=self.label1)
        axs.plot(self.x2, self.y2, color='green', label=self.label2)

        axs.legend()
        axs.set_title(self.title)
        axs.set_xlabel(self.xlabel, labelpad=10)
        axs.set_ylabel(self.ylabel, labelpad=10)
        axs.set_xlim(self.xlow, self.xhigh)
        axs.set_ylim(self.ylow, self.yhigh)
        #

        font_names = [f.name for f in fm.fontManager.ttflist]
        print('font list：', font_names)

        mpl.rcParams['font.family'] = 'Avenir'
        plt.rcParams['font.size'] = 18
        plt.rcParams['axes.linewidth'] = 2

        axs.spines['right'].set_visible(True)
        axs.spines['top'].set_visible(True)

        axs.xaxis.set_tick_params(which='major', size=15, width=2, direction='in', top='on')

        axs.xaxis.set_tick_params(which='minor', size=7, width=2, direction='in', top='on')
        axs.yaxis.set_tick_params(which='major', size=15, width=2, direction='in', right='on')
        axs.yaxis.set_tick_params(which='minor', size=7, width=2, direction='in', right='on')

        axs.xaxis.set_major_locator(mpl.ticker.MultipleLocator(50))
        axs.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(25))
        axs.yaxis.set_major_locator(mpl.ticker.MultipleLocator(50))
        axs.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(20))


        plt.show()

    def bar_1d(self):
        self.x1 = ["Vaccum Vessel", "Pressure Vessel", "HDPE Vessel", "Hydraulic Fluid", "Outer Jar", "Liquid Argon"]
        self.y1=[91.0, 7.95,0.866,0.010,0.011,0.121]
        distance = 0
        x_position = np.arange(len(self.x1))
        label1 = 'Correct'
        # label2 = ' Session'
        fig, axs = plt.subplots(1, 1, figsize=(3, 3))
        # axs.bar(x_position - distance, self.y1, 0.4, label=label1)
        axs.bar(x_position - distance, self.y1, 0.4)
        # axs.bar(x_position + distance, self.y2, 0.4, label=label2)
        # axs.legend()
        axs.set_xticks(x_position)
        axs.set_yscale("log")
        axs.set_xticklabels(self.x1)
        for i, v in enumerate(self.y1):
            axs.text(i-0.125, v, str(v), color='black')
        plt.xlabel("Physical Volume in SBC")
        plt.ylabel("Neutron Capture Ratio/%")
        plt.show()

    def color_hist_2d(self):

        bins = 40
        fig, axs = plt.subplots(1, 1, figsize=(3, 3))
        (h, xedge, yedge, image) = axs.hist2d(self.x1, self.y1, bins=bins, norm=colors.Normalize(), cmap='plasma')

        fig.colorbar(image, ax=axs)


        plt.show()

    def color_2d(self):

        self.array = [[5, 2], [7, 3]]
        self.xlow = 0
        self.xhigh = 10
        self.ylow = 0
        self.yhigh = 10
        fig, axs = plt.subplots(1, 1, figsize=(3, 3))
        image = axs.imshow(self.array, extent=(self.xlow, self.xhigh, self.ylow, self.yhigh), cmap='plasma')
        fig.colorbar(image, ax=axs)
        plt.show()

    def pandas_select_sample(self):

        self.x1 = self.df[self.df['Correct'] >= 6]['Name'].to_list()
        self.y1 = self.df[self.df['Correct'] >= 6]["Correct"].to_list()

        self.x2 = self.df[self.df['Correct'] >= 6]['Name'].to_list()
        self.y2 = self.df[self.df['Correct'] >= 6]["Session"].to_list()


if __name__ == "__main__":
    practice = custom_plot()
    # practice.read_Information()
    # practice.pandas_select_sample()
    # practice.plot_line_example()
    practice.bar_1d()
    # practice.x1 = practice.y2
    # practice.color_hist_2d()
    # practice.color_2d()


    # practice = custom_plot()
    # practice.read_Information()
    # practice.pandas_select()
    # practice.plot_line()
    # practice.bar_1d()
    # # practice.color_hist_2d()
    # # practice.color_2d()
    # 选择你需要的函数uncomment,如果有多个没有被comment, 那么会同时画出多个图片


