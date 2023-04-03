#!/usr/bin/env python3 

# plotrta.py
#
# Simple script to visualize Thermal conductivity from *.kl files.
#
#

import numpy as np
import matplotlib as mpl
from matplotlib.gridspec import GridSpec
#try:
#    mpl.use("Qt5agg")
#except:
#    pass
mpl.use('Pdf')
import matplotlib.pyplot as plt



import argparse    # 1. argparseをインポート
parser = argparse.ArgumentParser(description='ase wrapper  ')    # 2. パーサを作る
# 
# 3. parser.add_argumentで受け取る引数を追加していく
parser.add_argument('arg1'  , help='alamode *.kl files',nargs='+')    # 必須の引数を追加
parser.add_argument('--tmin', help="minimum Temperature[K]", type=float, default="0.0")    # オプション引数（指定しなくても良い引数）を追加(type=で型指定)
parser.add_argument('--tmax', help="maximum Temperature[K]", type=float, default="-1")    # オプション引数（指定しなくても良い引数）を追加(type=で型指定)
parser.add_argument('--ymin', help="minimum LTC[W/mK]", type=float, default="0.0")   
parser.add_argument('--ymax', help="maximum LTC[W/mK]", type=float, default="-1")    

parser.add_argument('-l', '--log', help="set logscale for xy axis", default="no" )
parser.add_argument('-a', '--axis', help="axis which is ploted. ", default="mean" )
#
#
args = parser.parse_args()    # 4. 引数を解析
#
# print('arg4='+args.arg4)




# font styles
mpl.rc('font', **{'family': 'Times New Roman', 'sans-serif': ['Helvetica']})
mpl.rc('xtick', labelsize=16)
mpl.rc('ytick', labelsize=16)
mpl.rc('axes', labelsize=16)
mpl.rc('lines', linewidth=1.5)
mpl.rc('legend', fontsize='small')
# line colors and styles
color = ['b', 'g', 'r', 'm', 'k', 'c', 'y', 'r']
lsty = ['-', '-', '-', '-', '--', '--', '--', '--']


# ファイルからデータを読み込み
def read_files(filenames):
    # 格納するデータ
    datas=[]
    # 
    for i,name in enumerate(filenames):
        datas.append(np.loadtxt(name))
    return datas
    

# 指定されなかった場合のtmaxを決定する．
def decide_tmax(datas):
    maxlist=[]
    for i in range(len(datas)):
        maxlist.append(np.max(datas[i][:,0]))
    return np.max(maxlist)

# yのmaxminを取り出す．
# この時，tmin,tmaxの範囲でのymin/ymaxを知る必要がある．
def get_y_minmax(datas,tmin,tmax):
    # 
    ymin = []
    ymax = []
    # 
    if plot_axis=="mean":
        for i in range(len(datas)):
            # まず，tmin,tmaxに応じたdatasのindexを取得
            min_indx=np.where(datas[i][:,0]>=tmin)[0][0]
            max_indx=np.where(datas[i][:,0]<=tmax)[0][-1]
            # print("tminより大きい要素", np.where(datas[i][:,0]>=tmin)[0][0]) #条件を満たす最初の要素
            # print("tmaxより小さい要素", np.where(datas[i][:,0]<=tmax)[0][-1]) #条件を満たす最後の要素
            #
            ymin.append(np.min((datas[i][min_indx:max_indx,1]+datas[i][min_indx:max_indx,5]+datas[i][min_indx:max_indx,9])/3))
            ymax.append(np.max((datas[i][min_indx:max_indx,1]+datas[i][min_indx:max_indx,5]+datas[i][min_indx:max_indx,9])/3))
        # 最終的に最小最大を出す    
        ymin=np.min(np.array(ymin))
        ymax=np.max(np.array(ymax))*1.1 # 1.1倍下駄を履かせておく
    return ymin, ymax



# plotする
def run_plot(datas,tmin,tmax,ymin,ymax):
    # tmaxの指定
    if tmax==-1:
        print(" ")
        print(" WARNING : tmax is not specified. ")
        print("       Use maximum values in *.kl ")
        print(" ")
        tmax=decide_tmax(datas)
    # ymaxの指定
    if ymax==-1:
        print(" ")
        print(" WARNING : ymax is not specified. ")
        print("       Automatic decision of ymax is used ")
        print(" ")
        ymin, ymax = get_y_minmax(datas,tmin,tmax)


    # datasにある複数のデータをプロットする，
    fig, ax = plt.subplots()
    ax.set_xlim(tmin,tmax) # plot range[K]
    ax.set_ylim(ymin,ymax) # plot range[K]


    # これはなんだ？
    # gs = GridSpec(nrows=1, ncols=nax )
    # gs.update(wspace=0.1)
    # 
    # 複数データのプロット
    for i in range(len(datas)):
        # ax = plt.subplot(gs[iax])
        #
        if plot_axis == "mean":
            ax.plot(datas[i][:,0], (datas[i][:,1]+datas[i][:,5]+datas[i][:,9])/3,linestyle=lsty[i], color=color[i], label=filenames[i])

        # 
        # 初めのデータの時にcaptionをつける
        if i == 0:
            ax.set_ylabel("LTC (W/mK)", labelpad=10)
            ax.set_xlabel("Temp (K)", labelpad=10)
        # else:
        #    ax.set_yticklabels([])
        #    ax.set_yticks([])

        # 
        # plt.axis([xmin_ax[iax], xmax_ax[iax], ymin, ymax])
        # ax.set_xticks(xticks_ax[iax])
        # ax.set_xticklabels(xticklabels_ax[iax])
        # ax.xaxis.grid(True, linestyle='-')

        #if options.print_key and iax == 0:
        #    ax.legend(loc='best', prop={'size': 10})
    ax.legend(loc='best', prop={'size': 10})
    # plt.tight_layout()
    # https://www.delftstack.com/ja/howto/matplotlib/save-figures-identical-to-displayed-figures-matplotlib/
    plt.savefig("ltc.pdf", bbox_inches='tight')
    plt.show()



if __name__ == '__main__':
    '''
    For details of available options, please type
    $ python plot_ltc.py -h
    '''
    #
    print("*****************************************************************")
    print("                        plot_ltc.py                              ")
    print("                      Version. 0.0.1                             ")
    print("*****************************************************************")
    print("")
    #
    args = parser.parse_args()    # 4. 引数を解析
    filenames = args.arg1
    tmin = args.tmin
    tmax = args.tmax
    ymin = args.ymin
    ymax = args.ymax
    plot_axis = args.axis


    if len(filenames) == 0:
        print("Usage: plotband.py [options] file1.bands file2.bands ...")
        print("For details of available options, please type\n$ python plotband.py -h")
        exit(1)
    else:
        print(" ## ====")
        print(" Number of files = %d" % len(filenames))
        print(" ")
    #
    if plot_axis == "mean":
        print(" ## ====")
        print(" axis=mean :: We plot (k_xx+k_yy+k_zz)/3 " )
        print(" ")


    # read *.kl
    datas=read_files(filenames)
    #
    #
    # decide ymin, ymax
    if ymax==-1:
        ymin, ymax = get_y_minmax(datas,tmin,tmax)

    # plot datas
    # とりあえずyminは0で固定
    run_plot(datas,tmin,tmax, 0, ymax)


    #nax, xticks_ax, xticklabels_ax, xmin_ax, xmax_ax, ymin, ymax, \
    #    data_merged_ax = preprocess_data(
    #        files, options.unitname, options.normalize_xaxis)

    #run_plot(nax, xticks_ax, xticklabels_ax,
    #         xmin_ax, xmax_ax, ymin, ymax, data_merged_ax)
    
