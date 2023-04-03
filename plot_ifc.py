#!/usr/bin/env python3

# IFC vs distanceをプロットする

class constant:
    ang_to_bohr=1.8897259886
    Ry_to_eV=13.605698066
    bohr_to_ang=0.529177249

# ==============
# こちらはfcsファイルからfcsの値と距離情報を読み込む.
# 

#クラス
class FcsDistance:
  def __init__(self,globalindex:int, fcs:float, distance:float):
    self.globalindex=globalindex
    self.fcs =fcs
    self.distance=distance


class Fcs_distance():
    def __init__(self,fcs_filename):
        self.__fcs_filename=fcs_filename # fcsファイル
        # self.maxorder=maxorder # maxorderを外部から入れる場合
    
    
    # fcsファイルのmaxorderを取得する．
    def get_order(self):
        import numpy as np
        self.check_order=np.zeros(5) # 2-6次があれば1，なければ0
        for i in range(6): #max6次までやってIFCsを調べる.
            # print("次数:: ", i+2)
            f = open(self.__fcs_filename, 'r')
            while True:
                data = f.readline()
                if data == "":
                    break
                if "*FC"+str(i+2) in data and "**FC"+str(i+2) not in data: 
                    self.check_order[i]=1
                    # print(data)
        # check_orderの最大次数を取得
        self.maxorder=int(np.sum(self.check_order))
        print(" ----------------------------------")
        print(" maxorder is ::", self.maxorder+1)
        print("")
    #
    # fcsの数(local index)を取得する．
    def get_localindex(self):
        # fcs用のcounter
        self.localindex=[-1 for i in range(self.maxorder)]
        # 1回目の読み込みでlocal indexを読み込む
        # 時間はかかるが,次数ごとに一回読み込むようにした方が確実．
        for i in range(self.maxorder):
            # print("次数:: ", i)
            f = open(self.__fcs_filename, 'r')
            while True:
                data = f.readline()
                if self.localindex[i]>=0:
                    if data =="\n":
                        break
                    self.localindex[i]+= 1
                    #print(data)
                    tmp=data.split()
                    # print(i, int(tmp[1]))
                # *FC?を見つけたらカウントを開始
                if "*FC"+str(i+2) in data and "**FC"+str(i+2) not in data: 
                    self.localindex[i]=0
        print(" ----------------------------------")
        print(" local index for each order ::" , self.localindex)
        print("")

    # fcsの値と原子間距離を取得する．
    def get_fcs(self):
        # fcs用のcounter
        fcs_count=[-1 for i in range(self.maxorder)]
        # 出力するfcs
        self.force_constant_with_distance=[[] for y in range(self.maxorder)]
        # 読み込み
        for i in range(self.maxorder):
            # print("次数:: ", i)
            f = open(self.__fcs_filename, 'r')
            while True:
                data = f.readline()
                #            
                if fcs_count[i]>=0:
                    if data == "\n":
                        break
                    tmp=data.split()
                    self.force_constant_with_distance[i].append(FcsDistance(globalindex=int(tmp[0]),fcs=float(tmp[2]),distance=float(tmp[6+i])))
                #
                #  *FC?を見つけたらカウントを開始
                if "*FC"+str(i+2) in data and "**FC"+str(i+2) not in data: 
                    fcs_count[i]= 0
        print(" ----------------------------------")
        print(" finish loading fcs and distances ")
        print("")
    #
    # 得られたfcsをプロット
    def make_plots(self):
        import matplotlib.pyplot as plt
        import numpy as np
        # プロットを作成するにあたり，横軸は最も長い距離までに制限しておいた方が便利に思う．
        # とりあえずは2次IFCsの最大値を取得してこれを利用する．
        self.__maxlength=max(np.array([self.force_constant_with_distance[0][j].distance*constant.bohr_to_ang for j in range(self.localindex[0]) ]))
        print(self.__maxlength)
        # 
        for i in range(self.maxorder):
            fig, ax = plt.subplots(figsize=(8,5),tight_layout=True) # figure, axesオブジェクトを作成
            x = np.array([self.force_constant_with_distance[i][j].distance*constant.bohr_to_ang for j in range(self.localindex[i]) ])
            y = np.array([self.force_constant_with_distance[i][j].fcs*constant.Ry_to_eV/pow(constant.bohr_to_ang,i+2) for j in range(self.localindex[i]) ])
            ax.set_xlim(0,10)
            ax.set_xlabel("distance A",fontsize=22)
            ax.set_ylabel("IFC[eV/A^"+str(i+2)+"]",fontsize=22)
            ax.scatter(x,np.abs(y),label=str(i+2)+"th IFC")
            ax.legend(loc="upper right",fontsize=15 )
            fig.show()
            fig.savefig(str(i+2)+"th_order_ifc.pdf")

    # 全ての関数をまとめる．
    def process(self):
        self.get_order()
        self.get_localindex()
        self.get_fcs()
        self.make_plots()
        return 0



# DEPRECATED :: now we use FCs_distance instead of load_fcs
def load_fcs(fcs_filename):
    # 
    maxorder=5
    
    # fcs用のcounter
    counter_fcs=[-1 for i in range(maxorder)]

    # 1回目の読み込みでlocal indexを読み込む
    # 時間はかかるが,次数ごとに一回読み込むようにした方が確実．
    for i in range(maxorder):
        print("次数:: ", i)
        f = open(fcs_filename, 'r')
        while True:
            data = f.readline()
            if counter_fcs[i]>=0:
                if data =="\n":
                    break
                counter_fcs[i]+= 1
                #print(data)
                tmp=data.split()
                # print(i, int(tmp[1]))
            # *FC?を見つけたらカウントを開始
            if "*FC"+str(i+2) in data and "**FC"+str(i+2) not in data: 
                counter_fcs[i]=0
                print(data)
    print("local index", counter_fcs)


    # fcs用のcounter
    counter_fcs=[-1 for i in range(maxorder)]
    # 出力するfcs
    force_constant_with_distance=[[] for y in range(maxorder)]

    # 2回目の読み込みでfcsを読み込む?
    for i in range(maxorder):
        print("次数:: ", i)
        f = open(fcs_filename, 'r')
        while True:
            data = f.readline()
            #            
            if counter_fcs[i]>=0:
                if data == "\n":
                    break
                tmp=data.split()
                force_constant_with_distance[i].append(FcsDistance(globalindex=int(tmp[0]),fcs=float(tmp[2]),distance=float(tmp[6+i])))
            #
            #  *FC?を見つけたらカウントを開始
            if "*FC"+str(i+2) in data and "**FC"+str(i+2) not in data: 
                counter_fcs[i]= 0
    return force_constant_with_distance




# fcs_filename="../plot_EPS/lasso_300_merged/TiO2224_anharm.fcs"



def parse_cml_args(cml):
    '''
    Command line parser
    '''
    import argparse

    arg = argparse.ArgumentParser(add_help=True)

    arg.add_argument('file',
                     help='ALAMODE ifc file (.fcs)')
    # arg.add_argument('-p', '--input', dest='poscar', action='store', type=str,
    #                  default='POSCAR',
    #                  help='POSCAR of equilibrium lattice. Default is POSCAR ')

    return arg.parse_args(cml)


def main(fcs_filename):
    print("")
    print(" input filename :: {0}".format(fcs_filename))
    print("")
    plot=Fcs_distance(fcs_filename)
    plot.process()
    return 0


if __name__ == '__main__':
    '''
       Simple script for plotting IFCs and Distances relations.
      Usage:
      $ python plot_ifc.py file1.fcs

      For details of available options, please type
      $ python plot_ifc.py -h
    '''
    print("*****************************************************************")
    print("                       plot_ifc.py                               ")
    print("                      Version. 0.0.1                             ")
    print("*****************************************************************")
    print("")
    
    import sys
    arg = parse_cml_args(sys.argv[1:])     # parse commandline
    fcs_filename=arg.file

    if fcs_filename == "":
        print("ERROR :: alm file(.fcs or .xml) is not specified ")
        print("For details of usage, please type\n$ plot_ifc.py -h")
        exit(1)
    if not fcs_filename.endswith("fcs"):
        print("WARNING :: a filename does not end with fcs" )

    main(fcs_filename)
