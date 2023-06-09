#!/usr/bin/env python3

'''
This script calculate the system energy for given atomic configurations.


usage 
plot_eps.py --alm=sample.xml --vasp=POSCAR file1 file2 ...

input
---------------
  --alm :: sample.xml :: IFCs file made by alm.
  --vasp:: POSCAR     :: VASP type atomic position file of PRIMITIVE cell without displacements.
  file1,file2,...     :: VASP type atomic position file of PRIMITIVE cell with displacements.
'''


def parse_cml_args(cml):
    '''
    Command line parser
    '''
    import argparse
    
    arg = argparse.ArgumentParser(add_help=True)

    arg.add_argument('files',nargs='+',
                     help='vasp displacement files')
    arg.add_argument('-p', '--input', dest='poscar', action='store', type=str,
                     default='POSCAR',
                     help='POSCAR of equilibrium lattice. Default is POSCAR ')
    arg.add_argument('-i', '--xml', dest='xml', action='store', type=str,
                     help='alm xml file to parse')
    arg.add_argument('-T', '--temperature', dest='temperature',
                     action='store', type=float,
                     default=300,
                     help='The temperature.')
    arg.add_argument('-f', '--maxorder', dest='maxorder',
                     action='store', type=int,
                     default=6,
                     help='maximum order in taylor expansion')

    return arg.parse_args(cml)


def main():
    import xml.etree.ElementTree as ET
    import sys

    # parse commandline
    arg = parse_cml_args(sys.argv[1:])

    #arg.unit == 'cm-1'
    #t0 = arg.frequency

    # read a xml file
    tree = ET.parse(arg.xml)
    
    # get top element
    root = tree.getroot()

    # load system information (maxorder::次数-1)
    fcs_phonon=Fcs_phonon(root=root,maxorder=arg.maxorder-1)

    # load fcs
    fcs_phonon.load_fcs_xml()

    # calculate energy of given files (arg.files)
    counter=0
    for filename in arg.files:
        print(" file number :: {0}  , file name :: {1}".format(counter,filename))
        # calculate displacement
        displacement=calc_displacement(arg.poscar,filename)
        # calculate max displacement (for graph)
        max_displacement=get_max_displace(arg.poscar,filename)
        energy[maxorder-1] # output energy
        for j in range(2,maxorder):
            energy[j-1]=fcs_phonon.calculate_energy(u0,j)        
        #
        print(filename, u_norm, end="")
        for j in range(2,maxorder):
            print(energy[j-1])
    #
    # for i in np.arange(1,21):
    #     i=str(i).zfill(2)
    #     calc_subtract(i, dir)
    #     filename=dir+"disp"+str(i)+".txt"
    #     u0=np.loadtxt(filename)
    #     u_norm=get_max_displace(filename)
    #     # print(u_norm)
    #     E2=fcs_phonon.calculate_energy(u0,2)
    #     E3=fcs_phonon.calculate_energy(u0,3)
    #     E4=fcs_phonon.calculate_energy(u0,4)
    #     E5=fcs_phonon.calculate_energy(u0,5)
    #     E6=fcs_phonon.calculate_energy(u0,6)
    #     print(u_norm,E2,E3,E4,E5,E6)
    return 0



if __name__ == '__main__':
    main()


# below are main codes
    
class constant:
    '''
    define nessesary physical constants
    '''
    ang_to_bohr=1.8897259886
    Ry_to_eV=13.605698066
    bohr_to_ang=0.529177249

    
def calc_displacement(poscar, disp):
    '''
    calculate atomic displacement from POSCAR and disp*.POSCAR.

    note
    ------------
    vasp uses angstrom, and alm uses bohr (atomic rydberg unit).
    '''
    import ase.io # !! caution: here I use ase
    
    # in angstrom (vasp)
    primitive=ase.io.read(poscar).get_positions()
    displace=ase.io.read(disp).get_positions()
    # ang to bohr (bohrで出力していることに注意．)
    subtract=(displace-primitive)*constant.ang_to_bohr
    # np.savetxt(dir+"disp"+str(i)+".txt", subtract)
    return subtract


def get_max_displace(poscar, disp):
    '''
     *.txtファイルを入力すると一番大きいdisplacement(これは全ての原子種の中で）を計算してくれる.
     これは実際にプロットを作る時に役に立つ.
    '''
    import numpy as np
    u_norm=np.linalg.norm(calc_displacement(poscar, disp),axis=1)
    return np.amax(u_norm)


# ---------------------------
# All classes below are from the alamode alm code.
#
class AtomCellSuper:
  def __init__(self,index:int, tran:int, cell_s:int):
    self.index=index
    self.tran =tran
    self.cell_s=cell_s

class FcsArrayWithCell:
  def __init__(self,fcs_val:float,pairs:AtomCellSuper):
    self.fcs_val=fcs_val
    self.pairs=pairs


class System:
    '''
    input
    ---------
     root  :: input xml file name
     
    '''

    def __init__(self,root):
        self.__root=root

    class Map: # for s2p
        def __init__(self,atom_num:int, tran_num:int):
            self.atom_num = atom_num
            self.tran_num = tran_num

    def load_system_info(self):
        '''
        load_system_info:
          - The number of atoms in supercell
          - The number of primitive cell
          - The number of atoms in primitive cell
        '''
        # 全原子数
        self.__nat:int=int(self.__root.find("Structure").find("NumberOfAtoms").text)
        # supercellの数
        self.__ntran=int(self.__root.find("Symmetry").find("NumberOfTranslations").text)
        #単位格子の原子数
        self.__natmin = self.__nat // self.__ntran
        #結果出力
        print(" nat,ntran,natmin=",self.__nat,self.__ntran,self.__natmin)

    def make_mapping(self):
        '''
        read map_p2s,map_s2p
        '''
        import numpy as np
        self.map_p2s=np.zeros([self.__natmin,self.__ntran])
        self.map_s2p=[[] for i in range(self.__nat)]
        for child in self.__root.find("Symmetry").find("Translations"):
            tran   =int(child.get("tran"))-1  # from 1-based index to 0-based index
            atom_p =int(child.get("atom"))-1  # from 1-based index to 0-based index
            atom_s =int(child.text)-1

            # primitive to super
            self.map_p2s[atom_p][tran] = atom_s
            # super to primitive
            self.map_s2p[atom_s]=self.Map(atom_num=atom_p,tran_num=tran)

    def __del__(self):
        del self.map_s2p
        del self.map_p2s
        del self.__nat
        del self.__ntran
        del self.__natmin
        del self.__root

class Fcs_phonon():
    def __init__(self,root,maxorder):
        self.__root=root # xml parse
        self.maxorder=maxorder
        self.system=System(root=self.__root) # we need System so that we can use map_s2p.
        self.system.load_system_info()
        self.system.make_mapping()
    
    def load_fcs_xml(self):
      from sympy.utilities.iterables import multiset_permutations
      # fcs
      self.force_constant_with_cell=[[] for y in range(self.maxorder)]
      # loop over IFC order
      for order in range(self.maxorder):
        print("order = ", order)
        if (order == 0):
          str_tag = "HARMONIC"
        else:
          str_tag = "ANHARM" + str(order + 2)
  
        # loop over IFC
        for child in self.__root.find("ForceConstants").find(str_tag):
          # get IFC
          fcs_val:float=float(child.text) 

          # initialization
          ivec_with_cell=[]
          ivec_tmp=[]
          ind=[] # for permutation. can be deprecated in future.

          # ======================
          # first get pair1
          tmp=[int(s) for s in child.get("pair1").split()]
          atmn:int =tmp[0]-1 # from 1-based index to 0-based index
          xyz:int  =tmp[1]-1 # from 1-based index to 0-based index
          ivec_with_cell.append(AtomCellSuper(index=3*self.system.map_p2s[atmn][0]+xyz,tran=0,cell_s=0)) #tran=0 is dummy,cel=0 in pair1
          ivec_pair1=AtomCellSuper(index=3*atmn+xyz,tran=0,cell_s=0) #pair1のみここで追加
          # after pair2
          for i in range(1,order+2):
            tmp=[int(s) for s in child.get("pair"+str(i+1)).split()]
            atmn:int   =tmp[0]-1
            xyz:int    =tmp[1]-1
            cell_s:int =tmp[2]-1
            ivec_with_cell.append(AtomCellSuper(index=3*atmn+xyz,tran=0,cell_s=cell_s)) #tran=0はdummy
            # ivec_tmp.append(AtomCellSuper(index=3*(map_s2p[atmn].atom_num - 1)+(xyz-1),tran=map_s2p[atmn].tran_num,cell_s=cell_s-1))
            ind.append(3*atmn+xyz)
          # ======================

          # ======================
          # permutation処理(次数だけのpairがある)
          # !! 現状cell_sの部分を同時に並び替えることができなかったので，代わりにpermutationの数を代入してある．
          # https://stackoverflow.com/questions/6284396/permutations-with-unique-values
          #for perm_list in multiset_permutations(ivec_with_cell[1:].index):
          #  for i in range(order+2-1):                #pair1を除いているので1引く
          #    atmn:int = int(perm_list[i].index / 3)  #3で割った商なので，atmnが出てくる.
          #    xyz:int  = perm_list[i].index % 3       #3でわった余なので，xyzが出てくる．
          #    ivec_tmp.append(AtomCellSuper(index=3*map_s2p[atmn].atom_num+xyz,tran=map_s2p[atmn].tran_num,cell_s=perm_list[i].cell_s))
          #
          #  force_constant_with_cell[order].append(FcsArrayWithCell(fcs_val=fcs_val,pairs=ivec_tmp))
          #
          multiplicity=len(list(multiset_permutations(ind)))

          for perm_list in multiset_permutations(ind):
            ivec_tmp=[ivec_pair1] #initialization
            # atmn_list= [ int(n /3) for n in perm_list]
            # xyz_list = [n % 3 for n in perm_list] 
            [ivec_tmp.append(AtomCellSuper(index=3*self.system.map_s2p[atmn].atom_num+xyz, tran=self.system.map_s2p[atmn].tran_num, cell_s=multiplicity)) for atmn,xyz in zip([ n//3 for n in perm_list], [n % 3 for n in perm_list])] 
            self.force_constant_with_cell[order].append(FcsArrayWithCell(fcs_val=fcs_val,pairs=ivec_tmp))
      #
      if not __debug__:
        for i in range(self.maxorder):
          print("  Number of non-zero IFCs for ", i + 2 , " order: (include permutation) ", len(self.force_constant_with_cell[i]))


    
    def calculate_energy_of_ith_order(self, u0, i:int):
      '''
      calculate energy of only (i+2)-th order contribution.
      '''
      import math
      length:int   = len(self.force_constant_with_cell[i])   #その次数に入っているIFCsの数  
      nelem:int    = i + 2;                             #次数,つまりIFCsにnelemだけの原子が関わっている.
      factorial    = math.factorial(nelem)              #エネルギーのtaylor展開につける係数.
      U:float      = 0
      #
      for j in range(length):   #loop over IFC at i-th order
        phi_val  = self.force_constant_with_cell[i][j].fcs_val
        dtmp     = 1.0/factorial * phi_val; #initialization
        #
        for k in range(nelem): #loop over IFC order
          (atmn,xyz) = divmod(self.force_constant_with_cell[i][j].pairs[k].index, 3)   #3で割った商はatmn(in primitive cell).余はxyz.
          #
          dtmp   *=  u0[atmn][xyz]  
        U += dtmp
      # energy in eV unit.
      return U*constant.Ry_to_eV #13.605698066

    def calculate_energy(self, u0, i:int):
      '''
      calculate_energy up to i-th order. 
      '''

      if (i<2):
        print("error:: i should be 2 or more. i corresponds i-th order energy.")
      if (self.maxorder<i-1):
        print("error::maxorder is too small::maxorder should be i or more")
        return 1
      if i==2: # 2nd order
        return self.calculate_energy_of_ith_order(u0, i-2)
      else:
        return self.calculate_energy_of_ith_order(u0, i-2)+self.calculate_energy(u0, i-1)

