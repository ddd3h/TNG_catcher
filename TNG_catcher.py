try:
    import illustris_python as il
    import numpy as np
    import pandas as pd
    import os
    import matplotlib.pyplot as plt
    from matplotlib.colors import SymLogNorm

            
except ModuleNotFoundError:
    print("The following modules is not installed at least:")
    print("\t- illustris_python")
    print("\t- numpy")
    print("\t- pandas")
    print("\t- os")
    print("\t- matplotlib")
    print('"llustris_python" is required to put in the same directory as this file.')
else:
    print("The following modules installed.")
    print("\t- import illustris_python as il")
    print("\t- import numpy as np")
    print("\t- import pandas as pd")
    print("\t- import os")
    print("\t- import matplotlib.pyplot as plt")
    print("\t- from matplotlib.colors import LogNorm")

class extract_subhalo_data:
    def __init__(self, basePath):
        self.basePath = basePath
        self.FoFnum = 0.5  # default
        self.plottype = 0  # default (plot type: gas)
        self.output_dirname = "cut_data"
        self.output_file_format = "npz"
        self.__is_subhalo = True
        self.subhalolist = [342447]
        self.fields = ("Coordinates", "Masses")
        self.__fields = ("x", "y", "z", "Masses")
        self.snapnum = 99
        
    def change_subhalo_mode(self, is_subhalo = True):
        self.__is_subhalo = is_subhalo
        if self.__is_subhalo:
            print("Mode: Subhalo")
            print("The only ID to set with `get_list` and `set_list` is the Subhalo ID; if you want to make it a Halo ID, change it to **False**.")
        else:
            print("Mode: Halo")
            print("The only ID to set with `get_list` and `set_list` is the Halo ID; if you want to make it a Subhalo ID, change it to **True**.")

    def get_list(self, filename):
        with open(filename, mode="r") as f:
            self.subhalolist = [int(i) for i in f.readlines()]

    def set_list(self, subhalolist):
        if type(subhalolist) != list:
            raise TypeError("Argument of set_list must be list.")
        else:
            self.subhalolist = subhalolist

    def set_fields(self, *fields):
        if isinstance(fields[0], list):
            self.fields = tuple(fields[0])
        else:
            self.fields = fields
        self.__in_coordinates = False
        _fields = []
        for i in self.fields:
            if i == "Coordinates":
                _fields.append("x")
                _fields.append("y")
                _fields.append("z")
                self.__in_coordinates = True
            elif i == "GFM_Metals":
                for j in ["H","He","C","N","O","Ne","Mg","Si","Fe","MetalTotal"]:
                    _fields.append(f"{j}")
            else:
                _fields.append(i)
        self.__fields = tuple(_fields)
        return self.__fields

    def set_snapshot(self, snapnum):
        self.snapnum = snapnum

    def set_FoF(self, FoFnum):
        self.FoFnum = FoFnum

    def set_plottype(self, plottype):
        self.plottype = plottype

    def set_output_dirname(self, output_dirname):
        self.output_dirname = output_dirname

    def set_output_file_format(self, output_file_format):
        if output_file_format not in ["npz", "pickel", "csv"]:
            raise ValueError("output_file_format must be 'npz' or 'pickel', 'csv'.")
        self.output_file_format = output_file_format

    def run(self):
        os.makedirs(self.output_dirname, exist_ok=True)
        for i in self.subhalolist:
            if self.__is_subhalo:
                pos = il.groupcat.loadSingle(self.basePath, self.snapnum, subhaloID=i)[
                    "SubhaloPos"
                ]
                dis = (
                    il.groupcat.loadSingle(self.basePath, self.snapnum, subhaloID=i)[
                        "SubhaloHalfmassRad"
                    ]
                    * self.FoFnum
                )
                sunhaloData = il.snapshot.loadSubhalo(
                    self.basePath, self.snapnum, i, self.plottype, fields=self.fields
                )
            else:
                # 少し変更が必要
                pos = il.groupcat.loadSingle(self.basePath, self.snapnum, haloID=i)[
                    "GroupPos"
                ]
                dis = (
                    il.groupcat.loadSingle(self.basePath, self.snapnum, haloID=i)[
                        "Group_R_Crit500"
                    ]
                    * self.FoFnum * 2
                )
                sunhaloData = il.snapshot.loadHalo(
                    self.basePath, self.snapnum, i, self.plottype, fields=self.fields
                )

            df_dict = {}
            for j in self.__fields:
                if j == "x":
                    df_dict[j] = np.array(sunhaloData["Coordinates"].T[0])
                elif j == "y":
                    df_dict[j] = np.array(sunhaloData["Coordinates"].T[1])
                elif j == "z":
                    df_dict[j] = np.array(sunhaloData["Coordinates"].T[2])
                elif j == "H":
                    df_dict[j] = np.array(sunhaloData["GFM_Metals"].T[0])
                elif j == "He":
                    df_dict[j] = np.array(sunhaloData["GFM_Metals"].T[1])
                elif j == "C":
                    df_dict[j] = np.array(sunhaloData["GFM_Metals"].T[2])
                elif j == "N":
                    df_dict[j] = np.array(sunhaloData["GFM_Metals"].T[3])
                elif j == "O":
                    df_dict[j] = np.array(sunhaloData["GFM_Metals"].T[4])
                elif j == "Ne":
                    df_dict[j] = np.array(sunhaloData["GFM_Metals"].T[5])
                elif j == "Mg":
                    df_dict[j] = np.array(sunhaloData["GFM_Metals"].T[6])
                elif j == "Si":
                    df_dict[j] = np.array(sunhaloData["GFM_Metals"].T[7])
                elif j == "Fe":
                    df_dict[j] = np.array(sunhaloData["GFM_Metals"].T[8])
                elif j == "MetalTotal":
                    df_dict[j] = np.array(sunhaloData["GFM_Metals"].T[9])
                else:
                    df_dict[j] = np.array(sunhaloData[j])

            df = pd.DataFrame(df_dict)
            df = df[
                (pos[0] - dis < df["x"])
                & (df["x"] < pos[0] + dis)
                & (pos[1] - dis < df["y"])
                & (df["y"] < pos[1] + dis)
                & (pos[2] - dis < df["z"])
                & (df["z"] < pos[2] + dis)
            ].reset_index(drop=True)

            if self.output_file_format == "npz":
                df_dict = {}
                for j in self.__fields:
                    df_dict[j] = np.array(df[j])
                if self.__is_subhalo:
                    np.savez(f"{self.output_dirname}/subhalo{i:06}_data.npz", **df_dict)
                else:
                    np.savez(f"{self.output_dirname}/halo{i:03}_data.npz", **df_dict)
            print("Complete output: ", i)
        print("#### Complete all output ####")

    def info(self):
        print("fields: ", self.fields)
        print("## You can set fields of here: https://www.tng-project.org/data/docs/specifications/#parttype0")
        print("")
        print("basePath: \t\t", self.basePath)
        print("snapnum: \t\t", self.snapnum)
        print("FoFnum: \t\t", self.FoFnum)
        print("plottype: \t\t", self.plottype)
        print("output_dirname: \t", self.output_dirname)
        print("output_file_format: \t", self.output_file_format)
        print("is_subhalo: \t\t", self.__is_subhalo)
        print("subhalolist: \t\t", self.subhalolist)


class plot_tools:
    def __init__(self, data):
        self.__data = dict(data)
        self.__files_list = data.files
        self.__rotation_ok = False

        if all(elem in self.__files_list for elem in ('x', 'y', 'z')):
            self.__rotation_ok = True
            self.__rotation_record = []

        # self.__file_format = data.split(".")[-1]
    
    def simple_plot(self,x,y,z,bins=301,log=True,fast=False, savefig=""):
        x,y,z = self.__data[x], self.__data[y], self.__data[z]
        xbins = np.linspace(x.min(), x.max(), bins)
        ybins = np.linspace(y.min(), y.max(), bins)
        hist, xed, yed = np.histogram2d(x, y, bins=(xbins, ybins), weights=z)
        if fast:
            if log:
                plt.imshow(hist.T, origin='lower', norm=SymLogNorm(linthresh=hist[hist>0].min(),vmin=hist[hist>0].min(), vmax=hist.max()))
                plt.axis("off")
                plt.colorbar()
                plt.show()
            else:
                plt.imshow(hist.T, origin='lower')
                plt.axis("off")
                plt.colorbar()
                plt.show()
        else:
            if log:
                plt.pcolormesh(xed, yed, hist.T, norm=SymLogNorm(linthresh=hist[hist>0].min(),vmin=hist[hist>0].min(), vmax=hist.max()))
                plt.colorbar()
                plt.show()
            else:
                plt.pcolormesh(xed, yed, hist.T)
                plt.colorbar()
                plt.show()

        if savefig:
            counter = 0
            if log:
                plt.pcolormesh(xed, yed, hist.T, norm=SymLogNorm(linthresh=hist[hist>0].min(),vmin=hist[hist>0].min(), vmax=hist.max()))
                plt.colorbar()
                while os.path.exists(savefig):
                    counter += 1
                    savefig = f"{counter}{savefig}"
                plt.savefig(savefig)
            else:
                plt.pcolormesh(xed, yed, hist.T)
                plt.colorbar()
                plt.savefig(savefig)
                while os.path.exists(savefig):
                    counter += 1
                    savefig = f"{counter}{savefig}"
                plt.savefig(savefig)


    def rot_x(self, theta):
        if self.__rotation_ok:
            theta_x = np.radians(theta)
            positions = np.array([self.__data['x'], self.__data['y'], self.__data['z']])
            rot_x = np.array([[ 1,                 0,                  0],
                          [ 0, np.cos(theta_x), -np.sin(theta_x)],
                          [ 0, np.sin(theta_x),  np.cos(theta_x)]])
            self.__data['x'], self.__data['y'], self.__data['z'] = rot_x @ positions
            self.__rotation_record.append(('x',theta))

    def rot_y(self, theta):
        if self.__rotation_ok:
            theta_y = np.radians(theta)
            positions = np.array([self.__data['x'], self.__data['y'], self.__data['z']])
            rot_y = np.array([[ np.cos(theta_y), 0,  np.sin(theta_y)],
                          [                 0, 1,                  0],
                          [-np.sin(theta_y), 0, np.cos(theta_y)]])
            self.__data['x'], self.__data['y'], self.__data['z'] = rot_y @ positions
            self.__rotation_record.append(('y',theta))

    def rot_z(self, theta):
        if self.__rotation_ok:
            theta_z = np.radians(theta)
            positions = np.array([self.__data['x'], self.__data['y'], self.__data['z']])
            rot_z = np.array([[ np.cos(theta_z), -np.sin(theta_z), 0],
                          [ np.sin(theta_z),  np.cos(theta_z), 0],
                          [                 0,                  0, 1]])
            self.__data['x'], self.__data['y'], self.__data['z'] = rot_z @ positions
            self.__rotation_record.append(('z',theta))
    
    def __rot_x(self,theta):
        theta = np.radians(theta)
        return np.array([[ 1,                 0,                  0],
                          [ 0, np.cos(theta), -np.sin(theta)],
                          [ 0, np.sin(theta),  np.cos(theta)]])
    
    def __rot_y(self,theta):
        theta = np.radians(theta)
        return np.array([[ np.cos(theta), 0,  np.sin(theta)],
                          [                 0, 1,                  0],
                          [-np.sin(theta), 0, np.cos(theta)]])
    
    
    def __rot_z(self,theta):
        theta = np.radians(theta)
        return np.array([[ np.cos(theta), -np.sin(theta), 0],
                          [ np.sin(theta),  np.cos(theta), 0],
                          [                 0,                  0, 1]])


    def save_rot(self,filename):
        if self.__rotation_ok:
            # with open(filename,'w') as f:
            #     for i in self.__rotation_record:
            #         f.write(f"{str(i)}\n")
            #     f.close()
            a = np.eye(3)
            for i in self.__rotation_record:
                if i[0] == 'x':
                    a = a @ self.__rot_x(i[1])
                elif i[0] == 'y':
                    a = a @ self.__rot_y(i[1])
                elif i[0] == 'z':
                    a = a @ self.__rot_z(i[1])
            np.save(filename,a)


    def save_npz(self,filename):
        if not os.path.exists(filename):
            np.savez(filename, **(self.__data))
            
            
    def get_temp(self):
        if 'Temperature' in self.__data.keys():
            return self.__data['Temperature']
        else:
            if 'InternalEnergy' in self.__data.keys() and 'ElectronAbundance' in self.__data.keys():
                return 