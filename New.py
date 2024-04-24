import illustris_python as il
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import h5py
import scienceplots
plt.style.use(["science","cjk-tc-font"])


subhalo_random_id = np.loadtxt("./random_subhalo.txt", dtype=int)
BASEPATH = "/home/nishihama/data/TNG50-1/output"
SNAPNUM = 99
Metal = {
    "H": 0,
    "He": 1,
    "C": 2,
    "N": 3,
    "O": 4,
    "Ne": 5,
    "Mg": 6,
    "Si": 7,
    "Fe": 8,
    "total": 9,
}

print("[DONE] IMPORTED parameters")

class SubhaloDataSet:
    def __init__(self, subhaloID):
        self.subhaloID = subhaloID
        self.DATA = {}
        self.basePath = "/home/nishihama/data/TNG50-1/output/"
        self.snapnum = 99
        with open("./eROSITAbubble_subhaloID.txt","r") as f:
            SubhaloIDList = f.readlines()
            self.SubhaloIDList = [int(i.replace("\n","")) for i in SubhaloIDList]
        self.abund = "aspl"
        self.Metals = [
            "H",
            "He",
            "C",
            "N",
            "O",
            "Ne",
            "Mg",
            "Si",
            "Fe",
            "total"
        ]

    
    def ChangeToCurrectCoordinates(self,x,y,z,cmx,cmy,cmz):
        HalfBoxSize = 35000/2
        if np.max(x) - np.min(x) > HalfBoxSize:
            x[x > HalfBoxSize] -= 35000
            if cmx > HalfBoxSize:
                cmx -= 35000
        if np.max(y) - np.min(y) > HalfBoxSize:
            y[y > HalfBoxSize] -= 35000
            if cmy > HalfBoxSize:
                cmy -= 35000
        if np.max(z) - np.min(z) > HalfBoxSize:
            z[z > HalfBoxSize] -= 35000
            if cmz > HalfBoxSize:
                cmz -= 35000
        return x,y,z,cmx,cmy,cmz
    
    
    def GetAllData(self):
        with h5py.File(
            self.basePath + f"snapdir_{self.snapnum:03}/snap_{self.snapnum:03}.0.hdf5",
            "r",
        ) as f:
            MassTable = f["Header"].attrs["MassTable"]
            
        cmx, cmy, cmz = il.groupcat.loadSingle(self.basePath, self.snapnum, subhaloID=self.subhaloID)["SubhaloCM"]
        self.DATA["cmx"] = cmx
        self.DATA["cmy"] = cmy
        self.DATA["cmz"] = cmz
        
        self.DATA[0] = {}
        self.DATA[1] = {}
        self.DATA[4] = {}
        self.DATA[5] = {}
        
        for i in range(len(MassTable)):
            
            if i in (0, 4, 5):
                _subhaloData = il.snapshot.loadSubhalo(
                    self.basePath,
                    self.snapnum,
                    self.subhaloID,
                    i,
                    fields=["Coordinates", "Masses"],
                )
                _x, _y, _z = _subhaloData["Coordinates"].T
                _m = _subhaloData["Masses"]
                
                _x, _y, _z, cmx, cmy, cmz = self.ChangeToCurrectCoordinates(_x, _y, _z, cmx, cmy, cmz)
                
                self.DATA[i]["x"] = _x
                self.DATA[i]["y"] = _y
                self.DATA[i]["z"] = _z
                self.DATA[i]["Masses"] = _m
                
                
            elif i == 1:
                _subhaloData = il.snapshot.loadSubhalo(
                    self.basePath, self.snapnum, self.subhaloID, i, fields=["Coordinates"]
                )
                _x, _y, _z = _subhaloData.T
                _m = np.full(len(_x), MassTable[i])
                
                _x, _y, _z, cmx, cmy, cmz = self.ChangeToCurrectCoordinates(_x, _y, _z, cmx, cmy, cmz)
                
                self.DATA[i]["x"] = _x
                self.DATA[i]["y"] = _y
                self.DATA[i]["z"] = _z
                self.DATA[i]["Masses"] = _m
            else:
                pass
        
        x,y,z,m = np.array([]),np.array([]),np.array([]),np.array([])
        
        for i in range(len(MassTable)):
            if i in (1, 0, 4, 5):
                _x = self.DATA[i]["x"]
                _y = self.DATA[i]["y"]
                _z = self.DATA[i]["z"]
                _m = self.DATA[i]["Masses"]
                x = np.hstack([x, _x])
                y = np.hstack([y, _y])
                z = np.hstack([z, _z])
                m = np.hstack([m, _m])
            else:
                pass
        
        self.DATA["x"] = x
        self.DATA["y"] = y
        self.DATA["z"] = z
        self.DATA["Masses"] = m
        
        return None
    
    def Get0Data(self):
        with h5py.File(
            self.basePath + f"snapdir_{self.snapnum:03}/snap_{self.snapnum:03}.0.hdf5",
            "r",
        ) as f:
            MassTable = f["Header"].attrs["MassTable"]
            
        cmx, cmy, cmz = il.groupcat.loadSingle(self.basePath, self.snapnum, subhaloID=self.subhaloID)["SubhaloCM"]
        self.DATA["cmx"] = cmx
        self.DATA["cmy"] = cmy
        self.DATA["cmz"] = cmz
        
        self.DATA[0] = {}
        i = 0
        _subhaloData = il.snapshot.loadSubhalo(
            self.basePath,
            self.snapnum,
            self.subhaloID,
            i,
            fields=["Coordinates", "Masses"],
        )
        _x, _y, _z = _subhaloData["Coordinates"].T
        _m = _subhaloData["Masses"]
        
        _x, _y, _z, cmx, cmy, cmz = self.ChangeToCurrectCoordinates(_x, _y, _z, cmx, cmy, cmz)
        
        self.DATA[i]["x"] = _x
        self.DATA[i]["y"] = _y
        self.DATA[i]["z"] = _z
        self.DATA[i]["Masses"] = _m
        
        return None
    
    
    def GetVirialRadius(self):
        df = pd.read_csv("../VirialTable.csv")
        return df.query("SubhaloID == @self.subhaloID")["VirialRadius"]
    
    def GetSubhaloID(self, field, parttype):
        if isinstance(field, list):
            for f in field:
                self.DATA[parttype][f] = il.snapshot.loadSubhalo(
                    self.basePath,
                    self.snapnum,
                    self.subhaloID,
                    parttype,
                    fields=f
                )
        else:
            self.DATA[parttype][field] = il.snapshot.loadSubhalo(
                self.basePath,
                self.snapnum,
                self.subhaloID,
                parttype,
                fields=field
            )
        return None
    
    def calc_rot(self):
        x = self.DATA["x"]
        y = self.DATA["y"]
        z = self.DATA["z"]
        m = self.DATA["Masses"]
        # Calculate the moment of inertia tensor
        I = np.zeros((3, 3))
        I[0][0] = np.sum(m * (y**2 + z**2))
        I[1][1] = np.sum(m * (x**2 + z**2))
        I[2][2] = np.sum(m * (x**2 + y**2))
        I[0][1] = -1 * np.sum(m * (x * y))
        I[0][2] = -1 * np.sum(m * (x * z))
        I[1][2] = -1 * np.sum(m * (y * z))
        I[1][0] = I[0][1]
        I[2][0] = I[0][2]
        I[2][1] = I[1][2]

        # Solve the eigenvalue equation and derive the rotation matrix
        eigen_values, rotation_matrix = np.linalg.eig(I)
        sort_inds = np.argsort(eigen_values)
        eigen_values = eigen_values[sort_inds]
        R = np.array(
            (
                rotation_matrix[sort_inds[0]],
                rotation_matrix[sort_inds[1]],
                rotation_matrix[sort_inds[2]],
            )
        )
        return R
    
    def calc_rot2(self,parttype):
        x = self.DATA[parttype]["x"]
        y = self.DATA[parttype]["y"]
        z = self.DATA[parttype]["z"]
        m = self.DATA[parttype]["Masses"]
        # Calculate the moment of inertia tensor
        I = np.zeros((3, 3))
        I[0][0] = np.sum(m * (y**2 + z**2))
        I[1][1] = np.sum(m * (x**2 + z**2))
        I[2][2] = np.sum(m * (x**2 + y**2))
        I[0][1] = -1 * np.sum(m * (x * y))
        I[0][2] = -1 * np.sum(m * (x * z))
        I[1][2] = -1 * np.sum(m * (y * z))
        I[1][0] = I[0][1]
        I[2][0] = I[0][2]
        I[2][1] = I[1][2]

        # Solve the eigenvalue equation and derive the rotation matrix
        eigen_values, rotation_matrix = np.linalg.eig(I)
        sort_inds = np.argsort(eigen_values)
        eigen_values = eigen_values[sort_inds]
        R = np.array(
            (
                rotation_matrix[sort_inds[0]],
                rotation_matrix[sort_inds[1]],
                rotation_matrix[sort_inds[2]],
            )
        )
        return R
    
    def calc_rot_x(self, theta):
        theta = np.radians(theta)
        return np.array(
            [
                [1, 0, 0],
                [0, np.cos(theta), -np.sin(theta)],
                [0, np.sin(theta), np.cos(theta)],
            ]
        )

    def calc_rot_y(self, theta):
        theta = np.radians(theta)
        return np.array(
            [
                [np.cos(theta), 0, np.sin(theta)],
                [0, 1, 0],
                [-np.sin(theta), 0, np.cos(theta)],
            ]
        )

    def calc_rot_z(self, theta):
        theta = np.radians(theta)
        return np.array(
            [
                [np.cos(theta), -np.sin(theta), 0],
                [np.sin(theta), np.cos(theta), 0],
                [0, 0, 1],
            ]
        )
    
    def import_abund(self,_print=True):
        abund_table = pd.read_csv("/home/nishihama/test_halo/solar_xspec_qiita.csv")
        if _print:
            print("[DONE] IMPORTED XSPEC Solar abundance table")
            print(f"[DONE] SET \033[31m{self.abund}\033[0m Solar abundance table")
        self.abund_table = {}
        for i in range(len(abund_table)):
            self.abund_table[abund_table["El"][i]] = abund_table[f"{self.abund}"][i]
        if _print:
            print(f"[DONE] IMPORTED \033[31m{self.abund}\033[0m Solar abundance table")
        self.atomicMass = {}
        for i in range(len(abund_table)):
            self.atomicMass[abund_table["El"][i]] = abund_table["atomicMass"][i]
        if _print:
            print(f"[DONE] CREATED atomicMass list")
            print(f"[DONE] All process is done.")
        
    def calc_abund(self):
        if "GFM_Metals" in self.DATA[0].keys():
            self.DATA[0]["abund"] = []
            h = self.DATA[0]["GFM_Metals"].T[0] / self.atomicMass[self.Metals[0]] / self.abund_table[self.Metals[0]]
            
            for i in range(9):
                fe = self.DATA[0]["GFM_Metals"].T[i] / self.atomicMass[self.Metals[i]] / self.abund_table[self.Metals[i]]
                
                self.DATA[0]["abund"].append(
                    fe / h
                )
    def plot(self,bins=301,parttype=0):
        plt.figure(figsize=(10,10),dpi=400)
        x = self.DATA[parttype]["x"]
        y = self.DATA[parttype]["y"]
        m = self.DATA[parttype]["Masses"]
        xbins = np.linspace(x.min(), x.max(), bins)
        ybins = np.linspace(y.min(), y.max(), bins)
        hist, xed, yed = np.histogram2d(x, y, bins=(xbins, ybins), weights = m)
        hist[hist == 0] = hist[hist > 0].min()
        plt.pcolormesh(xed, yed, hist.T, norm=LogNorm())
        plt.colorbar()
        plt.show()
        
    def plot2(self,bins=301,parttype=0):
        plt.figure(figsize=(10,10),dpi=400)
        x = self.DATA[parttype]["x"]
        y = self.DATA[parttype]["y"]
        z = self.DATA[parttype]["z"]
        
        x,y,z = self.calc_rot2(parttype) @ np.array([x,y,z])
        
        m = self.DATA[parttype]["Masses"]
        xbins = np.linspace(x.min(), x.max(), bins)
        ybins = np.linspace(y.min(), y.max(), bins)
        hist, xed, yed = np.histogram2d(x, y, bins=(xbins, ybins), weights = m)
        hist[hist == 0] = hist[hist > 0].min()
        plt.pcolormesh(xed, yed, hist.T, norm=LogNorm())
        plt.colorbar()
        plt.show()
        
    def Get4Data(self):
        with h5py.File(
            self.basePath + f"snapdir_{self.snapnum:03}/snap_{self.snapnum:03}.0.hdf5",
            "r",
        ) as f:
            MassTable = f["Header"].attrs["MassTable"]
            
        cmx, cmy, cmz = il.groupcat.loadSingle(self.basePath, self.snapnum, subhaloID=self.subhaloID)["SubhaloCM"]
        self.DATA["cmx"] = cmx
        self.DATA["cmy"] = cmy
        self.DATA["cmz"] = cmz
        
        self.DATA[4] = {}        
        i = 4
        _subhaloData = il.snapshot.loadSubhalo(
            self.basePath,
            self.snapnum,
            self.subhaloID,
            i,
            fields=["Coordinates", "Masses"],
        )
        _x, _y, _z = _subhaloData["Coordinates"].T
        _m = _subhaloData["Masses"]
        
        _x, _y, _z, cmx, cmy, cmz = self.ChangeToCurrectCoordinates(_x, _y, _z, cmx, cmy, cmz)
        
        self.DATA[i]["x"] = _x
        self.DATA[i]["y"] = _y
        self.DATA[i]["z"] = _z
        self.DATA[i]["Masses"] = _m
        
        x,y,z,m = np.array([]),np.array([]),np.array([]),np.array([])
        
        return None
    
    def calc_temp(
        self,
        gamma = 5.0/3.0,
        KB = 1.3807e-16,
        mp = 1.6726e-24,
        unit = "K"
        ):
        
        u = self.DATA[0]["InternalEnergy"]
        E = self.DATA[0]["ElectronAbundance"]
        mu = 4/(1 + 3*0.76 + 4*0.76*E) * mp
        temperature = (gamma-1)* (u/KB)* mu* 1e10
        
        if unit == "K":
            return temperature
        elif unit == "keV":
            return temperature * 8.61732814974493e-8