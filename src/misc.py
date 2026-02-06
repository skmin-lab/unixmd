from functools import wraps
import sys, time, os
import numpy as np

# Atomic weight
data = { "xx" : 1.00794, "H" : 1.00794, "He" : 4.00260, "Li" : 6.941, "Be" : 9.012187, "B" : 10.811,
        "C" : 12.0107, "N" : 14.00674, "O" : 15.9994, "F" : 18.99840, "Ne" : 20.1797, "Na" : 22.98977,
        "Mg" : 24.3050, "Al" : 26.98152, "Si" : 28.0855, "P" : 30.97376, "S" : 32.066, "Cl" : 35.4527,
        "Ar" : 39.948, "K" : 39.0983, "Ca" : 40.078, "Sc" : 44.95591, "Ti" : 47.867, "V" : 50.9415,
        "Cr" : 51.9961, "Mn" : 54.93805, "Fe" : 55.845, "Co" : 58.93320, "Ni" : 58.6934, "Cu" : 63.546,
        "Zn" : 65.39, "Ga" : 69.723, "Ge" : 72.61, "As" : 74.92160, "Se" : 78.96, "Br" : 79.904,
        "Kr" : 83.80, "Rb" : 85.4678, "Sr" : 87.62, "Y" : 88.90585, "Zr" : 91.224, "Nb" : 92.90638,
        "Mo" : 95.94, "Tc" : 98.0, "Ru" : 101.07, "Rh" : 102.90550, "Pd" : 106.42, "Ag" : 107.8682,
        "Cd" : 112.411, "In" : 114.818, "Sn" : 118.710, "Sb" : 121.760, "Te" : 127.60, "I" : 126.90477,
        "Xe" : 131.29, "Cs" : 132.90545, "Ba" : 137.327, "La" : 138.9055, "Ce" : 140.116, "Pr" : 140.90765,
        "Nd" : 144.24, "Pm" : 145.0, "Sm" : 150.36, "Eu" : 151.964, "Gd" : 157.24, "Tb" : 158.92534,
        "Dy" : 162.50, "Ho" : 164.93032, "Er" : 167.26, "Tm" : 168.93421, "Yb" : 173.04, "Lu" : 174.967,
        "Hf" : 178.49, "Ta" : 180.9479, "W" : 183.84, "Re" : 186.207, "Os" : 190.23, "Ir" : 192.217,
        "Pt" : 195.078, "Au" : 196.96655, "Hg" : 200.59, "Tl" : 204.3833, "Pb" : 207.2, "Bi" :208.98038,
        "Po" : 209.0, "At" : 210.0, "Rn" : 222.0, "Fr" :223.0, "Ra" : 226.0, "Ac" : 227.0,
        "Th" : 232.0381, "Pa" : 231.03588, "U" : 238.0289, "Np" : 237.0, "Pu" : 244.0, "Am" : 243.0,
        "Cm" : 247.0, "Bk" : 247.0, "Cf" : 251.0, "Es" : 252.0, "Fm" : 257.0, "Md" : 258.0,
        "No" : 259.0, "Lr" : 262.0, "Rf" : 261.0, "Db" : 262.0, "Sg" : 263.0, "Bh" : 264.0,
        "Hs" : 265.0, "Mt" : 268.0, "Ds" : 271.0, "Rg" : 272.0, "Uub" : 285.0, "Uut" : 284.0,
        "Uuq" : 289.0, "Uup" : 288.0, "Uuh" : 292.0}

# Conversion units
amu_to_au = 1822.888486192

au_to_A = 0.529177249
A_to_au = 1 / au_to_A

au_to_fs = 0.02418884344
fs_to_au = 1 / au_to_fs # = 41.34137304

au_to_K = 3.15774646E+5

au_to_eV = 27.2113961
eV_to_au = 1 / au_to_eV # = 0.03674931
au_to_kcalmol = 627.503
kcalmol_to_au = 1 / au_to_kcalmol # = 0.00159362

# Speed of light in atomic unit
c = 137.035999108108 

# Frequency unit
cm_to_au = 1.0E-8 * au_to_A * c

eps = 1.0E-12

data.update({n : amu_to_au * data[n] for n in data.keys()})

def elapsed_time(func):
    @wraps(func)
    def check(*args, **kwargs):
        tbegin = time.time()
        func(*args, **kwargs)
        tend = time.time()
        print (f"{func.__name__} : Elapsed time = {tend - tbegin} seconds", flush=True)
    return check

def call_name():
    return sys._getframe(1).f_code.co_name

class FileManager:
    """ Singleton class to manage persistent file handles for efficient I/O

        Instead of opening/closing files for each write operation,
        this class keeps file handles open and reuses them.
    """
    _instance = None

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
            cls._instance._handles = {}
        return cls._instance

    def write(self, string, dir_name, filename, mode):
        """ Write string to file using persistent handle

            :param string string: Text string for output file
            :param string dir_name: Directory of output file
            :param string filename: Filename of output file
            :param string mode: Fileopen mode ('w' or 'a')
        """
        tmp_name = os.path.join(dir_name, filename)

        # For write mode, close existing handle if any and create new file
        if mode == "w":
            if tmp_name in self._handles:
                self._handles[tmp_name].close()
                del self._handles[tmp_name]
            with open(tmp_name, mode) as f:
                f.write(string + "\n")
            return

        # For append mode, use persistent handle
        if tmp_name not in self._handles:
            self._handles[tmp_name] = open(tmp_name, "a")

        self._handles[tmp_name].write(string + "\n")
        self._handles[tmp_name].flush()  # Ensure data is written to disk

    def close_all(self):
        """ Close all open file handles """
        for handle in self._handles.values():
            handle.close()
        self._handles.clear()

    def close_dir(self, dir_name):
        """ Close all file handles for a specific directory

            :param string dir_name: Directory path
        """
        to_close = [path for path in self._handles if path.startswith(dir_name)]
        for path in to_close:
            self._handles[path].close()
            del self._handles[path]

# Global file manager instance
_file_manager = FileManager()

def typewriter(string, dir_name, filename, mode):
    """ Function to open/write any string in dir_name/filename

        Uses persistent file handles for append mode to improve I/O performance.

        :param string string: Text string for output file
        :param string dir_name: Directory of output file
        :param string filename: Filename of output file
        :param string mode: Fileopen mode
    """
    _file_manager.write(string, dir_name, filename, mode)

def close_files(dir_name=None):
    """ Close file handles managed by the FileManager

        :param string dir_name: If provided, close only files in this directory.
                                If None, close all files.
    """
    if dir_name is None:
        _file_manager.close_all()
    else:
        _file_manager.close_dir(dir_name)

def gaussian1d(x, const, sigma, x0):
    if (sigma < 0.0):
        return -1
    else:
        res = const / (sigma * np.sqrt(2. * np.pi)) * np.exp(- (x - x0) ** 2 / (2. * sigma ** 2))
        return res 
