
import numpy as np

class xyzdata:
    "Information contained in a single xyz file"
    numAtoms : int
    comments : str
    syms     : list[str]
    gem      : list[float]

    def __init__(self, n : int, c : str, s : list[str], data : list[list[float]]) -> None:
        self.numAtoms = n
        self.comments = c
        self.syms     = s
        self.gem      = data

    def __repr__(self) -> str:
        out : str = '{0:^5d}\n'.format(self.numAtoms)
        out += self.comments
        for atom_num in range(len(self.gem)):
            coord = self.gem[atom_num]
            out += '{0:^4s} {1:^14.10f} {2:^14.10f} {3:^14.10f}\n'.format(self.syms[atom_num], coord[0], coord[1], coord[2])
        return out

    def saveData(self, fileName : str) -> None:
        f = open(fileName+'.xyz', 'w')
        f.write(self.__repr__())
        f.close()
    
def read_xyz(filename : str) -> (list[str],list[list[float]]):
    atomic_symbols  = []
    xyz_coordinates = []
    num_atoms = -1
    with open(filename, "r") as file:
        for line_number, line in enumerate(file):
            if   line_number == 0: num_atoms = int(line)
            elif line_number == 1: comments  = line
            else:
                atomic_symbol, x, y, z = line.split()
                atomic_symbols.append(atomic_symbol)
                xyz_coordinates.append([float(x), float(y), float(z)])

    return (atomic_symbols, xyz_coordinates)

def read_xyz2(filename : str) -> xyzdata:
    atomic_symbols  = []
    xyz_coordinates = []
    num_atoms = -1
    with open(filename, "r") as file:
        for line_number, line in enumerate(file):
            if   line_number == 0: num_atoms = int(line)
            elif line_number == 1: comments  = line
            else:
                atomic_symbol, x, y, z = line.split()
                atomic_symbols.append(atomic_symbol)
                xyz_coordinates.append([float(x), float(y), float(z)])

    return xyzdata(num_atoms, comments, atomic_symbols, xyz_coordinates)

def rot_xyzdata(inmol : xyzdata, angle : float) -> xyzdata:

    theta : float = np.deg2rad(angle)
    
    rot = np.array([[np.cos(theta), -np.sin(theta), 0.0],
                    [np.sin(theta),  np.cos(theta), 0.0],
                    [          0.0,            0.0, 1.0]])
    
    rot_geom : List[float] = []
    for num in range(len(inmol.gem)):
        vec   = np.array(inmol.gem[num])
        rot_v = np.dot(rot,vec)
        rot_geom.append(rot_v.tolist())
        
    return xyzdata(inmol.numAtoms, inmol.comments, inmol.syms, rot_geom)

def translate(inmol : xyzdata, xdist : float, ydist : float) -> xyzdata:
    new_geom : List[float] = []
    for num in range(len(inmol.gem)):
        coord = inmol.gem[num]
        xyz = [coord[0]+xdist, coord[1]+ydist, coord[2]]
        new_geom.append(xyz)
    return xyzdata(inmol.numAtoms, inmol.comments, inmol.syms, new_geom)        
    
def overlay(mol1 : xyzdata, mol2 : xyzdata) -> xyzdata:
    newNum  = mol1.numAtoms  + mol2.numAtoms
    #newComm = mol1.comments  + mol2.comments
    newComm = 'hoge\n'
    newSyms = mol1.syms      + mol2.syms
    newGeom = mol1.gem       + mol2.gem
    return xyzdata(newNum, newComm, newSyms, newGeom)

if __name__ == '__main__':
    
#    c60 = read_xyz2('c60-cma.xyz')    
#    mc  = read_xyz2('gfn1-mc+au-au-ni-o-n-fixed.xyz')
#
#    rot_c60 = rot_xyzdata(c60, 25)
#    complex = overlay(mc, rot_c60)
#    complex.saveData('hoge')

    au       = read_xyz2('au-only-py.xyz')
    c60mc    = read_xyz2('mc+c60-py.xyz')
    
    newc60mc = translate(c60mc, 0.5, 0)
    complex  = overlay(au, newc60mc)
    complex.saveData('piyo2')
