import math
import meep as mp
from meep import mpb

#
num_bands = 1
k_vec_res = 16;

boundary_k_points = [mp.Vector3(),
            mp.Vector3(0.5),
            mp.Vector3(0.5, 0.5),
            mp.Vector3()]

k_points = mp.interpolate(16, boundary_k_points)

for k_vec in k_points:
    print(k_vec)

geometry = [mp.Cylinder(0.2, material=mp.Medium(epsilon=2.2))]
geometry_lattice = mp.Lattice(size=mp.Vector3(1,1))
resolution = 64


solver = mpb.ModeSolver(num_bands=num_bands,
                        k_points=k_points,
                        geometry=geometry,
                        geometry_lattice=geometry_lattice,
                        resolution=resolution)

# print_heading("Square lattice of rods: TE bands")
solver.run_te()
solver.run_tm()
