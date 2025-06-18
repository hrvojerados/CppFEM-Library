#!/usr/bin/env python3
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from mpl_toolkits.mplot3d import Axes3D

def main():
    header = sys.stdin.readline().split()
    Nn, Nt = map(int, header)
    node_data = np.loadtxt(sys.stdin, max_rows=Nn)
    x_coords = node_data[:,0]
    y_coords = node_data[:,1]
    sol      = node_data[:,2]

    triangles = np.loadtxt(sys.stdin, dtype=int)

    triangulation = tri.Triangulation(x_coords, y_coords, triangles)

    fig = plt.figure(figsize=(10, 8))
    ax  = fig.add_subplot(111, projection='3d')
    surf = ax.plot_trisurf(
        triangulation,
        sol,
        cmap='viridis',
        edgecolor='none'
    )
    fig.colorbar(surf, ax=ax, label='Solution')

    ax.set_title('FEM Solution - 3D Plot')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Solution Value')
    plt.show()

if __name__ == "__main__":
    main()
