import math
from pathlib import Path
from time import perf_counter

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from scipy.sparse import lil_matrix, eye
from scipy.sparse.linalg import factorized
from matplotlib.path import Path as MplPath

vertices_octagon = np.array([
    [2, 4], [4, 2], [4, -2], [2, -4],
    [-2, -4], [-4, -2], [-4, 2], [-2, 4], [2, 4]
])
vertices = vertices_octagon  

x_min, x_max = -10.0, 10.0
y_min, y_max = -10.0, 10.0
h = 0.5 

x_coords = np.arange(x_min, x_max + h, h)
y_coords = np.arange(y_min, y_max + h, h)

Dx = 7.0
Dy = 6.0
vx = 0.5
vy = -0.3

source_func = lambda x, y: 2.0 * np.sin(x * 0.5) + 3.0 * np.cos(y * 0.5)

# (t = 0)
initial_cond = 0.0

final_time = 2.0
dt = 0.05

def calculate_polygon_system(
    vertices, x_coords, y_coords, Dx, Dy, vx, vy, source
):
    """
    Збирає провідний оператор для рівняння дифузії-адвекції всередині довільного полігону.
    """
    vertices = np.asarray(vertices, dtype=float)
    X, Y = np.meshgrid(x_coords, y_coords)
    Ny, Nx = X.shape
    
    polygon_path = MplPath(vertices)
    points = np.column_stack((X.ravel(), Y.ravel()))
    inside = polygon_path.contains_points(points).reshape(X.shape)
    
    n_unknowns = int(np.count_nonzero(inside))
    if n_unknowns == 0:
        raise ValueError("Сітка не містить внутрішніх точок фігури.")

    node_map = np.full(X.shape, -1, dtype=int)
    node_map[inside] = np.arange(n_unknowns)

    if callable(source):
        forcing = np.asarray(source(X[inside], Y[inside]), dtype=float).reshape(-1)
    else:
        forcing = np.full(n_unknowns, source, dtype=float)

    h_x = h_y = h
    left = Dx / h_x**2 - vx / (2 * h_x)
    right = Dx / h_x**2 + vx / (2 * h_x)
    down = Dy / h_y**2 - vy / (2 * h_y)
    up = Dy / h_y**2 + vy / (2 * h_y)

    directions = (
        (0, -1, left, np.array([-h_x, 0.0])),
        (0, 1, right, np.array([h_x, 0.0])),
        (-1, 0, down, np.array([0.0, -h_y])),
        (1, 0, up, np.array([0.0, h_y])),
    )

    operator = lil_matrix((n_unknowns, n_unknowns), dtype=float)
    
    for row, col in np.argwhere(inside):
        k = node_map[row, col]
        operator[k, k] = -2.0 * (Dx / h_x**2 + Dy / h_y**2)

        for row_step, col_step, coefficient, _ in directions:
            nr, nc = row + row_step, col + col_step
            if 0 <= nr < Ny and 0 <= nc < Nx and inside[nr, nc]:
                neighbor_index = node_map[nr, nc]
                operator[k, neighbor_index] += coefficient
            else:
                # Граничні умови Діріхле для x < 0, Нейман для x >= 0
                if X[row, col] < 0:
                    # Умова Діріхле 
                    dirichlet_val = X[row, col]**2 + abs(Y[row, col])
                    forcing[k] -= coefficient * dirichlet_val
                else:
                    # Умова Неймана
                    operator[k, k] += coefficient

    return operator.tocsc(), forcing, X, Y, inside, Ny, Nx


def solve_polygon_over_time(
    vertices, x_coords, y_coords, Dx, Dy, vx, vy, source, initial_condition, final_time, dt
):
    """
    Розв'язує нестаціонарну задачу за допомогою неявної схеми Ейлера 
    """
    operator, forcing, X, Y, inside, Ny, Nx = calculate_polygon_system(
        vertices, x_coords, y_coords, Dx, Dy, vx, vy, source
    )
    
    number_of_steps = math.ceil(final_time / dt)
    step_size = final_time / number_of_steps
    
    step_matrix = eye(operator.shape[0], format="csc") - step_size * operator
    solve_step = factorized(step_matrix)

    if callable(initial_condition):
        solution = np.asarray(initial_condition(X[inside], Y[inside]), dtype=float).reshape(-1)
    else:
        solution = np.full(np.count_nonzero(inside), initial_condition, dtype=float)

    for step in range(1, number_of_steps + 1):
        solution = solve_step(solution + step_size * forcing)

    return solution, X, Y, inside, Ny, Nx

start_time = perf_counter()
solution, X, Y, inside, Ny, Nx = solve_polygon_over_time(
    vertices=vertices,
    x_coords=x_coords,
    y_coords=y_coords,
    Dx=Dx,
    Dy=Dy,
    vx=vx,
    vy=vy,
    source=source_func,
    initial_condition=initial_cond,
    final_time=final_time,
    dt=dt,
)
solve_seconds = perf_counter() - start_time
print(f"Розрахунок завершено за {solve_seconds:.3f} секунд.")

u_full = np.full(X.shape, np.nan)
u_full[inside] = solution


np.set_printoptions(threshold=np.inf, linewidth=150, suppress=True)

print("=" * 70)
print(f" МАТРИЦЯ ЗНАЧЕНЬ КОНЦЕНТРАЦІЇ u (Розмір сітки: {Nx}x{Ny}, крок h = {h}):")
print("=" * 70)

matrix_table = np.where(np.isnan(u_full), 0.0, np.round(u_full, 1))
print(matrix_table[::-1, :]) 
print("=" * 70)


u_plot = np.where(np.isnan(u_full), 0.0, u_full)

fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(X, Y, u_plot, cmap=cm.coolwarm, edgecolor='none', alpha=0.9)
fig.colorbar(surf, ax=ax, shrink=0.5, aspect=5, label='Концентрація (u)')
ax.set_title(f'модель (h = {h})', fontsize=12)
ax.set_xlabel('Вісь X', fontsize=10)
ax.set_ylabel('Вісь Y', fontsize=10)
ax.set_zlabel('u', fontsize=10)
ax.view_init(elev=30, azim=-60)
plt.show()

plt.figure(figsize=(8, 6))
plt.contourf(X, Y, u_plot, levels=20, cmap='coolwarm')
plt.colorbar(label='Концентрація (u)')
plt.plot(vertices[:, 0], vertices[:, 1], 'k-', linewidth=2, label='Контур фігури')
plt.title(f'Плоска карта розподілу (t = {final_time}, h = {h})', fontsize=12)
plt.xlabel('Вісь X')
plt.ylabel('Вісь Y')
plt.grid(True, linestyle=':', alpha=0.6)
plt.legend()
plt.axis('equal')
plt.show()
