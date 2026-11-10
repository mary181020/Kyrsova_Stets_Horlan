import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.cm as cm

vertices_octagon = np.array([
    [2, 4],
    [4, 2],
    [4, -2],
    [2, -4],
    [-2, -4],
    [-4, -2],
    [-4, 2],
    [-2, 4],
    [2, 4]  
])

vertices = vertices_octagon  
polygon_path = Path(vertices)

x_min, x_max = -10.0, 10.0
y_min, y_max = -10.0, 10.0
h = 0.5 

x = np.arange(x_min, x_max + h, h)
y = np.arange(y_min, y_max + h, h)
X, Y = np.meshgrid(x, y)
Ny, Nx = X.shape

points = np.column_stack((X.ravel(), Y.ravel()))
mask_flat = polygon_path.contains_points(points)
mask = mask_flat.reshape(X.shape)

boundary_mask = np.zeros((Ny, Nx), dtype=bool)
for j in range(1, Ny - 1):
    for i in range(1, Nx - 1):
        if mask[j, i]:  
            if not mask[j, i-1] or not mask[j, i+1] or not mask[j-1, i] or not mask[j+1, i]:
                boundary_mask[j, i] = True  


u = np.zeros((Ny, Nx))

source_matrix = np.full((Ny, Nx), 2.0)

Dx = 7
Dy = 6
vx = 0.5
vy = -0.3

tau = 0.005      
num_steps = 400   

for step in range(num_steps):
    u_new = u.copy()
    
    for j in range(1, Ny - 1):
        for i in range(1, Nx - 1):
            
            if not mask[j, i]:
                continue
                
            if boundary_mask[j, i]:
                if X[j, i] < 0:
                    u_new[j, i] = X[j, i]**2 + abs(Y[j, i])  # Діріхле 
                else:
                    u_new[j, i] = u[j, i-1]  # Нейман 
                continue
               
            # Сусіди всередині фігури
            left  = u[j, i-1]
            right = u[j, i+1]
            bottom = u[j-1, i]
            top = u[j+1, i]
            center = u[j, i]
            
            d2x = (right - 2.0 * center + left) / (h**2)
            d2y = (top - 2.0 * center + bottom) / (h**2)
            
            dx_conv = (right - left) / (2.0 * h)
            dy_conv = (top - bottom) / (2.0 * h)
            
            # Швидкість зміни концентрації в цій точці (права частина рівняння)
            # Баланс = Дифузія + Конвекція + Джерело(x, y)
            change_rate = (Dx * d2x + Dy * d2y) - (vx * dx_conv + vy * dy_conv) + source_matrix[j, i]
            
            # Явна схема Ейлера за часом: нове = старе + tau * швидкість зміни
            u_new[j, i] = center + tau * change_rate
            
    u = u_new

# ==========================================
# КРОК 4: ВИВІД ТАБЛИЦІ ТА ГРАФІКІВ
# ==========================================
np.set_printoptions(threshold=np.inf, linewidth=150, suppress=True)

print("=" * 70)
print(f" МАТРИЦЯ ЗНАЧЕНЬ КОНЦЕНТРАЦІЇ u ПІСЛЯ {num_steps} КРОКІВ ЧАСУ (крок h = {h}):")
print("=" * 70)

matrix_table = np.where(mask, np.round(u, 1), 0.0)
print(matrix_table[::-1, :])  
print("=" * 70)

u_plot = np.where(mask, u, 0.0)

# 3D Візуалізація
fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(X, Y, u_plot, cmap=cm.coolwarm, edgecolor='none', alpha=0.9)
fig.colorbar(surf, ax=ax, shrink=0.5, aspect=5, label='Концентрація (u)')
ax.set_title(f'Нестаціонарна модель (Часових кроків: {num_steps}, h = {h})', fontsize=12)
ax.set_xlabel('Вісь X', fontsize=10)
ax.set_ylabel('Вісь Y', fontsize=10)
ax.set_zlabel('u', fontsize=10)
ax.view_init(elev=30, azim=-60)
plt.show()

# Плоска карта розподілу
plt.figure(figsize=(8, 6))
plt.contourf(X, Y, u_plot, levels=20, cmap='coolwarm')
plt.colorbar(label='Концентрація (u)')
plt.plot(vertices[:, 0], vertices[:, 1], 'k-', linewidth=2, label='Контур фігури')
plt.title(f'Плоска карта розподілу (t > 0, h = {h})', fontsize=12)
plt.xlabel('Вісь X')
plt.ylabel('Вісь Y')
plt.grid(True, linestyle=':', alpha=0.6)
plt.legend()
plt.axis('equal')
plt.show()

"""import numpy as npimport matplotlib.pyplot as pltfrom matplotlib.path import Pathimport matplotlib.cm as cm# Зірка

vertices_star = np.array([

[0, 4],

[1, 1],

[4, 1],

[2, -1],

[3, -4],

[0, -2],

[-3, -4],

[-2, -1],

[-4, 1],

[-1, 1],

[0, 4]

])# Восьмикутник

vertices_octagon = np.array([

[2, 4],

[4, 2],

[4, -2],

[2, -4],

[-2, -4],

[-4, -2],

[-4, 2],

[-2, 4],

[2, 4]

])# квадрат

vertices_box = np.array([

[-3, -3],

[3, -3],

[3, 3],

[-3, 3],

[-3, -3]

])

vertices_circle = np.array([

[4, 0], [3, 3], [0, 4], [-3, 3],

[-4, 0], [-3, -3], [0, -4], [3, -3], [4, 0]

])

vertices_random = np.array([

[-4, 0], [-2, 3], [1, 2], [3, 4],

[4, -1], [2, -3], [-1, -2], [-3, -3], [-4, 0]

])



vertices = vertices_random



polygon_path = Path(vertices)#зшиваємо точки у єдиний замкнений контур



x_min, x_max = -10.0, 10.0

y_min, y_max = -10.0, 10.0

h = 0.25



x = np.arange(x_min, x_max + h, h)#створює одновимірний масив від мінімального значення до максимального з певним кроком

y = np.arange(y_min, y_max + h, h)

X, Y = np.meshgrid(x, y)#meshgrid перетворює дві окремі лінійки (x та y) на дві величезні матриці, які разом утворюють повноцінну сітку

Ny, Nx = X.shape



points = np.column_stack((X.ravel(), Y.ravel()))#Метод .ravel() бере двовимірну матрицю (нашу сітку) і «витягує» її в один довгий-довгий одновимірний список точок З'єднує ці два довгі списки в одну таблицю координат пар $(x, y)$, де в першому стовпчику йдуть ікси, а в другому — ігрики. # Перевіряємо кожну точку

mask_flat = polygon_path.contains_points(points)

mask = mask_flat.reshape(X.shape)# Межа — це такі точки (True), які мають хоча б одного сусіда зі значенням False (ззовні)

boundary_mask = np.zeros((Ny, Nx), dtype=bool)for j in range(1, Ny - 1):

for i in range(1, Nx - 1):

if mask[j, i]:


if not mask[j, i-1] or not mask[j, i+1] or not mask[j-1, i] or not mask[j+1, i]:

boundary_mask[j, i] = True



u = np.zeros((Ny, Nx))



source_term = -50.0 # Права частина рівняння (= 5)# Коефіцієнти

Dx = 1.2

Dy = 0.8

vx = 0.5

vy = -0.3



iterations = 2000for step in range(iterations):

u_new = u.copy()


for j in range(1, Ny - 1):

for i in range(1, Nx - 1):


# Якщо точка ззовні — пропускаємо

if not mask[j, i]:

continue


# Якщо точка на МЕЖІ: розділяємо умови Діріхле та Неймана

if boundary_mask[j, i]:

# Приклад розділення: якщо координата X < 0, то це умова Діріхле

if X[j, i] < 0:

u_new[j, i] = X[j, i]**2 + abs(Y[j, i]) # Умова Діріхле (функція)

else:

# Умова Неймана (нульовий потік: беремо значення з сусідньої внутрішньої точки)

u_new[j, i] = u[j, i-1]

continue


left = u[j, i-1]

right = u[j, i+1]

bottom = u[j-1, i]

top = u[j+1, i]


# Коефіцієнти для знаменника (щоб усе ділилося правильно)

coeff_x = Dx / (h**2)

coeff_y = Dy / (h**2)


# Вплив конвекції (перші похідні)

convection_x = vx * (right - left) / (2.0 * h)

convection_y = vy * (top - bottom) / (2.0 * h)


# Знаменник (центральний коефіцієнт)

denominator = 2.0 * (coeff_x + coeff_y)


# Збираємо чисельник (дифузія + конвекція + джерело)

numerator = (coeff_x * (right + left) +

coeff_y * (top + bottom) -

convection_x - convection_y + source_term)


u_new[j, i] = numerator / denominator


u = u_new





np.set_printoptions(threshold=np.inf, linewidth=150, suppress=True)print("="70)
print(f" МАТРИЦЯ ЗНАЧЕНЬ КОНЦЕНТРАЦІЇ u (Розмір сітки: {Nx}x{Ny}, крок h = {h}):")print("="70)


matrix_table = np.where(mask, np.round(u, 1), 0.0)print(matrix_table[::-1, :]) print("="70)


u_plot = np.where(mask, u, 0.0)



plt.figure(figsize=(8, 6))

plt.contourf(X, Y, u_plot, levels=20, cmap='coolwarm')

plt.colorbar(label='Концентрація (u)')

plt.plot(vertices[:, 0], vertices[:, 1], 'k-', linewidth=2, label='Контур фігури')

plt.title(f'Плоска карта розподілу (h = {h})', fontsize=12)

plt.xlabel('Вісь X')

plt.ylabel('Вісь Y')

plt.grid(True, linestyle=':', alpha=0.6)

plt.legend()

plt.axis('equal')

plt.show()  """