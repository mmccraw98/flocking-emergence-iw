import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sklearn.decomposition import PCA
import pandas as pd
from utils import get_velocity_order_parameter

# Assuming df_A is already defined and it has columns 'x', 'y', 'z', and 'time'

df = pd.read_csv('A_nodes.csv', index_col=0)

sdf = df[df.time == np.sort(df.time.unique())[700]]

order_parameter = get_velocity_order_parameter(sdf)

pca = PCA(n_components=3)
pca.fit(sdf[['x', 'y', 'z']])
v1 = pca.components_[0]  # First PCA axis
v2 = pca.components_[1]  # Second PCA axis
v3 = pca.components_[2]  # Third PCA axis (for offset)

normal = np.cross(v1, v2)
point = np.mean(sdf[['x', 'y', 'z']].values, axis=0)
d = -point.dot(normal)

xx, yy = np.meshgrid(np.linspace(min(sdf.x), max(sdf.x), 10), np.linspace(min(sdf.y), max(sdf.y), 10))
zz = (-normal[0] * xx - normal[1] * yy - d) / normal[2]

# Calculate offset for the third PCA axis
offset = np.linalg.norm(v3) * pca.explained_variance_[2] / 5

fig, ax = plt.subplots(subplot_kw={'projection': '3d'})
print(order_parameter)
plt.title(str(order_parameter.values[0]) + ' ' + str(pca.explained_variance_ratio_[2]))
# Plot the planes (one above and one below the swarm)
ax.plot_surface(xx, yy, zz + offset, alpha=0.5, rstride=100, cstride=100, color='lightblue')
ax.plot_surface(xx, yy, zz - offset, alpha=0.5, rstride=100, cstride=100, color='lightblue')

# Highlight the original vectors
ax.quiver(0, 0, 0, v1[0], v1[1], v1[2], color='r')
ax.quiver(0, 0, 0, v2[0], v2[1], v2[2], color='b')

# Plot the swarm
ax.scatter(sdf.x, sdf.y, sdf.z, s=0.1)

# Set axes to be equal
max_range = np.array([sdf.x.max() - sdf.x.min(), sdf.y.max() - sdf.y.min(), sdf.z.max() - sdf.z.min()]).max() / 2.0
mid_x = (sdf.x.max() + sdf.x.min()) * 0.5
mid_y = (sdf.y.max() + sdf.y.min()) * 0.5
mid_z = (sdf.z.max() + sdf.z.min()) * 0.5

ax.set_xlim(mid_x - max_range, mid_x + max_range)
ax.set_ylim(mid_y - max_range, mid_y + max_range)
ax.set_zlim(mid_z - max_range, mid_z + max_range)

plt.show()