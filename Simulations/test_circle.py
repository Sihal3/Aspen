import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

def plot_rotating_ellipse_3d(a, b, direction_vector, center_point, theta):
    # Normalize the direction vector to ensure it is a unit vector
    direction_vector = direction_vector / np.linalg.norm(direction_vector)
    
    # Convert theta from degrees to radians
    theta = np.radians(theta)
    
    # Generate points on the ellipse in 3D using parametric equations
    t = np.linspace(0, 2 * np.pi, 100)
    
    # Generate the normal vector to the plane of the ellipse
    n = np.cross(direction_vector, np.array([0, 0, 1]))
    n = n / np.linalg.norm(n)
    
    # Rotate the direction vector around its axis (theta)
    rotation_matrix = np.array([[np.cos(theta) + direction_vector[0]**2 * (1 - np.cos(theta)),
                                 direction_vector[0] * direction_vector[1] * (1 - np.cos(theta)) - direction_vector[2] * np.sin(theta),
                                 direction_vector[0] * direction_vector[2] * (1 - np.cos(theta)) + direction_vector[1] * np.sin(theta)],
                                [direction_vector[0] * direction_vector[1] * (1 - np.cos(theta)) + direction_vector[2] * np.sin(theta),
                                 np.cos(theta) + direction_vector[1]**2 * (1 - np.cos(theta)),
                                 direction_vector[1] * direction_vector[2] * (1 - np.cos(theta)) - direction_vector[0] * np.sin(theta)],
                                [direction_vector[0] * direction_vector[2] * (1 - np.cos(theta)) - direction_vector[1] * np.sin(theta),
                                 direction_vector[1] * direction_vector[2] * (1 - np.cos(theta)) + direction_vector[0] * np.sin(theta),
                                 np.cos(theta) + direction_vector[2]**2 * (1 - np.cos(theta))]])
    
    direction_vector = np.dot(rotation_matrix, direction_vector)
    
    # Generate the ellipse points in 3D space using the rotated direction vector and normal vector
    x = center_point[0] + a * (np.cos(t) * direction_vector[0] + np.sin(t) * n[0])
    y = center_point[1] + b * (np.cos(t) * direction_vector[1] + np.sin(t) * n[1])
    z = center_point[2] + a * direction_vector[2] * np.cos(t) + b * n[2] * np.sin(t)
    
    # Create a 3D plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    # Plot the ellipse
    ax.plot(x, y, z, color='b', alpha=0.6)
    
    # Set axis labels
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    
    # Set axis limits to ensure equal aspect ratio
    max_range = max(max(x) - min(x), max(y) - min(y), max(z) - min(z)) / 2.0
    mean_x = np.mean(x)
    mean_y = np.mean(y)
    mean_z = np.mean(z)
    
    ax.set_xlim(mean_x - max_range, mean_x + max_range)
    ax.set_ylim(mean_y - max_range, mean_y + max_range)
    ax.set_zlim(mean_z - max_range, mean_z + max_range)
    
    # Show the plot
    plt.show()

# Example usage
a = 5    # Radius along the major axis
b = 3    # Radius along the minor axis
direction_vector = np.array([1, 1, 1])  # Specify the direction vector here
center_point = np.array([2, 3, 4])      # Specify the center point here
theta = 45  # Rotation angle in degrees
plot_rotating_ellipse_3d(a, b, direction_vector, center_point, theta)