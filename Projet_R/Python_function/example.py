from TimeWarp import TimeWarp
import numpy as np
import matplotlib.pyplot as plt


# Génère une séquence de 0 à pi avec 467 valeurs
x = np.linspace(0, 25*np.pi, 467)
cos_values = np.cos(x)

cos_values = cos_values.reshape(-1, 1)

new_array = TimeWarp(cos_values, 1000, 365, 5, True)

print(new_array)

plt.plot(new_array)
plt.show()
