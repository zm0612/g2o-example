import numpy as np
import matplotlib.pyplot as plt

filename = './points.txt'
X, Y = [], []
with open(filename, 'r') as f:
    lines = f.readlines()
    for line in lines:
        value = [float(s) for s in line.split()]
        X.append(float(value[0]))
        Y.append(float(value[1]))

result_name = './result.txt'
with open(result_name, 'r') as r:
    lines = r.readlines()
    for line in lines:
        value = [float(s) for s in line.split()]
        a = float(value[0])
        b = float(value[1])
        my_lambda = float(value[2])

x = np.linspace(0, 10, 10)
y = a * np.exp(-my_lambda * x) + b

plt.plot(x, y, 'r')
plt.scatter(X, Y)
plt.show()
