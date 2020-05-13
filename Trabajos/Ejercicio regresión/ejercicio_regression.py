import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('fake_regression.txt')
x_dat = data[:, 0]
y_dat = data[:, 1]
x_training = x_dat[:70]
y_training = y_dat[:70]
x_test = x_dat[70:]
y_test = y_dat[70:]

def ECP(x_training, y_training, x_test, y_test, N):
    error_training = np.array([])
    error_test = np.array([])
    for i in range(1, N+1):
        poly = np.poly1d(np.polyfit(x_training, y_training, i))
        error = sum((y_training-poly(x_training))**2)/np.size(x_training)
        error_training = np.append(error_training, error)
        error = sum((y_test-poly(x_test))**2)/np.size(x_test)
        error_test = np.append(error_test, error)
    x = np.linspace(1, N, N)
    return x, error_training, error_test

x, error_training, error_test = ECP(x_training, y_training, x_test, y_test, 20)

plt.plot(x, error_training)
plt.plot(x, error_test)
plt.savefig('test_vs_train.jpg')
plt.show()
