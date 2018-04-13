#from eval import evaluate
import numpy as np


def hill_climbing(cost):
    epsilon = 0.00002
    dims = 10
    # initial_step_sizes = np.ones(dims) * 0.2
    initial_step_sizes = np.array([100,1000,0,0.05,0.05,0.05,0.1,0.1,0.1,0.2])
    some_acceleration = 1.2
    current_point = np.array([50, 100, -100, 0.9, 0.05, 0.4, 0.2,0.2, 0.2, 0])   # the zero-magnitude vector is common
    step_size = initial_step_sizes   # a vector of all 1's is common
    acceleration = some_acceleration  # a value such as 1.2 is common
    candidate = np.zeros(5)
    candidate[0] = -acceleration
    candidate[1] = -1 / acceleration
    candidate[2] = 0
    candidate[3] = 1 / acceleration
    candidate[4] = acceleration
    iter = 0
    while True:

        iter+=1
        print("##### iter", iter, "#####")
        before = cost(current_point)
        acc_lowered = False
        for i in range(len(current_point)):
            print("##### param", i, "#####")
            best = -1
            best_score = 100000
            for j in range(5):        # try each of 5 candidate locations
                current_point[i] = current_point[i] + step_size[i] * candidate[j]
                temp = cost(current_point)
                print(current_point)
                current_point[i] = current_point[i] - step_size[i] * candidate[j]
                if temp < best_score:
                    best_score = temp
                    best = j
            if candidate[best] == 0:
                step_size[i] = step_size[i] / acceleration
                acc_lowered = True
            else:
                current_point[i] = current_point[i] + step_size[i] * candidate[best]
                step_size[i] = step_size[i] * candidate[best]  # accelerate

        if abs(cost(current_point) - before) < epsilon and not acc_lowered:
            return current_point


def func(X):
    x = X[:,0]
    y = X[:,1]
    return (x)**2+(y)**2

def mov_func(X):
    x = X[:, 0]
    y = X[:, 1]
    return (x)**2+(y+2)**2

def rosenbrock(X):
    _a = 0
    _b = 10
    x = X[:,0]
    y = X[:,1]
    return pow(_a - y, 2) + _b * pow(y - x**2, 2)

#res = hill_climbing(mov_func)

#print(res)