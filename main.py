import numpy as np
import matplotlib.pyplot as plt

global_step = 0.5
goal_number = 2

number_of_iterations = int(goal_number/global_step)

def ziel_func(val):
    return np.exp(val)

start_value = ziel_func(0)

def global_func(val, time=1):
    return val

x = np.linspace(0,goal_number,200)
plt.plot(x,ziel_func(x))

def euler_verfahren_r(step,val):
    return val + step * global_func(val)

def euler_verfahren():
    res_list = [start_value]
    for i in range(number_of_iterations):
        res_list.append(euler_verfahren_r(global_step,res_list[i]))
    plt.plot(np.linspace(0,goal_number,number_of_iterations+1),res_list)

    
def zwei_schritt_Adams_Bashforth_verfahren_r(step,val_1,val_2):
    return val_1 + 3/2 * step * global_func(val_1) - 1/2 * step * global_func(val_2)

def zwei_schritt_Adams_Bashforth_verfahren():
    res_list = [start_value,euler_verfahren_r(global_step,start_value)]
    for i in range(number_of_iterations-1):
        res_list.append(zwei_schritt_Adams_Bashforth_verfahren_r(global_step,res_list[i+1],res_list[i]))
        print(global_step,res_list[i+1],res_list[i])
        print(res_list)
    plt.plot(np.linspace(0,goal_number,number_of_iterations+1),res_list)


euler_verfahren()
zwei_schritt_Adams_Bashforth_verfahren()
plt.show()