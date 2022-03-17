import numpy as np
import matplotlib.pyplot as plt
import esv
import msv
from numpy import sin,cos,tan


GLOBAL_STEP = 0.5
GOAL_NUMBER = 2

iter = int(GOAL_NUMBER/GLOBAL_STEP) # Number of iterations

def euler_function(val):
    return np.exp(val) 

def euler_f(val):
    return val

def log_function(val):
    return np.log(val+1)

def log_ableitung(val):
    return 1/(val+1)

def sin_func(val):
    return sin(val)

def cos_func(val):
    return cos(val)

def tang_func(val):
    return tan(val)

def tang_func_r(val):
    return 1 + pow(val,2)

ziel_func = log_function

ableitung = log_ableitung

START_VALUE = ziel_func(0)

if __name__ == "__main__":
    x = np.linspace(0,GOAL_NUMBER,200)
    plt.plot(x,ziel_func(x),label="Zielfunktion")

    euler_res = esv.generelle_einschritt_verfahren(number_of_iterations=iter,start_value=START_VALUE,step=GLOBAL_STEP,func=ableitung, verfahren=esv.Euler_verfahren_r)
    imp_euler_res = esv.generelle_einschritt_verfahren(number_of_iterations=iter,start_value=START_VALUE,step=GLOBAL_STEP,func=ableitung, verfahren=esv.verbessertes_Euler_verfahren_r)
    heun_res = esv.generelle_einschritt_verfahren(number_of_iterations=iter,start_value=START_VALUE,step=GLOBAL_STEP,func=ableitung, verfahren=esv.Heun_verfahren_r)
    msv_res = msv.zwei_schritt_Adams_Bashforth_verfahren(number_of_iterations=iter,start_value=START_VALUE,step=GLOBAL_STEP,func=ableitung)
    plt.plot(np.linspace(0,GOAL_NUMBER,iter+1),euler_res,label="Euler Verfahren")
    plt.plot(np.linspace(0,GOAL_NUMBER,iter+1),heun_res,label="Heun Verfahren")
    print(heun_res,imp_euler_res)
    plt.plot(np.linspace(0,GOAL_NUMBER,iter+1),imp_euler_res,label="Verbessertes Euler Verfahren")
    #plt.plot(np.linspace(0,GOAL_NUMBER,iter+1),msv_res,label="Zwei-Schritt Adams Bashforth verfahren")
    
    plt.legend()
    plt.grid()
    plt.show()