import numpy as np
import matplotlib.pyplot as plt
import esv
import msv



GLOBAL_STEP = 0.5
GOAL_NUMBER = 10

NUMBER_OF_ITERATIONS = int(GOAL_NUMBER/GLOBAL_STEP)

def euler_function(val):
    return np.exp(val) 

def euler_ableitung(val):
    return val

def log_function(val):
    return np.log(val+1)

def log_ableitung(val):
    return 1/(val+1)

ziel_func = log_function

ableitung = log_ableitung

START_VALUE = ziel_func(0)

if __name__ == "__main__":
    x = np.linspace(0,GOAL_NUMBER,200)
    plt.plot(x,ziel_func(x),label="Zielfunktion")

    esv_res = esv.euler_verfahren(number_of_iterations=NUMBER_OF_ITERATIONS,start_value=START_VALUE,step=GLOBAL_STEP,func=ableitung)
    msv_res = msv.zwei_schritt_Adams_Bashforth_verfahren(number_of_iterations=NUMBER_OF_ITERATIONS,start_value=START_VALUE,step=GLOBAL_STEP,func=ableitung)
    plt.plot(np.linspace(0,GOAL_NUMBER,NUMBER_OF_ITERATIONS+1),esv_res,label="euler_verfahren")
    plt.plot(np.linspace(0,GOAL_NUMBER,NUMBER_OF_ITERATIONS+1),msv_res,label="zwei_schritt_Adams_Bashforth_verfahren")

    plt.legend()
    plt.grid()
    plt.show()