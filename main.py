import numpy as np
import matplotlib.pyplot as plt
import esv
import msv



GLOBAL_STEP = 0.5
GOAL_NUMBER = 10

NUMBER_OF_ITERATIONS = int(GOAL_NUMBER/GLOBAL_STEP)

def ziel_func(val):
    return np.exp(val)

START_VALUE = ziel_func(0)

def global_func(val, time=1):
    return val

if __name__ == "__main__":
    x = np.linspace(0,GOAL_NUMBER,200)
    plt.plot(x,ziel_func(x),label="Zielfunktion")

    esv_res = esv.euler_verfahren(number_of_iterations=NUMBER_OF_ITERATIONS,start_value=START_VALUE,step=GLOBAL_STEP,func=global_func)
    msv_res = msv.zwei_schritt_Adams_Bashforth_verfahren(number_of_iterations=NUMBER_OF_ITERATIONS,start_value=START_VALUE,step=GLOBAL_STEP,func=global_func)
    plt.plot(np.linspace(0,GOAL_NUMBER,NUMBER_OF_ITERATIONS+1),esv_res,label="euler_verfahren")
    plt.plot(np.linspace(0,GOAL_NUMBER,NUMBER_OF_ITERATIONS+1),msv_res,label="zwei_schritt_Adams_Bashforth_verfahren")

    plt.legend()
    plt.grid()
    plt.show()