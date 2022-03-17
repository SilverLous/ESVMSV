import numpy as np
import matplotlib.pyplot as plt
import esv
import msv
from numpy import sin,cos,tan


GLOBAL_STEP = 0.5
GOAL_NUMBER = 5


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

def Lotka_Volterra_derivative(X, alpha, beta, delta, gamma):
    x, y = X[0],X[1]
    dotx = x * (alpha - beta * y)   # growth - predation
    doty = y * (-delta + gamma * x) # growth - mortality (or other)
    return np.array([dotx, doty])

def Lotka_temp(val):
    return Lotka_Volterra_derivative(val,1,1,1,1)

ziel_func = euler_function

ableitung = euler_f

START_VALUE = ziel_func(0) # ziel_func(0)

def normal():
    iter = int(GOAL_NUMBER/GLOBAL_STEP) # Number of iterations
    x = np.linspace(0,GOAL_NUMBER,200)
    plt.plot(x,ziel_func(x),label="Zielfunktion")

    euler_res = esv.generelle_einschritt_verfahren(      number_of_iterations=iter,start_value=START_VALUE,step=GLOBAL_STEP,func=ableitung, verfahren=esv.Euler_verfahren_r)
    imp_euler_res = esv.generelle_einschritt_verfahren(  number_of_iterations=iter,start_value=START_VALUE,step=GLOBAL_STEP,func=ableitung, verfahren=esv.verbessertes_Euler_verfahren_r)
    heun_res = esv.generelle_einschritt_verfahren(       number_of_iterations=iter,start_value=START_VALUE,step=GLOBAL_STEP,func=ableitung, verfahren=esv.Heun_verfahren_r)
    msv_res = msv.zwei_schritt_Adams_Bashforth_verfahren(number_of_iterations=iter,start_value=START_VALUE,step=GLOBAL_STEP,func=ableitung)
    plt.plot(np.linspace(0,GOAL_NUMBER,iter+1),euler_res,label="Euler")
    plt.plot(np.linspace(0,GOAL_NUMBER,iter+1),heun_res,label="Heun Verfahren")
    plt.plot(np.linspace(0,GOAL_NUMBER,iter+1),imp_euler_res,label="Verbessertes Euler Verfahren")
    plt.plot(np.linspace(0,GOAL_NUMBER,iter+1),msv_res,label="Zwei-Schritt Adams-Bashforth verfahren")
    
    plt.legend()
    plt.grid()
    plt.show()

def lotka_vol(func,abl,start,step,goal):
    p_goal = goal
    p_step = step
    p_start_value = start
    p_func = func
    p_ableitung = abl
    iter = int(goal/p_step) # Number of iterations
    x = np.linspace(0,goal,200)
    #plt.plot(x,ziel_func(x),label="Zielfunktion")

    euler_res = esv.generelle_einschritt_verfahren(number_of_iterations=iter,start_value=p_start_value,step=p_step,func=p_ableitung, verfahren=esv.Euler_verfahren_r)
    imp_euler_res = esv.generelle_einschritt_verfahren(number_of_iterations=iter,start_value=p_start_value,step=p_step,func=p_ableitung, verfahren=esv.verbessertes_Euler_verfahren_r)
    plt.plot(np.linspace(0,goal,iter+1),np.array(euler_res).T[0],label="Euler Beute")
    plt.plot(np.linspace(0,goal,iter+1),np.array(euler_res).T[1],label="Euler Räuber")
    plt.plot(np.linspace(0,goal,iter+1),np.array(imp_euler_res).T[0],label="Verbessertes Euler Verfahren Beute")
    plt.plot(np.linspace(0,goal,iter+1),np.array(imp_euler_res).T[1],label="Verbessertes Euler Verfahren Räuber")
    
    plt.legend()
    plt.grid()
    plt.show()

if __name__ == "__main__":
    normal()
    lotka_vol(func=None,abl=Lotka_temp,start=[4,2],step=0.125,goal = 20)