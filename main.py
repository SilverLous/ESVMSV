import numpy as np
import matplotlib.pyplot as plt
import esv
import msv
from numpy import sin,cos,tan
from scipy.integrate import odeint
import pandas as pd


step = 0.125
goal = 1.5


def euler_function(val,t=None):
    return np.exp(val) 

def euler_f(val,t=None):
    return val

def log_function(val):
    return np.log(val+1)

def log_ableitung(val,t=None):
    return 1/(val+1)

def sin_func(val):
    return sin(val)

def cos_func(val):
    return cos(val)

def tang_func(val):
    return tan(val)

def tang_func_r(val,t=None):
    return 1 + pow(val,2)

def Lotka_Volterra_derivative(X, alpha, beta, delta, gamma):
    x, y = X[0],X[1]
    dotx = x * (alpha - beta * y)   # growth - predation
    doty = y * (-delta + gamma * x) # growth - mortality (or other)
    return np.array([dotx, doty])

def Lotka_temp(val):
    return Lotka_Volterra_derivative(val,1,1,1,1)

ziel_func = tang_func

ableitung = tang_func_r

START_VALUE = ziel_func(0) # ziel_func(0)

def normal(p_ziel_func,p_GDL,p_step,p_goal_number):
    n_ziel_func = p_ziel_func
    Gdl = p_GDL
    step = p_step
    goal = p_goal_number
    start = n_ziel_func(0)
    einschritt_verfahren = [esv.Euler_verfahren_r,esv.verbessertes_Euler_verfahren_r,esv.Heun_verfahren_r]
    mehrschritt_verfahren = [msv.zwei_schritt_Adams_Bashforth_verfahren_r]
    num_verfahren = len(einschritt_verfahren)+len(mehrschritt_verfahren)
    sqrt_of_l = round((num_verfahren+1) ** 0.5 + 0.5)
    fig = plt.figure(figsize=(12, 12))
    index = 1
    iter = int(goal/step) # Number of iterations
    x = np.linspace(0,goal,200)
    x_array = np.linspace(0,goal,iter+1,endpoint=True)

    for verfahren in einschritt_verfahren:
        verfahren_name = str(verfahren)
        verfahren_name = verfahren_name[10:verfahren_name.find("at 0")]
        fig.add_subplot(sqrt_of_l, sqrt_of_l, index)
        plt.plot(x,n_ziel_func(x),label="Zielfunktion")
        plt.plot(x_array,esv.generelle_einschritt_verfahren(start,iter,step,Gdl,verfahren),label=verfahren_name)
        plt.legend()
        plt.grid()
        index+=1
        
    for verfahren in mehrschritt_verfahren:
        verfahren_name = str(verfahren)
        verfahren_name = verfahren_name[10:verfahren_name.find("at 0")]
        fig.add_subplot(sqrt_of_l, sqrt_of_l, index)
        plt.plot(x,n_ziel_func(x),label="Zielfunktion")
        plt.plot(x_array,msv.generelle_mehrschritt_verfahren(start,iter,step,Gdl,verfahren),label=verfahren_name)
        plt.legend()
        plt.grid()
        index+=1
    fig.add_subplot(sqrt_of_l, sqrt_of_l, index)
    plt.plot(x,n_ziel_func(x),label="Zielfunktion")
    plt.plot(x_array,odeint(Gdl,start,x_array),label="Isode")
    plt.legend()
    plt.grid()
    plt.show()

def lotka_vol(p_func,p_abl,p_start,p_step,p_goal):
    goal = p_goal
    step = p_step
    start = p_start
    func = p_func
    abl = p_abl
    iter = int(goal/step) # Number of iterations
    x = np.linspace(0,goal,200)
    #plt.plot(x,ziel_func(x),label="Zielfunktion")

    euler_res = esv.generelle_einschritt_verfahren(number_of_iterations=iter,start_value=start,step=step,func=abl, verfahren=esv.Euler_verfahren_r)
    imp_euler_res = esv.generelle_einschritt_verfahren(number_of_iterations=iter,start_value=start,step=step,func=abl, verfahren=esv.verbessertes_Euler_verfahren_r)
    plt.plot(np.linspace(0,goal,iter+1),np.array(euler_res).T[0],label="Euler Beute")
    plt.plot(np.linspace(0,goal,iter+1),np.array(euler_res).T[1],label="Euler Räuber")
    plt.plot(np.linspace(0,goal,iter+1),np.array(imp_euler_res).T[0],label="Verbessertes Euler Verfahren Beute")
    plt.plot(np.linspace(0,goal,iter+1),np.array(imp_euler_res).T[1],label="Verbessertes Euler Verfahren Räuber")
    
    plt.legend()
    plt.grid()
    plt.show()

if __name__ == "__main__":
    normal(euler_function,euler_f    ,0.125,15)
    normal(tang_func,tang_func_r     ,0.125,1.5)
    normal(log_function,log_ableitung,0.125,15)
    #lotka_vol(func=None,abl=Lotka_temp,start=[4,2],step=0.125,goal = 20)