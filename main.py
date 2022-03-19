from functools import wraps
import numpy as np
import matplotlib.pyplot as plt
import esv
import msv
from numpy import sin,cos,tan
from scipy.integrate import odeint
import scipy.integrate as inte
import pandas as pd


step = 0.125
goal = 1.5

def reversed_args(f):
    @wraps(f)
    def g(*args):
        return f(*args[::-1])
    return g

def euler_function(val,t=None):
    return np.exp(val) 

def euler_f(val,t=None):
    return val

def euler_f_temp(t,val):
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

def normal(p_ziel_func,p_GDL,p_step,p_goal_number):
    alle_verfahren = {}
    n_ziel_func = p_ziel_func
    Gdl = p_GDL
    step = p_step
    goal = p_goal_number
    start = n_ziel_func(0)
    einschritt_verfahren = [esv.Euler_verfahren,esv.verbessertes_Euler_verfahren,esv.Heun_verfahren]
    mehrschritt_verfahren = [msv.zwei_schritt_Adams_Bashforth_verfahren,msv.backward_differentiation_verfahren,msv.generell_Adams_Bashforth_Verfahren]
    scipy_verfahren = [inte.BDF,inte.LSODA,inte.Radau,inte.RK23,inte.RK45,inte.DOP853]
    num_verfahren = len(einschritt_verfahren)+len(mehrschritt_verfahren)+len(scipy_verfahren)
    sqrt_of_l = round((num_verfahren+2) ** 0.5 + 0.5)
    fig = plt.figure(figsize=(12, 12))
    index = 1
    iter = int(goal/step) # Number of iterations
    x = np.linspace(0,goal,200)
    x_array = np.linspace(0,goal,iter+1,endpoint=True)
    ziel_func_name =  str(n_ziel_func)
    ziel_func_name =  ziel_func_name[10:ziel_func_name.find("at 0")]

    for verfahren in einschritt_verfahren:
        verfahren_name = str(verfahren)
        verfahren_name = verfahren_name[10:verfahren_name.find("at 0")]
        fig.add_subplot(sqrt_of_l, sqrt_of_l, index)
        plt.plot(x,n_ziel_func(x),label="Zielfunktion")
        esv_res_list = esv.generelle_einschritt_verfahren(start,iter,step,Gdl,verfahren)
        plt.plot(x_array,esv_res_list,label=verfahren_name)
        diff_list = []
        for i,value in enumerate(esv_res_list):
            diff_list.append(n_ziel_func(x_array[i])-value)
        plt.plot(x_array,diff_list,label="Differenz")
        plt.title(verfahren_name)
        fig.suptitle(ziel_func_name)
        plt.legend()
        plt.grid()
        alle_verfahren[verfahren_name] = diff_list
        index+=1
        
    for verfahren in mehrschritt_verfahren:
        verfahren_name = str(verfahren)
        verfahren_name = verfahren_name[10:verfahren_name.find("at 0")]
        ziel_func_name =  str(n_ziel_func)
        ziel_func_name =  ziel_func_name[10:ziel_func_name.find("at 0")]
        fig.add_subplot(sqrt_of_l, sqrt_of_l, index)
        plt.plot(x,n_ziel_func(x),label="Zielfunktion")
        msv_res_list = msv.generelle_mehrschritt_verfahren(start,iter,step,Gdl,verfahren)
        plt.plot(x_array,msv_res_list,label=verfahren_name)
        diff_list = []
        for i,value in enumerate(esv_res_list):
            diff_list.append(n_ziel_func(x_array[i])-value)
        plt.plot(x_array,diff_list,label="Differenz")
        plt.title(verfahren_name)
        plt.legend()
        plt.grid()
        index+=1
        alle_verfahren[verfahren_name] = diff_list

    fig.add_subplot(sqrt_of_l, sqrt_of_l, index)
    diff_list = []
    plt.plot(x,n_ziel_func(x),label="Zielfunktion")
    isode_res_list = odeint(Gdl,start,x_array)
    plt.plot(x_array,isode_res_list,label="lsoda")
    for i,value in enumerate(isode_res_list):
        diff_list.append(n_ziel_func(x_array[i])-value)
    alle_verfahren["lsoda"] = diff_list
    plt.plot(x_array,diff_list,label="Differenz")
    plt.legend()
    plt.grid()
    index+=1

    for verfahren in scipy_verfahren:
        verfahren_name = str(verfahren)
        verfahren_name = verfahren_name[29:verfahren_name.find("at 0")]
        fig.add_subplot(sqrt_of_l, sqrt_of_l, index)
        diff_list = []
        plt.plot(x,n_ziel_func(x),label="Zielfunktion")
        start_arr = [start]
        start_arr.append(esv.Euler_verfahren(step,start,Gdl))
        scipy_res = verfahren(reversed_args(Gdl),0,start_arr,goal,max_step=step)
        for i in range(len(x_array)-1):
            scipy_res.step()
        scipy_res_list = scipy_res.dense_output().__call__(x_array)[0]
        #plt.plot(x_array,BDF_res_list,label=verfahren_name)
        scipy_res_list2 = scipy_res.dense_output().__call__(x_array)[1]
        #plt.plot(x_array,BDF_res_list2,label=verfahren_name)
        plt.fill_between(x_array,scipy_res_list,scipy_res_list2,color="r",alpha=0.5,label=verfahren_name)
        for i,value in enumerate(scipy_res_list):
            diff_list.append(n_ziel_func(x_array[i])-value)
        alle_verfahren[verfahren_name] = diff_list
        plt.legend()
        plt.grid()
        index+=1


    fig.add_subplot(sqrt_of_l, sqrt_of_l, index)
    for key in alle_verfahren.keys():
        plt.plot(x_array,alle_verfahren[key],label=key[0])
    plt.legend(fontsize=8,)
    plt.grid()
    plt.title("Differenzen zur Zielfunktion")
    plt.savefig("output.png", bbox_inches="tight")
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

    euler_res = esv.generelle_einschritt_verfahren(number_of_iterations=iter,start_value=start,step=step,func=abl, verfahren=esv.Euler_verfahren)
    imp_euler_res = esv.generelle_einschritt_verfahren(number_of_iterations=iter,start_value=start,step=step,func=abl, verfahren=esv.verbessertes_Euler_verfahren)
    plt.plot(np.linspace(0,goal,iter+1),np.array(euler_res).T[0],label="Euler Beute")
    plt.plot(np.linspace(0,goal,iter+1),np.array(euler_res).T[1],label="Euler Räuber")
    plt.plot(np.linspace(0,goal,iter+1),np.array(imp_euler_res).T[0],label="Verbessertes Euler Verfahren Beute")
    plt.plot(np.linspace(0,goal,iter+1),np.array(imp_euler_res).T[1],label="Verbessertes Euler Verfahren Räuber")
    
    plt.legend()
    plt.grid()
    plt.show()

if __name__ == "__main__":
    normal(euler_function,euler_f    ,0.125,2)
    normal(tang_func,tang_func_r     ,0.125,1.5)
    normal(log_function,log_ableitung,0.125,15)
    #lotka_vol(func=None,abl=Lotka_temp,start=[4,2],step=0.125,goal = 20)