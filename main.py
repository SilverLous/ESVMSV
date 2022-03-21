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

def euler_function(x, t=None):
    return np.exp(x)

def euler_f_ableitung(x, val, t=None):
    return val

def n_euler_function(x,t=None):
    return -1 * np.exp(-x)

def n_euler_f_ableitung(x, val, t=None):
    return -val

def log_function(x):
    return np.log(x+1)

def log_ableitung(x,val,t=None):
    return 1/(x+1)

def sin_func(x,val=None):
    return sin(x)

def cos_func(x,val=None,t=None):
    return cos(x)

def banalt(x,val=None):
    return x**2

def banal_abl(x,val,t=None):
    return 2*x

def tang_func(x,val=None):
    return tan(x)

def tang_func_ableitung(x, val, t=None):
    return 1 + pow(val,2)

def Lotka_Volterra_derivative(X, alpha, beta, delta, gamma):
    x, y = X[0],X[1]
    dotx = x * (alpha - beta * y)   # growth - predation
    doty = y * (-delta + gamma * x) # growth - mortality (or other)
    return np.array([dotx, doty])

def Lotka_temp(x,val,t=None):
    return Lotka_Volterra_derivative(val,1,1,1,1)

def Lsoda_Lotka_temp(val,t=None):
    return Lotka_Volterra_derivative(val,1,1,1,1)

def plot_details(title,y_ticks=None,zoom=None):
    plt.title(title,fontsize=9)
    plt.legend(fontsize=8)
    plt.xticks(fontsize=5)
    y_ticks=None
    if y_ticks is not None:
        plt.yticks(y_ticks,fontsize=5)
    if zoom is not None:
        plt.axis(zoom)
    plt.grid()

def plot_ziel_func(x_arr,zielfunc,custom_color="lightskyblue"):
    plt.plot(x_arr,zielfunc,label="Zielfunktion",lw=3,c=custom_color)
    


def normal(p_ziel_func,p_GDL,p_step,p_goal_number,var,overwrite_start=None):
    cm = plt.get_cmap('gist_rainbow')
    alle_verfahren = {}
    n_ziel_func = p_ziel_func
    Gdl = p_GDL
    step = p_step
    goal = p_goal_number
    start = n_ziel_func(0)
    ziel_func_name =  str(n_ziel_func)
    ziel_func_name =  ziel_func_name[10:ziel_func_name.find("at 0")-1]

    selbst_aufrufende_func = ["euler_function","n_euler_function","tang_func"]

    scipy_ver_erlaubt = ziel_func_name in selbst_aufrufende_func
    normale_ver_erlaubt = True

    if overwrite_start is not None:
        start = overwrite_start
    einschritt_verfahren =  []
    mehrschritt_verfahren = []
    scipy_verfahren = []
    if normale_ver_erlaubt:
        einschritt_verfahren = [esv.Euler_verfahren,esv.verbessertes_Euler_verfahren,esv.Heun_verfahren]
        mehrschritt_verfahren = [msv.generell_Adams_Bashforth_Verfahren,msv.mittelpunkt_verfahren]
    num_diff_colors = len(einschritt_verfahren)+len(mehrschritt_verfahren)+len(scipy_verfahren)+5
    if scipy_ver_erlaubt:
        scipy_verfahren = [inte.Radau,inte.RK23,inte.RK45,inte.DOP853,inte.BDF,inte.LSODA]
    num_verfahren = len(einschritt_verfahren)+len(mehrschritt_verfahren)+len(scipy_verfahren)+5
    sqrt_of_l = round((num_verfahren) ** 0.5 + 0.5)
    index = 1
    iter = int(goal/step) # Number of iterations
    x = np.linspace(0,goal,200)
    x_array = np.linspace(0,goal,iter+1,endpoint=True)
    ziel_func = n_ziel_func(x)
    ziel_func_arr = n_ziel_func(x_array)
    y_ticks = np.linspace(start-1,n_ziel_func(goal),5)
    fig = plt.figure(figsize=(16,10))
    fig.suptitle(ziel_func_name)
    c_zoom = [0,goal,min(ziel_func)-1,max(ziel_func)+1]

    if normale_ver_erlaubt:

        for verfahren in einschritt_verfahren:
            verfahren_name = str(verfahren)
            verfahren_name = verfahren_name[10:verfahren_name.find("at 0")-1]
            fig.add_subplot(sqrt_of_l, sqrt_of_l, index)
            esv_res_list,aufrufe = esv.generelle_einschritt_verfahren(start, x_array, iter, step, Gdl, verfahren)
            print(f"Das {verfahren_name} hat {aufrufe} Funktionsaufrufe über {iter} Schritten also eine Rate von {aufrufe/iter} Aufrufen pro Schritt")
            plot_ziel_func(x,ziel_func)
            plt.plot(x_array,esv_res_list,"--",label=verfahren_name,c="r")
            diff_list = []
            for i,value in enumerate(esv_res_list):
                diff_list.append(ziel_func_arr[i]-value)
            line = plt.plot(x_array,diff_list,label="Differenz")
            line[0].set_color(cm(index/3*3/num_diff_colors))
            plot_details(verfahren_name,y_ticks,c_zoom)
            alle_verfahren[verfahren_name] = diff_list
            index+=1

            
        for verfahren in mehrschritt_verfahren:
            verfahren_name = str(verfahren)
            verfahren_name = verfahren_name[10:verfahren_name.find("at 0")-1]
            if verfahren_name=="generell_Adams_Bashforth_Verfahren":
                for stufen in range(2,6):
                    verfahren_name = f"{stufen} Schritt AB Verfahren"
                    fig.add_subplot(sqrt_of_l, sqrt_of_l, index)
                    msv_res_list,aufrufe = msv.generelle_mehrschritt_verfahren(start, x_array, iter, step, Gdl, verfahren, stufen)
                    print(f"Das {verfahren_name} hat {aufrufe} Funktionsaufrufe über {iter} Schritten also eine Rate von {aufrufe/iter} Aufrufen pro Schritt")
                    plot_ziel_func(x,ziel_func)
                    plt.plot(x_array,msv_res_list,"--",label=verfahren_name,c="r")
                    diff_list = []
                    for i,value in enumerate(msv_res_list):
                        diff_list.append(ziel_func_arr[i]-value)
                    line = plt.plot(x_array,diff_list,label="Differenz")
                    line[0].set_color(cm(index/3*3/num_diff_colors))
                    plot_details(verfahren_name,y_ticks,c_zoom)
                    index+=1
                    alle_verfahren[verfahren_name] = diff_list
            else:
                fig.add_subplot(sqrt_of_l, sqrt_of_l, index)
                msv_res_list,aufrufe = msv.generelle_mehrschritt_verfahren(start, x_array, iter, step, Gdl, verfahren, 2)
                
                print(f"Das {verfahren_name} hat {aufrufe} Funktionsaufrufe über {iter} Schritten also eine Rate von {aufrufe/iter} Aufrufen pro Schritt")
                plot_ziel_func(x,ziel_func)
                plt.plot(x_array,msv_res_list,"--",label=verfahren_name,c="r")
                diff_list = []
                for i,value in enumerate(msv_res_list):
                    diff_list.append(ziel_func_arr[i]-value)
                line = plt.plot(x_array,diff_list,label="Differenz")
                line[0].set_color(cm(index/3*3/num_diff_colors))
                plot_details(verfahren_name,y_ticks,c_zoom)
                index+=1
                alle_verfahren[verfahren_name] = diff_list

        fig.add_subplot(sqrt_of_l, sqrt_of_l, index)
        diff_list = []
        isode_res_list = odeint(reversed_args(Gdl),start,x_array)
        
        plot_ziel_func(x,ziel_func)
        plt.plot(x_array,isode_res_list,"--",label="lsoda",c="r") #"lsoda Räuber"
        for i,value in enumerate(isode_res_list):
            diff_list.append(ziel_func_arr[i]-value[0])
        alle_verfahren["lsoda"] = diff_list
        line = plt.plot(x_array,diff_list,label="Differenz")
        line[0].set_color(cm(index/3*3/num_diff_colors))
        plot_details("lsoda",y_ticks,c_zoom)
        index+=1
    
    if scipy_ver_erlaubt:
        for verfahren in scipy_verfahren:
            verfahren_name = str(verfahren)
            verfahren_name = verfahren_name[29:verfahren_name.find("at 0")]
            fig.add_subplot(sqrt_of_l, sqrt_of_l, index)
            diff_list = []
            plot_ziel_func(x,ziel_func)
            start_arr = [start]
            res,_ = esv.Euler_verfahren(step,x_array[0],start,Gdl)
            start_arr.append(res)
            scipy_res = verfahren(fun=Gdl,t0=0,y0=start_arr,t_bound=goal,max_step=step)
            for i in range(len(x_array)-1):
                scipy_res.step()
            scipy_res_list = scipy_res.dense_output().__call__(x_array)[0]
            #plt.plot(x_array,BDF_res_list,label=verfahren_name)
            scipy_res_list2 = scipy_res.dense_output().__call__(x_array)[1]
            #plt.plot(x_array,BDF_res_list2,label=verfahren_name)
            plt.fill_between(x_array,scipy_res_list,scipy_res_list2,color="r",alpha=0.5,label=verfahren_name)
            #for i,value in enumerate(scipy_res_list):
            #    diff_list.append(n_ziel_func(x_array[i])-value)
            #alle_verfahren[verfahren_name] = diff_list
            plot_details(verfahren_name,y_ticks,c_zoom)
            index+=1
    

    fig.add_subplot(sqrt_of_l, sqrt_of_l, index)
    for i,key in enumerate(alle_verfahren.keys()):
        line = plt.plot(x_array,alle_verfahren[key],label = i)
        line[0].set_color(cm((i+1)/3*3/num_diff_colors))
    #plt.subplots_adjust(hspace=1)
    plt.legend()
    plot_details("Differenzen zur Zielfunktion")
    # plt.tight_layout()
    if overwrite_start is not None:
        plt.savefig(f"output_{ziel_func_name}_mit_fehler.png", bbox_inches="tight")
    else:
        plt.savefig(f"output_{ziel_func_name}.png", bbox_inches="tight")
    plt.show()
    return x_array,alle_verfahren

def lotka_vol(p_func,p_abl,p_start,p_step,p_goal):
    goal = p_goal
    step = p_step
    start = p_start
    func = p_func
    abl = p_abl
    iter = int(goal/step) # Number of iterations
    x = np.linspace(0,goal,200)
    x_array = np.linspace(0,goal,iter+1,endpoint=True)
    #plt.plot(x,ziel_func(x),label="Zielfunktion")
    

    cm = plt.get_cmap('gist_rainbow')
    alle_verfahren = {}
    n_ziel_func = p_func
    Gdl = p_abl
    step = p_step
    goal = p_goal
    start = p_start
    einschritt_verfahren =  []
    mehrschritt_verfahren = []
    einschritt_verfahren = [esv.Euler_verfahren,esv.verbessertes_Euler_verfahren,esv.Heun_verfahren]
    mehrschritt_verfahren = [msv.generell_Adams_Bashforth_Verfahren,msv.mittelpunkt_verfahren]
    num_diff_colors = len(einschritt_verfahren)+len(mehrschritt_verfahren)+3
    num_verfahren = len(einschritt_verfahren)+len(mehrschritt_verfahren)+3
    sqrt_of_l = round((num_verfahren) ** 0.5 + 0.5)
    index = 1
    iter = int(goal/step) # Number of iterations
    x = np.linspace(0,goal,200)
    x_array = np.linspace(0,goal,iter+1,endpoint=True)
    ziel_func = odeint(Lsoda_Lotka_temp,start,x_array)
    ziel_func = ziel_func.T
    ziel_func_arr = ziel_func[0]
    ziel_func_name =  str(n_ziel_func)
    ziel_func_name =  "Lotka_Volterra"
    y_ticks = np.linspace(start[0]-1,max(ziel_func_arr),5)
    
    c_zoom = [0,goal,ziel_func_arr[goal]*-0.1,ziel_func_arr[goal]]
    fig = plt.figure(figsize=(16,10))
    fig.suptitle(ziel_func_name)

    for verfahren in einschritt_verfahren:
        verfahren_name = str(verfahren)
        verfahren_name = verfahren_name[10:verfahren_name.find("at 0")-1]
        fig.add_subplot(sqrt_of_l, sqrt_of_l, index)
        esv_res_list,aufrufe = esv.generelle_einschritt_verfahren(start, x_array, iter, step, Gdl, verfahren)
        print(f"Das {verfahren_name} hat {aufrufe} Funktionsaufrufe über {iter} Schritten also eine Rate von {aufrufe/iter} Aufrufen pro Schritt")
        plot_ziel_func(x_array,ziel_func_arr)
        plot_ziel_func(x_array,ziel_func[1],"b")
        plt.plot(x_array,np.array(esv_res_list).T[0],"--",label=verfahren_name+" Räuber",c="r")
        plt.plot(x_array,np.array(esv_res_list).T[1],"--",label=verfahren_name+" Beute" ,c="g")            
        plot_details(verfahren_name,y_ticks)
        index+=1
        #plt.show()
        
    for verfahren in mehrschritt_verfahren:
        verfahren_name = str(verfahren)
        verfahren_name = verfahren_name[10:verfahren_name.find("at 0")-1]
        if verfahren_name=="generell_Adams_Bashforth_Verfahren":
            for stufen in range(2,6):
                fig.add_subplot(sqrt_of_l, sqrt_of_l, index)
                msv_res_list,aufrufe = msv.generelle_mehrschritt_verfahren(start,x_array, iter, step, Gdl, verfahren, stufen)
                print(f"Das {verfahren_name} hat {aufrufe} Funktionsaufrufe über {iter} Schritten also eine Rate von {aufrufe/iter} Aufrufen pro Schritt")
                plot_ziel_func(x_array,ziel_func_arr)
                plot_ziel_func(x_array,ziel_func[1],"b")
                plt.plot(x_array,np.array(msv_res_list).T[0],"--",label=f"{stufen} Schritt AB Verfahren Räuber",c="r")
                plt.plot(x_array,np.array(msv_res_list).T[1],"--",label=f"{stufen} Schritt AB Verfahren Beute",c="g")                   
                plot_details(f"{stufen} Schritt AB Verfahren",y_ticks)
                index+=1
                #plt.show()
        else:
            fig.add_subplot(sqrt_of_l, sqrt_of_l, index)
            msv_res_list,aufrufe = msv.generelle_mehrschritt_verfahren(start, x_array, iter, step, Gdl, verfahren, 2)
            print(f"Das {verfahren_name} hat {aufrufe} Funktionsaufrufe über {iter} Schritten also eine Rate von {aufrufe/iter} Aufrufen pro Schritt")
            plot_ziel_func(x_array,ziel_func_arr)                    
            plot_ziel_func(x_array,ziel_func[1],"b")
            plt.plot(x_array,np.array(msv_res_list).T[0],"--",label=verfahren_name+" Räuber",c="r")
            plt.plot(x_array,np.array(msv_res_list).T[1],"--",label=verfahren_name+" Beute" ,c="g")               
            plot_details(verfahren_name,y_ticks)
            index+=1
    plt.savefig(f"output_{ziel_func_name}.png", bbox_inches="tight")
    plt.show()

    """
    euler_res = esv.generelle_einschritt_verfahren(number_of_iterations=iter,start_value=start,step=step,func=abl, verfahren=esv.Euler_verfahren)
    imp_euler_res = esv.generelle_einschritt_verfahren(number_of_iterations=iter,start_value=start,step=step,func=abl, verfahren=esv.verbessertes_Euler_verfahren)
    plt.plot(np.linspace(0,goal,iter+1),np.array(euler_res).T[0],label="Euler Beute")
    plt.plot(np.linspace(0,goal,iter+1),np.array(euler_res).T[1],label="Euler Räuber")
    plt.plot(np.linspace(0,goal,iter+1),np.array(imp_euler_res).T[0],label="Verb. Euler Verfahren Beute")
    plt.plot(np.linspace(0,goal,iter+1),np.array(imp_euler_res).T[1],label="Verb. Euler Verfahren Räuber")
    isode_res_list = odeint(abl,start,x_array)
    plt.plot(np.linspace(0,goal,iter+1),np.array(isode_res_list).T[0],label="Lsoda Verfahren Beute")
    plt.plot(np.linspace(0,goal,iter+1),np.array(isode_res_list).T[1],label="Lsoda Verfahren Räuber")
    
    plt.legend(fontsize=8,loc="upper right")
    plt.grid()"""




if __name__ == "__main__":
    #x_array,euler_with_errors = normal(euler_function, euler_f_ableitung, 0.125, 2, 5, 0.8)
    #_,euler_without_errors = normal(euler_function, euler_f_ableitung, 0.125, 2, 5)
    #for i,key in enumerate(euler_with_errors.keys()):
    #    diff_list = []
    #    for value1,value2 in zip(euler_with_errors[key],euler_without_errors[key]):
    #        diff_list.append(value1-value2)
    #    plt.plot(x_array,diff_list,label = key)
    #plt.legend()
    #plot_details("Differenzen zur Zielfunktion")
    #plt.savefig("output_differenzen_bei_fehler")
    #plt.show()
    #normal(n_euler_function,n_euler_f_ableitung    ,0.125,2,  5)
    #normal(tang_func,tang_func_ableitung     ,0.125,1.5,5)
    #normal(log_function,log_ableitung,0.125,15, 5)
    #normal(sin_func,cos_func,0.125,15, 5)
    #normal(banalt,banal_abl,0.125,15, 5)
    lotka_vol(p_func=None,p_abl=Lotka_temp,p_start=[4,2],p_step=0.125,p_goal = 20)