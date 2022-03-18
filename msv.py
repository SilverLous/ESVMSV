import esv
from scipy.integrate import odeint
import numpy as np
import main
def zwei_schritt_Adams_Bashforth_verfahren_r(step,values,func,number_of_iterations):
    return values[0] + 3/2 * step * func(values[0]) - 1/2 * step * func(values[1])

def generelle_mehrschritt_verfahren(start_value, number_of_iterations, step,func,verfahren):
    res_list = [start_value,esv.Euler_verfahren_r(step,start_value,func)]
    for i in range(number_of_iterations-1):
        parameter_array = [res_list[i+1],res_list[i]]
        res_list.append(verfahren(step,parameter_array,func,number_of_iterations))
    return res_list
