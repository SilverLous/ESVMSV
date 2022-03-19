import esv
from scipy.integrate import odeint
import numpy as np
import main

k_Adams_Bashforth = {
    '2': [1/2, 3, -1]
}
def zwei_schritt_Adams_Bashforth_verfahren(step,values,func,number_of_iterations):
    return values[0] + 3/2 * step * func(values[0]) - 1/2 * step * func(values[1])

def generelle_mehrschritt_verfahren(start_value, number_of_iterations, step,func,verfahren):
    res_list = [start_value,esv.Euler_verfahren(step,start_value,func)]
    for i in range(number_of_iterations-1):
        parameter_array = [res_list[i+1],res_list[i]]
        res_list.append(verfahren(step,parameter_array,func,number_of_iterations))
    return res_list

def backward_differentiation_verfahren(step,values,func,number_of_iterations):
     if len(values)>1:
         return 4/3 * values[0] - 1/3 * values[1] + 2/3*step*func(values[1])

def generell_Adams_Bashforth_Verfahren(step,values,func,number_of_iterations):
    iter = str(number_of_iterations)
    koeffizienten = k_Adams_Bashforth[iter]
    res = 0
    for i in range(len(koeffizienten)-1):
        res += koeffizienten[i+1] * func(values[i+1])
    return values[0] + koeffizienten[0] * step * res
