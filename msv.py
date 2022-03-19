import esv
from scipy.integrate import odeint
import numpy as np
import main

k_Adams_Bashforth = {
    '1': [1, 1],
    '2': [1/2, 3, -1],
    '3': [1/12, 23, -16, 5],
    '4': [1/24, 55, -59, 37,-9],
    '5': [1/720, 1901, -2774, 2616, -1274, 251]
}
def zwei_schritt_Adams_Bashforth_verfahren(step,values,func,number_of_iterations,var):
    return values[0] + 3/2 * step * func(values[0]) - 1/2 * step * func(values[1])

def generelle_mehrschritt_verfahren(start_value,number_of_iterations,step,func,verfahren,var):
    start_value = [start_value]
    #res_list = [start_value,esv.Euler_verfahren(step,start_value,func)]
    res_list = generate_starting_values(start_value,var,step,func)
    #print(res_list)
    #for i in range(number_of_iterations-1):
    #print(number_of_iterations)
    for i in range(number_of_iterations - (var-1)):
        #parameter_array = [res_list[i+1],res_list[i]]
        parameter_array = []
        for j in range(var):#res_list[var-1-i]
            parameter_array.append(res_list[i+(var-j-1)])
        #print(parameter_array)
        res_list.append(verfahren(step,parameter_array,func,number_of_iterations,var))
        #print(res_list)
    return res_list

# def backward_differentiation_verfahren(step,values,func,number_of_iterations,var):
#      if len(values)>1:
#          return 4/3 * values[0] - 1/3 * values[1] + 2/3*step*func(values[1])

def generell_Adams_Bashforth_Verfahren(step,values,func,number_of_iterations,var):
    iter = str(var)
    koeffizienten = k_Adams_Bashforth[iter]
    res = 0
    for i in range(len(koeffizienten)-1):
        # print(i)
        res += koeffizienten[i+1] * func(values[i])

    return values[0] + koeffizienten[0] * step * res

def generate_starting_values(curr_list,iter,step,func):
    curr_list.append(esv.Euler_verfahren(step,curr_list[len(curr_list)-1],func))
    if len(curr_list) == iter:
        return curr_list
    return generate_starting_values(curr_list,iter,step,func)




def mittelpunkt_verfahren(step,values,func,number_of_iterations,var):
    return values[1] + 2 * step*func(values[0])

