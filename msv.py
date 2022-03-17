import esv
    
def zwei_schritt_Adams_Bashforth_verfahren_r(step,val_1,val_2,func):
    return val_1 + 3/2 * step * func(val_1) - 1/2 * step * func(val_2)

def zwei_schritt_Adams_Bashforth_verfahren(start_value, number_of_iterations, step,func):
    res_list = [start_value,esv.Euler_verfahren_r(step,start_value,func)]
    for i in range(number_of_iterations-1):
        res_list.append(zwei_schritt_Adams_Bashforth_verfahren_r(step,res_list[i+1],res_list[i],func))
    return res_list
    