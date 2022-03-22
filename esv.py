def Euler_verfahren(step,x,val,func,anzahl_aufrufe=0):
    return val + step * func(x,val),anzahl_aufrufe+1

def generelle_einschritt_verfahren(start_value, x_values, number_of_iterations, step, func, verfahren,steile_abl = 1):
    anzahl_aufrufe = 0
    res_list = [start_value]
    for i in range(number_of_iterations):
        result,anzahl_aufrufe = verfahren(step,x_values[i],res_list[i],func,anzahl_aufrufe)
        res_list.append(result*steile_abl)
    return res_list,anzahl_aufrufe

def Heun_verfahren(step,x,val,func,anzahl_aufrufe=0):
    #1/2 * ( func(x,val) + func(x+step, val + step * func(x,val) )
    temp_val = val + step * func(x,val)
    temp_x = x + step
    return val + 1/2 * step *(func(x,val) + func(temp_x,temp_val)),anzahl_aufrufe+3


    #temp = Euler_verfahren(step,x,val,func,anzahl_aufrufe=0)
    #return 1/2 * (func(val) + func(temp) )
    #return val + 1/2 * step * (func(val) + func(temp))

    #return 1/2 * val + 1/2 * (temp + step * func(tem))

def verbessertes_Euler_verfahren(step,x,val,func,anzahl_aufrufe=0):
    temp_x = x + 1/2 * step
    temp_val = val + 1/2 * step * func(x,val)
    return val + step * func(temp_x,temp_val),anzahl_aufrufe+2
    #temp = val + step/2 * func(val) # 1/2 * ( func(x,val) + func(x + step, val + step * func(x,val)) )
    return val + step * func(temp)

def mittelpunkt_verfahren(step,values,func,number_of_iterations,var,anzahl_aufrufe=0):
    return values[0] + step*func(values[0] + (step/2.0)*func(values[0])),anzahl_aufrufe+2