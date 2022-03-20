def Euler_verfahren(step,x,val,func):
    return val + step * func(x,val)

def generelle_einschritt_verfahren(start_value, x_values, number_of_iterations, step, func, verfahren):
    res_list = [start_value]
    for i in range(number_of_iterations):
        res_list.append(verfahren(step,x_values[i],res_list[i],func))
    return res_list

def Heun_verfahren(step,x,val,func):
    #1/2 * ( func(x,val) + func(x+step, val + step * func(x,val) )
    temp_val = val + step * func(x,val)
    temp_x = x + step
    return val + 1/2 * step *(func(x,val) + func(temp_x,temp_val))


    #temp = Euler_verfahren(step,x,val,func)
    #return 1/2 * (func(val) + func(temp) )
    #return val + 1/2 * step * (func(val) + func(temp))

    #return 1/2 * val + 1/2 * (temp + step * func(tem))

def verbessertes_Euler_verfahren(step,x,val,func):
    temp_x = x + 1/2 * step
    temp_val = val + 1/2 * step * func(x,val)
    return val + step * func(temp_x,temp_val)
    #temp = val + step/2 * func(val) # 1/2 * ( func(x,val) + func(x + step, val + step * func(x,val)) )
    return val + step * func(temp)

def mittelpunkt_verfahren(step,values,func,number_of_iterations,var):
    return values[0] + step*func(values[0] + (step/2.0)*func(values[0]))