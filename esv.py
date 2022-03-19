def Euler_verfahren(step,val,func):
    return val + step * func(val)

def generelle_einschritt_verfahren(start_value,number_of_iterations,step,func,verfahren):
    res_list = [start_value]
    for i in range(number_of_iterations):
        res_list.append(verfahren(step,res_list[i],func))
    return res_list

def Heun_verfahren(step,val,func):
    temp = Euler_verfahren(step,val,func)
    #return 1/2 * (func(val) + func(temp) )
    return val + 1/2 * step * (func(val) + func(temp))

    return 1/2 * val + 1/2 * (temp + step * func(tem))

def verbessertes_Euler_verfahren(step,val,func):
    temp = val + step/2 * func(val)
    return val + step * func(temp)

def mittelpunkt_verfahren(step,values,func,number_of_iterations,var):
    return values[0] + step*func(values[0] + (step/2.0)*func(values[0]))