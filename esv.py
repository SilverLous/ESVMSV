def euler_verfahren_r(step,val,func):
    print(step,val,func(val))
    return val + step * func(val)

def euler_verfahren(start_value, number_of_iterations, step,func ):
    res_list = [start_value]
    for i in range(number_of_iterations):
        res_list.append(euler_verfahren_r(step,res_list[i],func))
    return res_list