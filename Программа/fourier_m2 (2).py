import numpy as np
import math
import locale
import sys

if sys.platform == "win32":
    locale.setlocale(locale.LC_ALL, 'rus')
else:
    locale.setlocale(locale.LC_ALL, ('RU','UTF8'))


def ff2(y_data, m, extra_step):
    n = len(y_data)
    t_data = ff2_GetT(n)
    print("ff2_GetT: "+str(len(t_data)))   
    k = ff2_GetK(y_data, t_data)
    print("ff2_GetK: "+str(k))   
    c = ff2_GetC(y_data, t_data, k)
    print("ff2_GetC: "+str(c))   
    a = ff2_GetCoefficientsForA(y_data, t_data, m)
    print("ff2_GetCoefficientsForA")   
    b = ff2_GetCoefficientsForB(y_data, t_data, m)
    print("ff2_GetCoefficientsForB")   
    
    
    results = ff2_GetRes(a, b, t_data, m, k, c, extra_step)
    print("ff2_GetRes: "+str(len(results[0])))   
    
    ff2_WriteResToFile(m, results, k, c, a, b)
    print("ff2_WriteResToFile")   

    return

def ff2_GetT(n):
    return [2*math.pi*i/n for i in range(n)]

def ff2_GetK(y_data, t_data):
    n = len(y_data)
    yt = [y*t for y, t in zip(y_data,t_data)]
    yt_sum = np.sum(yt)
    y_sum = np.sum(y_data)

    tt = [t*t for t in t_data]
    tt_sum = np.sum(tt)
    t_sum = np.sum(t_data)
    
    return (n*yt_sum - t_sum*y_sum)/(n*tt_sum-t_sum*t_sum)

def ff2_GetC(y_data, t_data, k):
    n = len(y_data)
    y_sum = np.sum(y_data)
    t_sum = np.sum(t_data)
    return y_sum/n - t_sum*k/n

def ff2_GetCoefficientsForA(y_data, t_data, m):
    a = []
    for i in range(1, m+1):  
        a.append(ff2_GetA(y_data, t_data, i))
    return a

def ff2_GetA(y_data, t_data, i):
    a_data = [y*math.cos(i*t) for y, t in zip(y_data,t_data)]
    return 2*np.average(a_data)

def ff2_GetCoefficientsForB(y_data, t_data, m):
    b = []
    for k in range(1, m+1):  
        b.append(ff2_GetB(y_data, t_data, k))
    return b

def ff2_GetB(y_data, t_data, i):
    b_data = [y*math.sin(i*t) for y, t in zip(y_data,t_data)]
    return 2*np.average(b_data)

def ff2_GetRes(a, b, t_data, m, k, c, extra_step):
    results = []
    for i in range(m):
        harmonic = ff2_GetHarmonic(a, b, t_data, i, k, c, extra_step)
        results.append(harmonic)
    return results

def ff2_GetHarmonic(a, b, t_data, i, k, c, extra_step):
    harmonic = []
    for t in t_data:
        sin_cos_sum = ff2_GetSinCosSum(a, b, t, i)
        harmonic.append(c+ k*t + sin_cos_sum*(1 + k*t/c))
    #extrapolate for next t value
    t_data_len = len(t_data)
    for j in range(t_data_len, t_data_len + extra_step):
        t_next = 2*math.pi*j/t_data_len
        sin_cos_sum = ff2_GetSinCosSum(a, b, t_next, i)
        harmonic.append(c + k*t_next + sin_cos_sum*(1 + k*t_next/c))
        
    return harmonic   
        
def ff2_GetSinCosSum(a, b, t, i):
    sin_cos_sum = 0
    for k in range(1, i+2):
        sin_cos_sum += a[k-1]*math.cos(k*t) + b[k-1]*math.sin(k*t)
    return sin_cos_sum        

def ff2_WriteResToFile(m, results, k, c, a, b):
    with open('fourier_results_m2.txt', 'w') as out_file:
        out_file.write("k="+locale.format_string("%f", k)+"\n")
        out_file.write("c="+locale.format_string("%f", c)+"\n")
        # a0, a1, ...
        for k in range(len(a)):
            out_file.write("a"+str(k+1)+"="+locale.format_string("%f", a[k])+"\n")
            out_file.write("b"+str(k+1)+"="+locale.format_string("%f", b[k])+"\n")
            
        out_file.write("\n")
        for i, y_i in enumerate(results):
            out_file.write("n="+str(i+1)+"\n")
            # y(t) values    
            for y in y_i:
                s = locale.format_string("%f", y)
                out_file.write(s+"\n")
            out_file.write("\n")
        
def ff2_readFile(file_name):
    y_data = []
    with open(file_name) as in_file:
        lines = in_file.readlines()
        for line in lines:
            try:
                if line != '':
                    y_data.append(locale.atof(line.strip()))
            except Exception as inst:
                print(inst)
    return y_data
        
try:
    y_data=ff2_readFile('fourier_m2_data.txt')
    m=int(input("Enter m:"))
    extra_step=int(input("Enter step count:"))
    ff2(y_data, m, extra_step)
except Exception as er:
    print(er)   