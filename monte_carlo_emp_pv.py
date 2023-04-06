import numpy as np
import os



def __monte_carlo__(T,rvs,n,k,other=None):
    t=np.zeros(k)
    if(other!=None):
        for i in range(k):
            x=rvs(n)
            t[i]=T(x, other)
    else:
        for i in range(k):
            x=rvs(n)
            t[i]=T(x)
    return t


monte_carlo_folder_name="monte carlo"

def monte_carlo_emp_pv(T,rvs,n,k,other=None,name=None):
    """
    T- функция вычисляющая статистику критерия
    rvs(n) - функция генерации выборки
    n - размер выборок
    k - число повторов в методе монте карло
    other - дополнительные параметры для вычисления статистики ( T(x,other) )
    O(n*k*log(k) )
    """
    
    file_name=""
    if(monte_carlo_folder_name!=None):
        file_name=monte_carlo_folder_name+"/"
    file_name+=name+" "+str(n)
    
    if not os.path.isdir(monte_carlo_folder_name):
        os.mkdir(monte_carlo_folder_name)
    
    t=np.array()
    if(name!=None and os.path.exists(file_name) ):
        file=open(file_name,"r")
        t=[float(i) for i in file.read().split()]
        t=np.array(t)
        file.close()
        del file
    
    if(k>len(t)):
        if(name!=None):
            file=open(file_name,"w")
        t0=t
        t1=__monte_carlo__(T,omp,rvs,n,k-len(x),other)
        t=np.zeros(k)
        for i in range(t0):
            t[i]=t0[i]
        for i in range(len(t1)):
            t[i+len(t0)]=t1[i]
        t=np.sort(t)
        if(name!=None):
            for i in t:
                file.write(str(i)+" ")
            file.close()
            del file
    
    def pv(X):
        "O(n*log k + вычисление T)"
        t_in=T(X,other)
        if(t[len(t)-1]<=t_in):
            return 0
        l=0
        r=len(t)-1
        mid=int((l+r)/2)
        while(r-l>1):
            if(t[mid]>t_in):
                r=mid
            else:
                l=mid
            mid=int((l+r)/2)
        return 1- float(l)/len(t)
    return pv




