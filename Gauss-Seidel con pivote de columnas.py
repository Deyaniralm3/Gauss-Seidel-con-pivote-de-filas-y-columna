# -*- coding: utf-8 -*-
"""
Created on Fri Oct  8 10:17:52 2021

@author: Iraide
"""

import numpy as np

def read_input (text_basis):
    a_temp, b_temp =text_basis.strip().split('=')
    b_temp=eval(b_temp.strip(' '))
    a_temp=a_temp[:-1]
    a_temp=[eval(i) for i in a_temp.split()]
    return a_temp,b_temp

def read_file(ruta):
    A=[]
    b=[]
    with open (ruta, 'r') as f:
        flag=0
        for line in f:
            if line.strip()!='X':
                pass
            else:
                flag=1
                continue
            if line == 'Y':
                break
            if flag==1:
                aux_1, aux_2=read_input(line)
                A.append(aux_1)
                b.append(aux_2)
    return A,b

ruta1='C:/Users/Iraide/Documents/sistema4.txt'
            
print ("A=",read_file(ruta1)[0])  
print ("b=",read_file(ruta1)[1])  

def pivoteo_col(A,b):
    tam=len(b)
    A_n= np.array(A)
    l=[0.0 for x in range(tam)]
    s=[0.0 for x in range(tam)] 
   
    for i in range (tam):
        l[i]=i
        smax=0.0
        for j in range(tam):
            if abs(A[j][i])>smax:
                smax=abs(A[j][i])
        s[i]=smax
        
    for i in range (tam-1):
        rmax=0.0
        for j in range(i,tam):
            r=abs(A[i][l[j]])/s[l[j]]
            if r>rmax:
                rmax=r
                rindex=j
        if(rindex!=0):       
            temp=np.copy(A_n[:,i])
            A_n[:,i]=A_n[:,rindex]
            A_n[:,rindex]=temp
            
            temp1=b[i]
            b[i]=b[rindex]
            b[rindex]=temp1
        
        
    return A_n,b
    

print ("A2=",pivoteo_col(read_file(ruta1)[0],read_file(ruta1)[1])[0])  
print ("b2=",pivoteo_col(read_file(ruta1)[0],read_file(ruta1)[1])[1])
        
def gauss_seidelnumphy_conlect_pivocol(ruta, umbral, max_iter):
    A=read_file(ruta)[0]
    b=read_file(ruta)[1]
    m=len(A[0])
    n=len(A)
    tam=len(b)
    A2=pivoteo_col(A,b)[0]
    b2=pivoteo_col(A,b)[1]
    if m==n:
        A_np = A2
        b_np = np.array(b2)
        x_np = np.zeros(tam) 
        aux_np = np.ones(tam)
        for ite in range(max_iter):
            for i in range(tam):
                aux_np[i] = 0.0
                x_np[i] = (b_np[i] - np.sum(x_np*aux_np*A_np[i,:]))/A_np[i][i]
                aux_np[i] = 1.0
            
            current_b = np.dot(A_np,x_np)
            error = np.sum(np.abs(current_b-b_np))
            
            if error < umbral:
                return x_np
                break
    else:
        return ("La matriz no es cuadrada, no puedo encontrar solución")    


print("Las soluciones son:", gauss_seidelnumphy_conlect_pivocol(ruta1, 0.0001, 100000))