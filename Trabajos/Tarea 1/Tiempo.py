import time as tm
import os 
import matplotlib.pyplot as plt
import numpy as np

def get_sec(time_str):
    m, s = time_str.split(':')
    return int(m) * 60 + float(s)

tm1 = os.popen('time ./tarea1.out 1 >data.txt 2> tiempo1.txt')
tm.sleep(100)
temp1 = open('tiempo1.txt')
tmp1 = temp1.read()
tm1.close()

tm2 = os.popen('time ./tarea1.out 2 >borrar.txt 2> tiempo2.txt')
tm.sleep(100)
temp2 = open('tiempo2.txt')
tmp2 = temp2.read()
tm2.close()

tm4 = os.popen('time ./tarea1.out 4 >borrar.txt 2> tiempo4.txt')
tm.sleep(100)
temp4 = open('tiempo1.txt')
tmp4 = temp4.read()
tm4.close()




t1=get_sec(tmp1.partition("system")[2].partition("elapsed")[0])
t2=get_sec(tmp2.partition("system")[2].partition("elapsed")[0])
t3=get_sec(tmp4.partition("system")[2].partition("elapsed")[0])
temp=[t1,t2,t3]
proc=[1,2,4]

plt.scatter(proc,temp,marker='^', s=40, c='black')
plt.xlabel("$numero\ de\ procesadores $", size=14)
plt.ylabel("$tiempo\  (s)$",size=14)
plt.xlim(0,5)
plt.savefig("tiempo_vs_procesadores.pdf")

