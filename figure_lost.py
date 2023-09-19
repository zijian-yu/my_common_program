#!/usr/bin/env python
# -*- coding: utf-8 -*-
from sympy import * #调用处理数字符号的sympy库
from matplotlib.pyplot import * #调用绘图matplotlib库
from numpy import *
from scipy.optimize import curve_fit
import sys
from scipy import stats
figure(figsize=(8,5))
root = axes([0, 0, 1, 1])
f1=open(sys.argv[1])
data1=[]
data2=[]
data3=[]
num=0
for line in f1.readlines():
    a=[]
    a=line.split()
    if(0 < num < 16):
        data1.append(int(a[0]))
        data2.append(int(a[1]))
    num+=1
f1.close()
print(data1)
print(data2)
def func(x,p,k):
    return (1-p)**(x-1)*p*k
popt,pcov = curve_fit(func,data1,data2)
print(popt[0])
def rsquared(x, y):
    """ Return R^2 where x and y are array-like."""
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    return r_value**2
yhat=[]
for i in data1:
    yhat.append(func(i,popt[0],popt[1]))                         # or [p(z) for z in x]
print(yhat)
res=rsquared(np.array(yhat),np.array(data2))
print(res)
for i in data1:
    data3.append(func(i,float(popt[0]),float(popt[1])))
geom=[]
top=0
X=[]
Y1=[]
Y2=[]
for i in range(0,100000):
    if(max(data2)>top):
        top+=200
    else:
        break
print(top)
x_self=0.8
y_self=0.7
x_s=0.11
y_s=0.15
y_other=0.02
p_y=y_self/top
for i in range(0,int(top/100)+1):
    if i%2==0:
        text(x_s-0.05,y_s+i*100*p_y-10*p_y,int(i)*100,color='k',fontsize=9)
    plot([x_s-0.005,x_s],[y_s+i*100*p_y,y_s+i*100*p_y],color='k',lw=2)

p_x=x_self/(len(data1))
print(min(data1),max(data1))
align = dict(family='Times New Roman',style='italic')
for i in range(0,len(data1)):
    #print sum(data)
    #print float(data2[i])/sum(data)
    x_true=data1[i]*p_x+x_s-0.5*p_x
    y_true=data2[i]*p_y+y_s
    yline=data3[i]*p_y+y_s
    print(x_true,y_true,yline)
    X.append(x_true+p_x/6)
    Y1.append(y_true)
    Y2.append(yline)
    for j in range(0,20):
	    x0=p_x/60
	    plot([x_true+j*x0,x_true+j*x0],[y_s,y_true],'b',linewidth = 1)
    text(x_true,y_s-0.03,data1[i],color="black",fontsize=10)	
plot([x_s,x_s],[y_s,y_s+y_self],'k',linewidth = 2)
plot([x_s,x_s+x_self],[y_s,y_s],'k',linewidth = 2)
plot([x_s,x_s+x_self],[y_s+y_self,y_s+y_self],'k',linewidth = 2)
plot([x_s+x_self,x_s+x_self],[y_s,y_s+y_self],'k',linewidth = 2)
# use pylab to plot x and y : Give your plots names
plot(X,Y2,linewidth = 2,color = 'red')
#title('Plot of y vs. x')# give plot a title
align = dict(family='Times New Roman',style='italic',fontsize=14)
text(0.45,0.07,"",color = "black",**align)# make axis labels
text(0.03,0.45,'',color = "black", rotation = 90,**align)
#text(0.35,0.5,'p = 0.27',color='black',**align)
root.set_xlim(0,1)
root.set_ylim(0,1)
root.set_axis_off()
savefig("lost.png", dpi=600)
rcdefaults()
show()
