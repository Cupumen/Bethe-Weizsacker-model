#!/usr/bin/env python
# coding: utf-8

# In[1]:


#Muhamad Rizki Septiawan (G74160006)
#The code was made for my task in Nuclear Physisc Lecture

#importing the libraries

import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# In[2]:


df = pd.read_csv('Periodic Table of Elements.csv') #load data
df.head() #show 5 data
#df #to show all data


# In[3]:


#Bethe-Weizsacker Parametrization for Nuclear Models

def B_BW(Z,N):
    A = Z+N
    if (A%2==1):
        #print('case1')
        B = 15.753*A-17.804*math.pow(A, 2/3)-0.7103*(math.pow(Z,2)/math.pow(A,1/3))-23.69*math.pow((N-Z),2)/A+0
    #print (B/A)
    elif (A%2==0):
        if (N%2==0 and Z%2==0):
            #print('case2')
            B = 15.753*A-17.804*math.pow(A, 2/3)-0.7103*(math.pow(Z,2)/math.pow(A,1/3))-23.69*math.pow((N-Z),2)/A +33.6*math.pow(A,(-3/4))
        elif (N%2==1 and Z%2==1):
            #print('case3')
            B = 15.753*A-17.804*math.pow(A, 2/3)-0.7103*(math.pow(Z,2)/math.pow(A,1/3))-23.69*math.pow((N-Z),2)/A -33.6*math.pow(A,(-3/4))
    return (B/A);


# In[4]:


def BE(Z,N):
    A=Z+N;
    mn = 938.27
    mp = 939.56;
    Be = -((Z*mp+N*mn)-A*mp)
    return 10*Be/A


# In[5]:


#make array
be=[]; nucleon=[];         #for Bethe-Weizsacker model
be2=[]; nucleon2=[];       #for data
    
for x in range(1,117):
    Z = df.loc[x,'AtomicNumber'];
    N = df.loc[x,'NumberofNeutrons'];
    
    be.insert(x, B_BW(Z, N));
    nucleon.insert(x,(Z+N));
    be2.insert(x, BE(Z, N));
    nucleon2.insert(x,(Z+N))
#make numpy array    
x1=np.array(nucleon); y1=np.array(be);   
x2=np.array(nucleon2); y2=np.array(be2);


# In[6]:


#plotting model
#plt.plot(nucleon2, be2, '.r', nucleon, be, '.b');
plt.plot(x1, y1, '.b');
#plt.legend(['Data', 'Bethe-Weizsacker Models'])
plt.xlabel('A');
plt.ylabel('B/A (MeV)');
plt.show()


# In[7]:


plt.plot(x2, y2, '.r');
plt.xlabel('A');
plt.ylabel('B/A (MeV)');
plt.show(1)


# In[8]:


#plotting model
plt.plot(nucleon2, be2, '.r', nucleon, be, '.b');
plt.legend(['Data', 'Bethe-Weizsacker Models'])
plt.xlabel('A');
plt.ylabel('B/A (MeV)');
plt.show(2)


# In[9]:


fig, axs = plt.subplots(3,sharex=True)
fig.suptitle('')
axs[0].plot(x1, y1, '.b')
axs[1].plot(x2, y2, '.r')
axs[2].plot(x1, y1, '.b', x2, y2, '.r')
axs[1].set(ylabel='B/A (MeV)')
axs[2].set(xlabel='A')
plt.show(3)

