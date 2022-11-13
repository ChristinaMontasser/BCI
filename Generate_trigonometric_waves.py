#Task1 _ BCI

import math
import numpy as np
import turtle
import matplotlib.pyplot as plt

#any trigmotric signal is on the form of a*trig(2*pi*f*t+phase_shift)
#time is sampled from the continous real time, a vector.
a = int(input("Enter the magntitude: ")) 
fs= int(input("Input the sampling rate: "))
maximun_time= int(input("Input maximum time: "))
ti= np.arange(0, maximun_time+(1/fs), 1/fs)
f = int(input("Input frequency: ")) 
phi = int(input("Input the shifted number (if any): ")  )

signal_sine = [a*math.sin(math.radians((2*180*f*t_i + phi))) for t_i in ti]
signal_cosine = [a*math.cos(math.radians((2*180*f*t_i + phi))) for t_i in ti]


#An animated plot  
win = turtle.Screen()
win.bgcolor("white")

#Setting the window, vertically and horizontally 
                #x_neg, y_neg, x_pos, y_pos
win.setworldcoordinates(0, -2*a, maximun_time, 2*a)
tur = turtle.Turtle()
 
tur.goto(0, 2*a)
tur.goto(0, -2*a)
tur.goto(0, 0)
tur.goto(maximun_time, 0)
tur.penup()

tur.goto(0, signal_sine[0])
tur.pendown()
tur.pencolor("blue")
tur.pensize(4)

#drawing sine
for x in range(maximun_time*fs):
    tur.goto(ti[x], signal_sine[x])

tur.penup()
tur.goto(0, signal_cosine[0])
tur.pendown()
tur.pencolor("red")

#drawing cosine
for x in range(maximun_time*fs):
    tur.goto(ti[x], signal_cosine[x])


#A fixed plot 
plt.plot(ti, signal_sine, label='sin', color='blue')
plt.plot(ti, signal_cosine, label='cos', color ='red')

plt.title('Sine and Cosine waves')
plt.xlabel('Time')
plt.ylabel('Amplitude')

plt.grid(True, which='both')
plt.legend(title="", labelcolor='linecolor')
plt.axhline(y=0, color='k')
plt.show()





