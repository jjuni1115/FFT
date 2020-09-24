from PIL import Image
import numpy as np
from matplotlib import pyplot as plot
import random
import math

def FFT(pixel,row,col,dir):  #fft funtion
    x=0
    y=0
    index=0
    M=0
    N=0
    num=0
    col_data = [0 for i in range(col)]
    row_data = [0 for i in range(row)]
    num=col
    while num>=2:
        num>>=1
        M=M+1
    num=row
    while num>=2:
        num>>=1
        N=N+1
    for y in range(0,row):
        index=y*col
        for x in range(0,col):
            row_data[x]=complex(pixel[index].real,pixel[index].imag)
            index=index+1
        fft(row_data,M,col,dir)
        index=y*col
        for x in range(0,col):
            pixel[index]=complex(row_data[x].real,row_data[x].imag)
            index=index+1

    for x in range(0,512):
        index=x
        for y in range(0,row):
            col_data[y]=complex(pixel[index].real,pixel[index].imag)
            index+=col
        fft(col_data,N,row,dir)
        index=x
        for y in range(0,row):
            pixel[index]=complex(col_data[y].real,col_data[y].imag)
            index+=col
    if dir==1:   #forward  |H(u,v)| operation

        dct_pixel = np.ones((512, 512), dtype=int)

        for i in range(0, 512):
            for j in range(0, 512):
                dct_pixel[i][j] = N*math.log(abs(1+abs(complex(pixel.real[i*512+j],pixel.imag[i*512+j]))))
        return dct_pixel




def fft(f,logN,numpoints,dir):
    scramble(numpoints,f)
    butterflies(numpoints,logN,dir,f)

def scramble(numpoints,f):
    i=0
    j=0
    m=0
    temp=0.0
    for i in range(0,numpoints):
        if(i>j):
            temp=complex(f[j].real,0)
            f[j]=complex(f[i].real,f[j].imag)
            f[i]=complex(temp.real,f[i].imag)
            temp=complex(temp.real,f[j].imag)
            f[j]=complex(f[j].real,f[i].imag)
            f[i]=complex(f[i].real,temp.imag)
        m=numpoints>>1
        while j>=m and m>=2:
            j-=m
            m=m>>1
        j+=m


def butterflies(numpoints,logN,dir,f):
    angle=0.0
    w=0+0j
    wp=0+0j
    temp=0+0j
    i=0
    j=0
    k=0
    offset=0
    N=0
    half_N=0
    wtemp=0.0

    N=1
    for k in range(0,logN):
        half_N=N
        N<<=1
        angle = -2.0*math.pi/N*dir
        wtemp=math.sin(0.5*angle)
        wp=complex(-2.0*wtemp*wtemp,math.sin(angle))
        w=complex(1.0,0.0)
        for offset in range(0,half_N):
            for i in range(offset,numpoints,N):
                j=i+half_N
                temp=complex((w.real*f[j].real)-(w.imag*f[j].imag),(w.imag*f[j].real)+(w.real*f[j].imag))
                f[j]=complex(f[i].real-temp.real,f[i].imag-temp.imag)
                f[i]=complex(f[i].real+temp.real,f[i].imag+temp.imag)
            wtemp=w.real
            w=complex(wtemp*wp.real-w.imag*wp.imag+w.real,w.imag*wp.real+wtemp*wp.imag+w.imag)
    if dir==-1:   #inverse
        for i in range(0,numpoints):
            f[i]=complex(f[i].real/numpoints,f[i].imag/numpoints)


def mean_square_error(pixel,newpixel):    #mean square error function
    sum=0
    mse=0.0
    for i in range(0,512):
        for j in range(0,512):
            sum+=(newpixel[i][j]-pixel[i][j])**2
    mse=sum/(512*512)
    return mse



#load raw file
infile=open("BOAT512.raw","rb")
image=np.fromfile(infile,dtype=np.uint8,count=512*512)
pixel=Image.frombuffer('L',[512,512],image,'raw','L',0,1)
pixel=np.array(pixel,dtype=np.complex)  #행렬에 저장

#call fft funtion
pixel=pixel.reshape(-1,)
new_pixel = np.ones((512, 512), dtype=int)
new_pixel=FFT(pixel,512,512,1)
print("FFT complete")
plot.imshow(new_pixel.real,cmap='gray')
plot.show()

#inverse
FFT(pixel,512,512,-1)
print("inverse complete")
pixel=pixel.reshape(512,512)
plot.imshow(pixel.real,cmap='gray')
plot.show()

#mean square error
infile=open("BOAT512.raw","rb")
image1=np.fromfile(infile,dtype=np.uint8,count=512*512)
original_pixel=Image.frombuffer('L',[512,512],image1,'raw','L',0,1)
original_pixel=np.array(original_pixel,dtype=int)  

mse=mean_square_error(original_pixel,pixel.real)
print("Mean square error: ",mse)

