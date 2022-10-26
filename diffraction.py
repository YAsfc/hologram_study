from email.headerregistry import ContentTransferEncodingHeader
import math
import cmath
from site import makepath
import wave
from PIL import Image
import numpy as np
import latexify
from decimal import Decimal, ROUND_HALF_UP, ROUND_HALF_EVEN


@latexify.with_latex
def normal_fresnel_diffraction(p:float,ramuda:float,z:float,u1):
    """
    フーリエ変換を用いないフレネル回折の実装

    Parameters
    ----------
    p : float
        画素ピッチ(1~100um)
    ramuda : float
        光の波長(400~800nm)
    z:float
        伝搬面と伝搬先の距離
        フレネルなのでz>>0.2m
    N:int
        伝搬面と伝搬先の格子数
    u1: complex[N][N]
        伝搬面の波

    Returns
    -------
    u2: complex[N][N]
        伝搬先の波
    """
    if len(u1)!=len(u1[0]):
        raise(NotSquareError("u1 size is not square"))

    N=len(u1)

    u2=np.full((N,N),0+0j)
    coefficient=1j*math.pi/(ramuda*z)

    for e in range(N):
        for f in range(N):
            for g in range(N):
                for h in range(N):
                    dx=((e-N/2)-(g-N/2))*p  
                    dy=((f-N/2)-(h-N/2))*p
                    arg=coefficient*(dx*dx+dy*dy)
                    temp=cmath.exp(arg)
                    u2[g,h]+=temp

    return u2

def wave_to_intensity(u2):
    """
    波形データから強度を返す関数
    """
    if len(u2)!=len(u2[0]):
        raise(NotSquareError("u2 size is not square"))
    
    size=len(u2)
    u3=np.full((size,size),0+0j)
    for l in range(size):
        for m in range(size):
            temp=abs(u2[l,m])*abs(u2[l,m])
            # print(temp)
            u3[l,m]=Decimal(str(temp)).quantize(Decimal('0.00001'), rounding=ROUND_HALF_UP) 
            print(u3[l,m])
            # print(round(temp.real,-6))
    u3=u3-np.min(u3)
    if((np.max(u3)-np.min(u3))!=0):
        u3=(u3/(np.max(u3)-np.min(u3)))*256.0
    else:
        for l in range(size):
            for m in range(size):
                u3[l,m]=256
    return u3

def fourier_fresnel_diffraction(ramuda:float,z:float,p:float,u1):
    """
    与えられたデータに対して、フーリエ変換を用いた用いたフレネル回折を返す関数

    Parameters
    ----------
    p : float
        画素ピッチ(1~100um)
    ramuda : float
        光の波長(400~800nm)
    z:float
        伝搬面と伝搬先の距離
        フレネルなのでz>>0.2m
    u1: complex[N][N]
        伝搬面の波
    Returns
    -------
    u2: complex[N][N]
        伝搬先の波
    """

    if len(u1)!=len(u1[0]):
        raise(NotSquareError("u1 size is not square"))
    N=len(u1)

    #fft_shiftとフーリエ変換の実行を合わせて行う
    
    u1_zeropadded=zero_padding(u1)
    u1_shifted=np.fft.fftshift(u1_zeropadded)
    u1_fouried=np.fft.fft2(u1_shifted)
    
    h2=impulse_response(ramuda,z,p,zero_padding(u1))
    h2_shifted=np.fft.fftshift(h2)
    h2_fouried=np.fft.fft2(h2_shifted)

    mul=u1_fouried*h2_fouried
    u2_reversed=np.fft.ifft2(mul)
    u2_unshifted=np.fft.fftshift(u2_reversed)
    makePhoto(u2_unshifted)

    return u2_unshifted


def impulse_response(ramuda:float,z:float,p:float,u1):
    """
    与えられたデータに対して、インパルス応答を返す関数

    Parameters
    ----------
    p : float
        画素ピッチ(1~100um)
    ramuda : float
        光の波長(400~800nm)
    z:float
        伝搬面と伝搬先の距離
        フレネルなのでz>>0.2m
    u1: complex[N][N]
        伝搬面の波
    Returns
    -------
    h2: complex[N][N]
        伝搬先の波
    """
    if len(u1)!=len(u1[0]):
        raise(NotSquareError("u1 size is not square"))

    N=len(u1)
    h2=np.full((N,N),0+0j)
    
    coefficient=(1j*math.pi)/(ramuda*z)
    for a in range(N):
        for b in range(N):
            dx=(a-N/2)*p
            dy=(b-N/2)*p
            dis=(dx*dx+dy*dy)
            phase=coefficient*dis
            res=cmath.exp(phase)
            h2[a,b]=res
    return h2
            
def zero_padding(u1):
    """
    ゼロパディングを行う関数
    """
    if len(u1)!=len(u1[0]):
        raise(NotSquareError("u1 size is not square"))
    
    N=len(u1)
    u2=np.full((2*N,2*N),0+0j)
    for l in range(N):
        for m in range(N):
            u2[l+int(N/2)][m+int(N/2)]=u1[l][m]
            
    return u2

def makePhoto(u1):
    plate=Image.new("L",(len(u1),len(u1)),256)
    size=len(u1)
    for l in range(size):
        for k in range(size):
            plate.putpixel((l,k),int(u1[l,k]))
    plate.show()


class NotSquareError(Exception):
    """正方形で無いときに発火させるエラーです"""
    pass




u1=np.full((100,100),0+0j)

u1[49][49]=255
u2=fourier_fresnel_diffraction(6*pow(10,-7),10,pow(10,-4),u1)


makePhoto(wave_to_intensity(u2))
# u3=wave_to_intensity(u2)