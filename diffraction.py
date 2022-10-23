from email.headerregistry import ContentTransferEncodingHeader
import math
import cmath
from PIL import Image
import numpy as np
import latexify

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

def waveToIntensity(u2):
    """
    波形データから強度を返す関数
    """
    if len(u2)!=len(u2[0]):
        raise(NotSquareError("u2 size is not square"))
    
    size=len(u2)
    u3=np.full((size,size),0+0j)
    for l in range(size):
        for m in range(size):
            u3[l,m]=abs(u2[l,m])*abs(u2[l,m])
    u3=u3-np.min(u3)
    u3=(u3/np.max(u3))*256
    return u3

def impulse_response(ramuda:float,z:float,p:float,u1):
    """
    波形からフーリエ変換を行う関数

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
            dx=(a-2/N)*p
            dy=(b-2/N)*p
            phase=coefficient*(dx*dx+dy*dy)
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
            u2[l+int(2/N)][m+int(2/N)]=u1[l][m]
            
    return u2

    


class NotSquareError(Exception):
    """正方形で無いときに発火させるエラーです"""
    pass


print(zero_padding([[2,1],[3,1]]))


#ローカルでは動かせないけど、500*500を保存しておく
# u1=np.full((500,500),0+0j)
# u1[250][250]=1

# u2=normal_fresnel_diffraction(pow(10,-5),6*pow(10,-7),10,500,u1)
# u3=waveToIntensity(u2)

# boxelPlate=Image.new("L",(500,500),256)
# for l in range(500):
#     for k in range(500):
#         boxelPlate.putpixel((l,k),int(u3[l,k]))

# boxelPlate.show()