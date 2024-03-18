import numpy as np
import matplotlib.pyplot as plt
import scipy 
import ROOT as root

# Funcion Uniformidad Integral.
# Input: N_pix, C_med.
# Output: unif_integral (escalar)



def field_of_view(N_pix, filas_eliminadas, columnas_eliminadas):
    fov = [N_pix - filas_eliminadas, N_pix - columnas_eliminadas]
    return fov


def matriz_uniform(N_pix, C_med, field_of_view):                                                   
    matriz_cuentas = np.zeros((N_pix,N_pix))
    filas_ini = int(N_pix - field_of_view[0])
    filas_fin = int(field_of_view[0])
    col_ini = int(N_pix - field_of_view[1])
    col_fin =int(field_of_view[1])
    
    for i in range(filas_ini, filas_fin):
        for j in range(col_ini, col_fin):
            matriz_cuentas[i,j] = np.random.normal(C_med, np.sqrt(C_med))
    return matriz_cuentas

def filtro(matriz_cuentas,field_of_view,N_pix):
    matriz_filtrada = scipy.ndimage.uniform_filter(matriz_cuentas)
    filas_ini = int(N_pix - field_of_view[0]+1)
    filas_fin = int(field_of_view[0]-1)
    
    col_ini = int(N_pix - field_of_view[1]+1)
    col_fin =int(field_of_view[1]-1)
    matriz_reducida = np.zeros([filas_fin-filas_ini,col_fin-col_ini])
    for i in range(filas_fin-filas_ini):
        for j in range(col_fin-col_ini):
            matriz_reducida[i,j]=matriz_filtrada[filas_ini + i,col_ini + j]
    return matriz_reducida

def filtro_9p(matriz_cuentas,N_pix):
    matriz_filtrada = np.zeros((N_pix,N_pix))
    for i in range(4,N_pix-4):
        for j in range(12,N_pix-12):
            P1=matriz_cuentas[i-1,j-1]
            P2=matriz_cuentas[i-1,j]
            P3=matriz_cuentas[i-1,j+1]
            P4=matriz_cuentas[i,j-1]
            P5=matriz_cuentas[i,j]
            P6=matriz_cuentas[i,j+1]
            P7=matriz_cuentas[i+1,j-1]
            P8=matriz_cuentas[i+1,j]
            P9=matriz_cuentas[i+1,j+1]
            matriz_filtrada[i,j]=(P1*1 + P2*2 + P3*1 +2*P4 + 4*P5 + 2*P6 + 1*P7 + 2*P8 + 1*P9)/(16)
    "Aqui extraigo la submatriz para continuar con el analisis"
    start_row = 4
    size_row = N_pix-8
    start_col=12
    size_col=N_pix-24
    submatrix = matriz_filtrada[start_row:start_row + size_row, start_col:start_col + size_col]
    return submatrix


    


def uniformidad_integral(matriz_cuentas):
    max = np.max(matriz_cuentas)
    min = np.min(matriz_cuentas)
    unif_integral = 100*(max - min)/(max + min)
    return unif_integral

def analisis(N_pix, C_med,filas_eliminadas, columnas_eliminadas):
    fov = field_of_view(N_pix, filas_eliminadas, columnas_eliminadas)
    matriz_cuentas = matriz_uniform(N_pix,C_med,fov)
    matriz_filtrada = filtro_9p(matriz_cuentas,N_pix)
    unif_integral = uniformidad_integral(matriz_filtrada)
    return unif_integral

#Generar distribucion Histograma sin defecto
def dist(N_pix,C_med,num_cuentas,filas_eliminadas,columnas_eliminadas):
    d = np.zeros(num_cuentas)
    for i in range(num_cuentas):
        d[i] = analisis(N_pix, C_med,filas_eliminadas, columnas_eliminadas)
    distribucion = np.sort(d)
    return distribucion



#Introduciendo los defectos en las matrices sin defectos 
def matriz_defectos(t,k,N_pix,C_med,fov):
    matriz = matriz_uniform(N_pix,C_med,fov)
    n=int(t*(N_pix/64))
    for i in range(int(N_pix/2),int(N_pix/2)+n):
        for j in range(int(N_pix/2),int(N_pix/2)+n):
             matriz[i,j] = np.random.normal((1+k/100)*C_med, np.sqrt((1+k/100)*C_med))

    return matriz

def analisis_defectos(t,k,N_pix, C_med,filas_eliminadas, columnas_eliminadas):
    fov = field_of_view(N_pix, filas_eliminadas, columnas_eliminadas)
    matriz_cuentas = matriz_defectos(t,k,N_pix,C_med,fov)
    matriz_filtrada = filtro_9p(matriz_cuentas,N_pix)
    unif_integral = uniformidad_integral(matriz_filtrada)
    return unif_integral

def dist_defectos(t,k,N_pix,C_med,num_cuentas,filas_eliminadas,columnas_eliminadas):
    d = np.zeros(num_cuentas)
    for i in range(num_cuentas):
        d[i] = analisis_defectos(t,k,N_pix,C_med,filas_eliminadas,columnas_eliminadas)
    distribucion = np.sort(d)
    return distribucion

#incertidumbres de un histograma
"400 bins a 0.025 de 0.0 a 10.0"

def incertidumbre(histograma,num_cuentas): #histograma debe estar normalizado
    freq=np.zeros(400)
    sigma=np.zeros(400)
    for i in range(400):
        freq[i]=histograma.GetBinContent(i)
        if freq[i]<1.0:
           sigma[i]=np.sqrt((1/num_cuentas)*freq[i]*(1-freq[i]))
        else:
            sigma[i]=0.0
    return sigma

def name(t,k,N_pix,nct,C_med):
    return "d_"+"t="+str(t)+"_k="+str(k)+"_"+str(N_pix)+"P_"+str(nct)+"Mc_"+"Cmed_"+str(C_med)

def in_defectos(t,k,C_med,matriz_sana):
    N_pix = np.shape(matriz_sana)[0]
    matriz = matriz_sana
    n=int(t*(N_pix/64))
    for i in range(int(N_pix/2),int(N_pix/2)+n):
        for j in range(int(N_pix/2),int(N_pix/2)+n):
             matriz[i,j] = np.random.normal((1+k/100)*C_med, np.sqrt((1+k/100)*C_med))
    return matriz

'''funcion para nombrar a los hitogramas: xxx(num_pixel)-xx(mill.cuentas)-fx(filtro, 0=no filtro)-xx(tamaño, 0=no hay defecto)-xx(contraste*10)
    nombre = 064-80-f1-2-10 significa 64pix 80Mc filtro9p t=2 k=1.0
    Cmed no es un parametro independiente, se calcula a partir de 80Mc y 64pix'''

def histogram_name(N_pix,nct,filtro,t=0,k=0):
    N_pix = str(int(N_pix))
    if len(N_pix)<3:
        N_pix = '0'+N_pix
    
    nct = str(int(nct))
    if (nct=="5"): nct="0"+nct
    filtro = str(int(filtro))
    
    if t==0:k=0
    if k==0:t=0
    t = str(t)
    if int(k)==10:
        k="100"
    else:
        k= str(int(k*10))
    name = N_pix +'-'+nct+ '-f' + filtro + '-' + t + '-' + k + ".root"
    return name
