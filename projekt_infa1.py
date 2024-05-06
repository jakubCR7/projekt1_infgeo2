# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 11:56:02 2024

@author: Uzytkownik
"""
from numpy import *
import sys
def is_float(string):
    try:
        float(string)
        return True
    except ValueError:
        return False
class Transformacje:
    def __init__(self, model: str = 'wgs80'):
        if model == "wgs84":
            self.a = 6378137.0 # semimajor_axis
            self.b = 6356752.31424518 # semiminor_axis
        elif model == "grs80":
            self.a = 6378137.0
            self.b = 6356752.31414036
        elif model == "Krasowski":
            self.a = 6378245.0
            self.b = 4
        else:
            raise NotImplementedError(f"{model} model not implemented")
        self.flat = (self.a - self.b) / self.a
        self.ecc = sqrt(2 * self.flat - self.flat ** 2) # eccentricity  WGS84:0.0818191910428 
        self.ecc2 = (2 * self.flat - self.flat ** 2) # eccentricity**2
        
    def hirvonen(self, X, Y, Z, output = 'dec_degree'):

            # promien rownoleznika
        p = sqrt(X**2 + Y**2)
            # przyblizona wartosc phi
        phi = arctan(Z/(p*(1-self.ecc2)))

        N = self.a/(1 - self.ecc2 * sin(phi)**2)**(1/2) 
            # obliczamy wartosc N
        while True:
            phi_poprzednie = phi
                #print(phi)
                # obliczenie wysokosci
            h = p/cos(phi_poprzednie) - N
                #obliczenie poprawionej wartosci phi
            phi = arctan(Z/p*((N*(1-self.ecc2)+h)/(N+h))**(-1))
            #warunek
            if abs(phi -phi_poprzednie)<(0.000001/3600*pi/180):
                break
            N= self.a/(1-self.ecc2*sin(phi)**2)**(1/2)
            h= p/cos(phi)-N
            lam= arctan2(Y,X)
        if output == "dec_degree":
            return degrees(phi), degrees(lam), h 
        elif output == "dms":
            phi = self.deg2dms(degrees(phi))
            lam = self.deg2dms(degrees(lam))
            return f"{lat[0]:02d}:{lat[1]:02d}:{lat[2]:.2f}", f"{lon[0]:02d}:{lon[1]:02d}:{lon[2]:.2f}", f"{h:.3f}"
        else:
            raise NotImplementedError(f"{output} - output format not defined")
        
    def phi_lam_XYZ(self, phi, lam, h, output = 'dec_degree'):
        phi = phi/180 * pi
        lam = lam/180* pi
        N = self.a/(1 - self.ecc2 * sin(phi)**2)**(1/2)
        X = (N+h) * cos(phi) * cos(lam)
        Y = (N+h) * cos(phi) * sin(lam)
        Z = (N*(1-self.ecc2)+h) * sin(phi)
        return X, Y, Z
    
    def transformacja_neu(self, X, Y, Z, phi, lam):
        R = array([[-sin(phi) * cos(lam), -sin(lam), cos(phi) * cos(lam)],
                      [-sin(phi) * sin(lam), cos(lam), cos(phi) * sin(lam)],
                      [cos(phi), 0 , sin(phi)]]) 
        XYZ = array([X, Y, Z])
        Rneu1 = R
        neu = Rneu1.T @ XYZ
        
        return neu
    
    def Gk2000(self, lam, phi, lam_0, nr):
        lam = lam/180 * pi
        phi = phi/180 * pi
        lam_0 = lam_0/180 * pi
        
        b2 = self.a**2 * (1-self.ecc2)
        e2p = (self.a **2 - b2)/b2
        n2 = e2p * cos(phi)**2
        t = tan(phi)
        dlam = lam - lam_0
        N = self.a/(1 - self.ecc2 * sin(phi)**2)**(1/2)
        A0 = 1 - (self.ecc2/4) - ((3 * self.ecc2**2)/64) - ((5 * self.ecc2 * self.ecc2**2)/256)
        A2 = 3/8 * (self.ecc2 + self.ecc2**2/4 + 15 * self.ecc2**3/128)
        A4 = 15/256 * (self.ecc2**2 + (3 * self.ecc2 * self.ecc2**2)/4)
        A6 = (35 * self.ecc2 * self.ecc2**2)/3072
        sigma = self.a * (A0 * phi - A2 * sin(2 * phi) + A4 * sin(4 * phi) - A6 * sin(6 * phi))
        xgk = sigma + (dlam**2)/2 * N * sin(phi) * cos(phi) * (1 + (dlam**2)/12 * cos(phi)**2 * (5 - t**2 + 9 * n2 + 4 * n2**2) + (dlam**4)/360 * cos(phi)**4 * (61 - 58 * t**2 + t**4 + 270 * n2 - 330 * n2 * t**2))
        ygk = dlam * N * cos(phi) * (1 + (dlam**2)/6 * cos(phi)**2 * (1 - t**2 + n2) + (dlam**4)/120 * cos(phi)**4 * (5 - 18 * t**2 + t**4 + 14 * n2 - 58 * n2 * t**2))
        x2000 = xgk * 0.999923
        y2000 = ygk * 0.999923 + 500000 + nr * 1000000
        return x2000, y2000
    
    def Gk92(self, lam, phi):
        lam_0 = 19
        lam_0 = lam_0/180 *pi
        lam = lam/180 *pi
        phi = phi /180 * pi
        dlam = lam - lam_0
        b2 = self.a**2 * (1-self.ecc2)
        e2p = (self.a **2 - b2)/b2
        n2 = e2p * cos(phi)**2
        t = tan(phi)
        N = self.a/(1 - self.ecc2 * sin(phi)**2)**(1/2)
        A0 = 1 - (self.ecc2/4) - ((3 * self.ecc2**2)/64) - ((5 * self.ecc2 * self.ecc2**2)/256)
        A2 = 3/8 * (self.ecc2 + self.ecc2**2/4 + 15 * self.ecc2**3/128)
        A4 = 15/256 * (self.ecc2**2 + (3 * self.ecc2 * self.ecc2**2)/4)
        A6 = (35 * self.ecc2 * self.ecc2**2)/3072
        sigma = self.a * (A0 * phi - A2 * sin(2 * phi) + A4 * sin(4 * phi) - A6 * sin(6 * phi))
        xgk = sigma + (dlam**2)/2 * N * sin(phi) * cos(phi) * (1 + (dlam**2)/12 * cos(phi)**2 * (5 - t**2 + 9 * n2 + 4 * n2**2) + (dlam**4)/360 * cos(phi)**4 * (61 - 58 * t**2 + t**4 + 270 * n2 - 330 * n2 * t**2))
        ygk = dlam * N * cos(phi) * (1 + (dlam**2)/6 * cos(phi)**2 * (1 - t**2 + n2) + (dlam**4)/120 * cos(phi)**4 * (5 - 18 * t**2 + t**4 + 14 * n2 - 58 * n2 * t**2))
        x92 = xgk * 0.9993 - 5300000
        y92 = ygk * 0.9993 + 500000
        return x92, y92
    
        
if __name__ == "__main__":
    # utworzenie obiektu
    print('podaj nazwe elipsoidy')
    print('wybierz jedno z: \n grs80\n wgs84 \n Krasowski')
    geo = Transformacje(model = input(str))
    print('wybierz rodzaj transformacji (wpisz odpowiednia cyfre):\n 1 = X, Y, Z --> phi, lam, h \n 2 = phi, lam, h --> X, Y, Z \n 3 = X, Y, Z --> neu \n 4 = BL --> X2000, Y2000 \n 5 = BL --> X92, Y92')
    info = input(str)
    f = open("phi_lam_h.txt")
    dane = f.read()
    dane = dane.split('\n')
    lista =[0,0,0]
    for i in dane:
        w = i.split(',')
        listaX = []
        for j in w:
            if is_float(j.strip('')) == True:
                listaX.append(float(j))
            else:
                listaX = [0,0,0]
        lista = vstack((lista,listaX))
    lista = array(lista)
    f.close()
    if info == '1':
        X = lista[:,0]; Y = lista[:,1]; Z = lista[:,2]
        t=0
        lista_phi = []
        lista_lam = []
        lista_h = []
        while t < len(X):
            if X[t] != 0:
                phi, lam, h = geo.hirvonen(X[t], Y[t], Z[t])
                lista_h.append(h)
                lista_lam.append(lam)
                lista_phi.append(phi)
            t = t + 1
        lista_phi = array(lista_phi)
        lista_lam = array(lista_lam)
        lista_h = array(lista_h)
        dane_wyjsc = column_stack((lista_phi,lista_lam,lista_h))
        f = open("phi_lam_h.txt",'w')
        f.write(f'  phi     lam     h \n')
        f.write(f'------------------------\n')
        for i in dane_wyjsc:
            for j in i:
                f.write(f'{j:.9f}')
                if j == i[-1]:
                    break
                f.write(f',')
            f.write(f'\n')
        f.close()
        print('program zapisał plik o nazwie "phi_lam_h.txt" ze współrzędnymi geodezyjnymi w bierzącym folderze')
    elif info == '2':
        if info == '2':
            phi = lista[:,0]; lam = lista[:,1]; h = lista[:,2]
            t=0
            lista_x = []
            lista_y = []
            lista_z = []
            while t < len(phi):
                if phi[t] != 0:
                    x, y, z = geo.phi_lam_XYZ(phi[t],lam[t],h[t])
                    lista_x.append(x)
                    lista_y.append(y)
                    lista_z.append(z)
                t = t + 1
            lista_x = array(lista_x)
            lista_y = array(lista_y)
            lista_z = array(lista_z)
            dane_wyjsc = column_stack((lista_x,lista_y,lista_z))
            r = open("X_Y_Z.txt",'w')
            r.write(f'   X            Y            Z \n')
            r.write(f'------------------------------------\n')
            for i in dane_wyjsc:
                for j in i:
                    r.write(f'{j:.3f}')
                    if j == i[-1]:
                        break
                    r.write(f',')
                r.write(f'\n')
            r.close()
    elif info == '3':
        if info == '3':
            X = lista[:,0]; Y = lista[:,1]; Z = lista[:,2]
            t=0
            lista_phi = []
            lista_lam = []
            while t < len(phi):
                if phi[t] != 0:
                    x, y, z =geo.hirvonen(X[t],Y[t],Z[t])
                    lista_x.append(x)
                    lista_y.append(y)
                    lista_z.append(z)
                t = t + 1
            lista_x = array(lista_x)
            lista_y = array(lista_y)
            lista_z = array(lista_z)
            dane_wyjsc = column_stack((lista_x,lista_y,lista_z))
            r = open("X_Y_Z.txt",'w')
            r.write(f'   X            Y            Z \n')
            f.write(f'------------------------------------\n')
            for i in dane_wyjsc:
                for j in i:
                    r.write(f'{j:.3f}')
                    if j == i[-1]:
                        break
                    r.write(f',')
                r.write(f'\n')
            r.close()
        X = 3664940.500; Y = 1409153.590; Z = 5009571.170; phi= 52.6278; lam=21.029876
        neu = geo.transformacja_neu(X,Y,Z, phi, lam)
        print(neu)
    elif info == '4':
        if info == '4':
            phi = lista[:,0]; lam = lista[:,1]
            t=0
            lista_x2000 = []
            lista_y2000 = []
            while t < len(phi):
                lam_0 = 21
                nr = 7
                if phi[t] != 0:
                    if float(lam[t])<16.5:
                        lam_0 = 15
                        nr=5
                    elif float(lam[t])<19.5:
                        lam_0 = 18
                        nr=6
                    elif float(lam[t])<24.5:
                        lam_0=21
                        nr = 7
                    elif float(lam[t])<25.5:
                        lam_0 = 24
                        nr = 8
                    x, y = geo.Gk2000(lam[t], phi[t], lam_0, nr)
                    lista_x2000.append(x)
                    lista_y2000.append(y)
                t = t + 1
            lista_x2000 = array(lista_x2000)
            lista_y2000 = array(lista_y2000)
            dane_wyjsc = column_stack((lista_x2000,lista_y2000))
            r = open("X_Y_2000.txt",'w')
            r.write(f'   X            Y  \n')
            r.write(f'-----------------------------\n')
            for i in dane_wyjsc:
                for j in i:
                    r.write(f'{j:.3f}')
                    if j == i[-1]:
                        break
                    r.write(f',')
                r.write(f'\n')
            r.close()
            
    elif info == '5':
        if info == '5':
            phi = lista[:,0]; lam = lista[:,1]
            t=0
            lista_x92 = []
            lista_y92 = []
            while t < len(phi):
                if phi[t] != 0:
                    x, y = geo.Gk92(lam[t], phi[t])
                    lista_x92.append(x)
                    lista_y92.append(y)
                t = t + 1
            lista_x92 = array(lista_x92)
            lista_y92 = array(lista_y92)
            dane_wyjsc = column_stack((lista_x92,lista_y92))
            r = open("X_Y_92.txt",'w')
            r.write(f'   X            Y  \n')
            for i in dane_wyjsc:
                for j in i:
                    r.write(f'{j:.3f}')
                    if j == i[-1]:
                        break
                    r.write(f',')
                r.write(f'\n')
            r.close()






















