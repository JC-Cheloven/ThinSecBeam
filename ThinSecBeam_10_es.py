# -*- coding: utf-8 -*-

#  ThinSecBeam 1.0
#  
#  Copyright 2018: Juan Carlos del Caño
#                  Prof. at the "Escuela de Ingenierias Industriales"
#                  University of Valladolid (Spain)
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  

# Importacion de librerias 
import numpy as np
from tkinter import *
from tkinter import ttk, messagebox, filedialog

import matplotlib
from matplotlib import pyplot as plt
matplotlib.use('TkAgg')  # mac requiere hacer esto explicitamente
from matplotlib.patches import Arc
from mpl_toolkits import mplot3d # Axes3D esta dentro de aqui
from matplotlib import gridspec
from os import path, kill, getppid, _exit
import signal

#global n_f, nfcompleto
#global ptos,tramos,pto_base,Gvec,A_secc,eps_L,Iy,Iz,Iyz,det,esferico
#global thetaPR,psi,eta,Wpsi,Weta,ppsi,peta,E,J_secc,I_a,cerrado,theta_secc
#global uT_medio,nucleo,C_tensor
#global n_ptos,n_tramos,yzsonPR,Ipsi,Ieta,v_psi,v_eta

ptos={}
tramos={}
n_f=''
elige=''

# definicion de clases para el calculo

class Pto:
    def __init__(self, y=0.,z=0., yG=0.,zG=0., inciP=[],
        uVy=0., uVz=0., uT=0., sxx=0.):

        # y,z= coordenadas en ejes iniciales, yG,zG en ejes por G
        # inciP[]= codigos de tramos conectados al pto
        # uVy, uVz, uT= despls ux para Vy=1, Vz=1, T=1
        # sxx= tension normal para los Mz, My, N dados

        self.y=y
        self.z=z
        self.yG=yG
        self.zG=zG
        self.inciP=inciP
        self.uVy=uVy
        self.uVz=uVz
        self.uT=uT
        self.sxx=sxx
        
    def pon_yGzG(self):
        # cambia coordenadas a ejes por G (requiere saberlo)
        global Gvec
        self.yG= self.y - Gvec[0]
        self.zG= self.z - Gvec[1]



class  Tramo:
    def __init__(self, tipoT=-1, inciT=[0,0], e=0.,   
                        # esto eran datos
        L=0.,alfa=0.,signo=1.,alfai=0.,R=0.,yC=0.,zC=0.,GT=[0.,0.],
                        # geometria, alguna dato
        Qy=0., Qz=0., iQy=0., iQz=0., Asect=0., 
                        # propiedades estaticas (las I_tr decidi no guardarlas)
        qVy=[0.,0.],qVz=[0.,0.],qT=[0.,0.], 
                        # flujos 'unitarios' en ptoi & ptoj
        tau_bi=0.,      # es la tension bidireccional de torsion
        tics=[], uVy_plot=[], uVz_plot=[], uT_plot=[],
        qVy_plot=[], qVz_plot=[]):
        
        # tipoT= 0 si es recto, =1 si es circular 
        # inciT= [ptoi, ptoj] del tramo 
        # e, L= espesor & longitud (no es dato)
        # alfa= angulo con eje y si recto, arco abarcado si circular
        # signo= signo de alfa
        # alfai= angulo de ptoi con el eje y (solo circular)
        # R, yC, zC= radio & coord del centro (solo circ)
        # GT[yGT,zGT]= centro de Gravedad del Tramo
        # Qy, Qz, iQy, iQz= mtos estaticos y sus integrales (ejes por G)
        # Asect= area sectorial del elemento (polo E, que se requiere)
        # qVy[0.,0.], qVz[0.,0.], qT[0.,0.]= flujos unitarios en i & j
        # tau_bi= maxima tension bidireccional (de ida y vuelta), torsion
        # tics =[y,z,s] sobre el elemento, para el trazado
        # uVy_plot =[ux con Vy=1] para los tics, & uVz_plot analogo & uT_plot
        # qVy_plot=[], qVz_plot=[] lo mismo con los flujos V
        # Para torsion qT=cte en el elemento. No hace falta muestreo. 
        
        self.tipoT=tipoT
        self.inciT=inciT

        self.e, self.L, self.alfa, self.signo, self.alfai= e,L,alfa,signo,alfai
        self.R, self.yC, self.zC, self.GT= R, yC, zC, GT
        self.Qy, self.Qz, self.iQy, self.iQz, self.Asect= Qy,Qz,iQy,iQz,Asect
        self.qVy, self.qVz, self.qT, self.tau_bi= qVy, qVz, qT, tau_bi
        self.tics, self.uVy_plot, self.uVz_plot= tics, uVy_plot, uVz_plot
        self.uT_plot, self.qVy_plot, self.qVz_plot = uT_plot, qVy_plot, qVz_plot


    # funciones para momentos estaticos:
    
    def Qz_s(self,s):
        if self.tipoT ==0:
            i,j=self.inciT
            return self.e*s*(ptos[i].yG+s*np.cos(self.alfa)/2.)
        else:
            theta= np.sign(self.alfa)* s/self.R
            a= theta*self.yC+ self.R*(np.sin(self.alfai+theta)
                                - np.sin(self.alfai))
            return(a*self.R*self.e*np.sign(self.alfa))

    def iQz_s(self,s):
        if self.tipoT == 0:
            i,j=self.inciT
            return self.e*s*s*(ptos[i].yG/2+s*np.cos(self.alfa)/6)
        else:
            theta= np.sign(self.alfa)* s/self.R
            a= (-np.cos(self.alfai+theta)+ np.cos(self.alfai)- 
                theta*np.sin(self.alfai))*self.e*self.R**3
            a += self.yC*self.e*(theta*self.R)**2 /2.
            return(a)
        
    def Qy_s(self,s):
        if self.tipoT == 0:
            i,j=self.inciT
            return self.e*s*(ptos[i].zG+s*np.sin(self.alfa)/2.)
        else:
            theta= np.sign(self.alfa)* s/self.R
            a= theta*self.zC+ self.R*(-np.cos(self.alfai+theta)
                                    + np.cos(self.alfai))
            return(a*self.R*self.e*np.sign(self.alfa))     
        
    def iQy_s(self,s):
        if self.tipoT == 0:
            i,j=self.inciT
            return self.e*s*s*(ptos[i].zG/2+s*np.sin(self.alfa)/6)
        else:
            theta= np.sign(self.alfa)* s/self.R
            a  = (-np.sin(self.alfai+theta)+ np.sin(self.alfai)
                +theta*np.cos(self.alfai))*self.e*self.R**3
            a += self.zC*self.e*(theta*self.R)**2 /2.
            return(a)

    # calculo de flujo debido al cortante, qV(s)
    def qV_s(self, s, Vy, Vz):
        global Iy,Iz,Iyz,det
        qi= self.qVy[0]* Vy + self.qVz[0]*Vz
        return(-qi - ( (self.Qz_s(s)*Iy-self.Qy_s(s)*Iyz)*Vy + 
                       (self.Qy_s(s)*Iz-self.Qz_s(s)*Iyz)*Vz ) / det)

    # calculo del ux debido al cortante, uV(s)
    def uV_s(self, s, Vy,Vz):
        global Iy,Iz,Iyz,det
        i,j = self.inciT
        a = -(self.iQz_s(s)*Iy-self.iQy_s(s)*Iyz)*Vy
        a -= (self.iQy_s(s)*Iz-self.iQz_s(s)*Iyz)*Vz
        a /= det
        a -= (self.qVy[0]*Vy+self.qVz[0]*Vz)*s
        a /= self.e
        a += ptos[i].uVy*Vy + ptos[i].uVz*Vz
        return(a)


    # calculo del Mto respecto de G para Vy=1 & para Vz=1

    def MG(self):
    
        if self.tipoT ==0:
            esy, esz = 0., 0.  
            L=self.L
            i,j= self.inciT
            yi,yj, zi,zj= ptos[i].yG, ptos[j].yG, ptos[i].zG, ptos[j].zG
            esy, esz = (yj-yi)/L, (zj-zi)/L   # vec unitario i -> j
            a = yi*esz - zi*esy
            MG_Vy = a*(-self.qVy[0]*L- (self.iQz*Iy-self.iQy*Iyz)/det)
            MG_Vz = a*(-self.qVz[0]*L- (self.iQy*Iz-self.iQz*Iyz)/det)

        else:
            # Tramo circular. Resultado de symPy con cse() y collect()
            yC,zC,e,R = self.yC, self.zC, self.e, self.R
            alfai,alfa,signo = self.alfai, self.alfa, np.sign(self.alfa)
            # Iy,Iz,Iyz,det -> se toman del global 
            
            ### para Vy=1: ###
            qi=self.qVy[0]

            x0 = np.cos(alfai)
            x1 = qi*zC
            x2 = np.sin(alfai)
            x3 = qi*yC
            x4 = 1/det
            x5 = yC*zC
            x6 = zC**2
            x7 = R**3
            x8 = yC**2
            x9 = R**2
            x10 = x2**2
            x11 = x10/2
            x12 = x0**2
            x13 = x0*x2
            x14 = alfa + alfai
            x15 = np.cos(x14)
            x16 = np.sin(x14)
            x17 = R*qi
            x18 = x16**2
            x19 = alfa*x18
            x20 = x15**2
            x21 = alfa*x20
            x22 = alfa*x15
            x23 = alfa*x16
            x24 = x18/2
            x25 = alfa**2
            x26 = x20*x25/2

            c1=(
            # termino en Vy:
            (-2*Iy*R*e*signo*x0*x4*x8 + 2*Iy*R*e*signo*x15*x4*x8 + 
            2*Iy*R*e*signo*x16*x4*x5 - 2*Iy*R*e*signo*x2*x4*x5 - 
            2*Iy*R*e*signo*x22*x4*x5 + 2*Iy*R*e*signo*x23*x4*x8 + 
            2*Iy*e*signo*x0**3*x4*x7 + 2*Iy*e*signo*x0*x10*x4*x7 + 
            2*Iy*e*signo*x11*x4*x9*yC - Iy*e*signo*x13*x4*x9*zC - 
            2*Iy*e*signo*x15**3*x4*x7 - Iy*e*signo*x15*x16*x4*x9*zC - 
            2*Iy*e*signo*x15*x18*x4*x7 + 2*Iy*e*signo*x15*x2*x4*x9*zC - 
            2*Iy*e*signo*x16*x2*x4*x9*yC - 2*Iy*e*signo*x19*x2*x4*x7 + 
            Iy*e*signo*x19*x4*x9*zC - 2*Iy*e*signo*x2*x21*x4*x7 + 
            Iy*e*signo*x21*x4*x9*zC + 2*Iy*e*signo*x24*x25*x4*x9*yC + 
            2*Iy*e*signo*x24*x4*x9*yC + 2*Iy*e*signo*x26*x4*x9*yC + 
            2*Iyz*R*e*signo*x0*x4*x5 - 2*Iyz*R*e*signo*x15*x4*x5 - 
            2*Iyz*R*e*signo*x16*x4*x6 + 2*Iyz*R*e*signo*x2*x4*x6 + 
            2*Iyz*R*e*signo*x22*x4*x6 - 2*Iyz*R*e*signo*x23*x4*x5 + 
            2*Iyz*e*signo*x0*x15*x4*x9*zC - 2*Iyz*e*signo*x0*x16*x4*x9*yC - 
            2*Iyz*e*signo*x0*x19*x4*x7 - 2*Iyz*e*signo*x0*x21*x4*x7 - 
            2*Iyz*e*signo*x11*x4*x9*zC - 2*Iyz*e*signo*x12*x2*x4*x7 - 
            2*Iyz*e*signo*x12*x4*x9*zC + Iyz*e*signo*x13*x4*x9*yC + 
            Iyz*e*signo*x15*x16*x4*x9*yC + 2*Iyz*e*signo*x16**3*x4*x7 + 
            2*Iyz*e*signo*x16*x20*x4*x7 + Iyz*e*signo*x19*x4*x9*yC - 
            2*Iyz*e*signo*x2**3*x4*x7 + Iyz*e*signo*x21*x4*x9*yC - 
            2*Iyz*e*signo*x24*x25*x4*x9*zC + 2*Iyz*e*signo*x24*x4*x9*zC - 
            2*Iyz*e*signo*x26*x4*x9*zC) +
            # termino ni Vy ni Vz:
            2*x0*x1 - 2*x1*x15 + 2*x16*x3 + 2*x17*x19 + 2*x17*x21 - 2*x2*x3
             )
            
            MG_Vy = -R*c1 /2



            ### para Vz=1: ###
            qi=self.qVz[0]
            
            x0 = np.cos(alfai)
            x1 = qi*zC
            x2 = np.sin(alfai)
            x3 = qi*yC
            x4 = 1/det
            x5 = Iyz*yC
            x6 = x5*zC
            x7 = Iyz*yC**2
            x8 = Iz*yC
            x9 = x8*zC
            x10 = R**3
            x11 = Iz*zC**2
            x12 = x2**2
            x13 = R**2
            x14 = x0**2
            x15 = x14/2
            x16 = x0*x2/2
            x17 = alfa + alfai
            x18 = np.cos(x17)
            x19 = np.sin(x17)
            x20 = x19**2
            x21 = R*alfa
            x22 = qi*x21
            x23 = x18**2
            x24 = x23/2
            x25 = x20/2
            x26 = alfa*x2
            x27 = alfa*x0
            x28 = alfa**2
            x29 = x25*x28
            
            c1=(
            # termino en Vz:
            (2*Iyz*alfa*e*signo*x13*x24*x4*zC + 
            2*Iyz*alfa*e*signo*x13*x25*x4*zC + 2*Iyz*e*signo*x0**3*x10*x4 + 
            2*Iyz*e*signo*x0*x10*x12*x4 - 2*Iyz*e*signo*x10*x18**3*x4 - 
            2*Iyz*e*signo*x10*x18*x20*x4 - 2*Iyz*e*signo*x10*x20*x26*x4 - 
            2*Iyz*e*signo*x10*x23*x26*x4 - 2*Iyz*e*signo*x13*x16*x4*zC - 
            Iyz*e*signo*x13*x18*x19*x4*zC + 2*Iyz*e*signo*x13*x18*x2*x4*zC + 
            2*Iz*e*signo*x0*x13*x18*x4*zC - 2*Iz*e*signo*x10*x14*x2*x4 + 
            2*Iz*e*signo*x10*x19**3*x4 + 2*Iz*e*signo*x10*x19*x23*x4 - 
            2*Iz*e*signo*x10*x2**3*x4 - 2*Iz*e*signo*x10*x20*x27*x4 - 
            2*Iz*e*signo*x10*x23*x27*x4 - 2*Iz*e*signo*x13*x15*x4*zC - 
            2*Iz*e*signo*x13*x24*x28*x4*zC - 2*Iz*e*signo*x13*x24*x4*zC - 
            2*Iz*e*signo*x13*x29*x4*zC - 2*R*e*signo*x0*x4*x7 + 
            2*R*e*signo*x0*x4*x9 - 2*R*e*signo*x11*x19*x4 + 
            2*R*e*signo*x11*x2*x4 + 2*R*e*signo*x18*x4*x7 - 
            2*R*e*signo*x18*x4*x9 + 2*R*e*signo*x19*x4*x6 - 
            2*R*e*signo*x2*x4*x6 + 2*alfa*e*signo*x13*x24*x4*x8 + 
            2*alfa*e*signo*x13*x25*x4*x8 - 2*e*signo*x0*x13*x19*x4*x8 + 
            2*e*signo*x11*x18*x21*x4 + 2*e*signo*x12*x13*x4*x5 + 
            2*e*signo*x13*x15*x4*x5 + 2*e*signo*x13*x16*x4*x8 + 
            e*signo*x13*x18*x19*x4*x8 - 2*e*signo*x13*x19*x2*x4*x5 + 
            2*e*signo*x13*x24*x28*x4*x5 - 2*e*signo*x13*x24*x4*x5 + 
            2*e*signo*x13*x29*x4*x5 - 2*e*signo*x18*x21*x4*x6 + 
            2*e*signo*x19*x21*x4*x7 - 2*e*signo*x19*x21*x4*x9)
            # termino ni Vy ni Vz: 
            - 2*x0*x1 + 2*x1*x18 - 2*x19*x3 + 2*x2*x3 - 2*x20*x22 - 2*x22*x23
              )

            MG_Vz= R*c1/2

        return (MG_Vy, MG_Vz)


    # calculo de la integral (u(s)-u_med )**2 para el I_a
    
    def int2(self):
        i,j = self.inciT
        L= self.L
        yi,zi = ptos[i].yG, ptos[i].zG
        yE, zE = E

        if self.tipoT == 0:
            alfa=self.alfa
            a= (yi-yE)*np.sin(alfa)-(zi-zE)*np.cos(alfa)
            c1 = ptos[i].uT - uT_medio
            c2 = self.qT[1]/self.e - theta_secc * a
            return( L*(c1**2+ (c2*L)**2 /3+ c1*c2*L) )
        else:
            yC,zC,e,R = self.yC, self.zC, self.e, self.R
            alfai,alfa,signo = self.alfai, self.alfa, np.sign(self.alfa)
            
            c1 = ptos[i].uT - uT_medio
            q = self.qT[1]
            tt=theta_secc

            x0 = np.cos(alfai)
            x1 = 2*x0
            x2 = c1*tt
            x3 = R*x2
            x4 = x3*yC
            x5 = x3*yE
            x6 = np.sin(alfai)
            x7 = 2*x6
            x8 = x3*zC
            x9 = x3*zE
            x10 = alfa*x7
            x11 = alfa*x1
            x12 = R**2
            x13 = alfa**2
            x14 = x12*x13
            x15 = alfa + alfai
            x16 = np.cos(x15)
            x17 = 2*x16
            x18 = np.sin(x15)
            x19 = 2*x18
            x20 = tt**2
            x21 = alfa**3
            x22 = x21/3
            x23 = R**3
            x24 = x20*yC
            x25 = x23*x24
            x26 = x20*yE
            x27 = x23*x26
            x28 = x20*x23
            x29 = x28*zC
            x30 = x28*zE
            x31 = q*signo/e
            x32 = x12*x6**2
            x33 = x24*zE
            x34 = x26*zC
            x35 = x13*x6
            x36 = x0*x13
            x37 = x24*zC
            x38 = x0**2
            x39 = 2*x38
            x40 = x12*x39
            x41 = x12*x33
            x42 = x26*zE
            x43 = alfa*x17
            x44 = alfa*x19
            x45 = x24*yE
            x46 = x12*x45
            x47 = x0*x6
            x48 = 3*x47
            x49 = x12*x20
            x50 = zC*zE
            x51 = x49*x50
            x52 = tt*x31
            x53 = x12*x52
            x54 = x53*yC
            x55 = x53*yE
            x56 = x53*zC
            x57 = x53*zE
            x58 = yC**2
            x59 = alfa*x20
            x60 = x32*x59
            x61 = yE**2
            x62 = zC**2
            x63 = x12*x59
            x64 = x38*x63
            x65 = zE**2
            x66 = x18**2
            x67 = x12*x66
            x68 = x49*x58
            x69 = 3*x47/2
            x70 = x49*x61
            x71 = x49*x62
            x72 = x49*x65
            x73 = alfa*x45
            x74 = x12*x37
            x75 = x0*x10
            x76 = x12*x34
            x77 = x12*x42
            x78 = x18*x7
            x79 = x0*x17
            x80 = x14*x52
            x81 = x6*x80
            x82 = x0*x80
            x83 = x16*x7
            x84 = x0*x19
            x85 = x16**2
            x86 = x50*x63
            x87 = x16*x18
            x88 = x63/2
            x89 = x58*x88
            x90 = x61*x88
            x91 = x62*x88
            x92 = x65*x88
            x93 = x87/2
        
            a=(
            R**4*x20*x22 + R*c1*x13*x31 + alfa*c1**2 - alfa*x46*x85 + 
            4*x0*x18*x51 + x1*x29 - x1*x30 - x1*x4 + x1*x5 - x1*x56 + 
            x1*x57 + x10*x4 - x10*x5 - x11*x8 + x11*x9 - x14*x2 - 
            4*x16*x46*x6 - x17*x29 + x17*x30 + x17*x4 - x17*x5 + x17*x56 - 
            x17*x57 + x19*x25 - x19*x27 - x19*x54 + x19*x55 + x19*x8 - 
            x19*x9 - 2*x21*x23*x52/3 - x25*x35 - x25*x43 - x25*x7 + 
            x27*x35 + x27*x43 + x27*x7 + x29*x36 - x29*x44 - x30*x36 + 
            x30*x44 + x32*x33 + x32*x34 - x32*x37 - x32*x42 - 2*x32*x73 + 
            x33*x67 - x34*x40 + x34*x67 + x37*x40 - x37*x67 - x39*x41 + 
            x40*x42 + x41*x75 - x41*x78 + x41*x79 - x42*x67 + x43*x54 - 
            x43*x55 + x44*x56 - x44*x57 + x46*x48 + x46*x87 - x48*x51 - 
            2*x50*x64 - x51*x87 + x54*x7 - x55*x7 + x58*x60 + x60*x61 + 
            x62*x64 + x64*x65 - x66*x86 + x66*x89 + x66*x90 + x66*x91 + 
            x66*x92 - x67*x73 - x68*x69 + x68*x83 - x68*x93 - x69*x70 + 
            x69*x71 + x69*x72 - x7*x8 + x7*x9 + x70*x83 - x70*x93 - 
            x71*x84 + x71*x93 - x72*x84 + x72*x93 - x74*x75 + x74*x78 - 
            x74*x79 + x75*x76 - x75*x77 - x76*x78 + x76*x79 + x77*x78 - 
            x77*x79 + x81*yC - x81*yE - x82*zC + x82*zE - x85*x86 + 
            x85*x89 + x85*x90 + x85*x91 + x85*x92 + 
            q**2*signo**2*x12*x22/e**2   )
            
            a *= R*signo
            return(a)
        
        

    # calculo del area sectorial (s). OJO *2 porque asi sale siempre
    
    def Asect_s (self, s):
        yE,zE=E
        if self.tipoT == 0:  # tramo recto
            i,j = self.inciT
            yi,zi = ptos[i].yG, ptos[i].zG
            alfa=self.alfa
            a= (yi-yE)*np.sin(alfa)-(zi-zE)*np.cos(alfa)
            a *= s
            return(a)
        else:
            alfa, alfai = self.alfa, self.alfai
            R, yC,zC = self.R, self.yC, self.zC
            ECy, ECz = (yC-yE), (zC-zE)
            theta = s*np.sign(alfa)/R
            a=( ECy* (np.sin(alfai+theta)-np.sin(alfai)) +
                ECz* (-np.cos(alfai+theta)+np.cos(alfai)) + R* theta )
            a *= R
            return(a)


    # calculo del desplazamiento ux de torsion (s) en el tramo

    def uT_s(self,s):
        i,j = self.inciT
        a  = - self.Asect_s(s) * theta_secc
        a += ptos[i].uT + self.qT[1]*s/self.e
        return (a)
    
    # calculo de [y,z] (coord iniciales) dado s
    def yz_s(self,s):
        i,j=self.inciT
        if self.tipoT ==0:
            y= ptos[i].y + s*np.cos(self.alfa)
            z= ptos[i].z + s*np.sin(self.alfa)
        else:
            theta= s*self.signo/self.R
            y= Gvec[0]+ self.yC+ self.R*np.cos(self.alfai+theta)
            z= Gvec[1]+ self.zC+ self.R*np.sin(self.alfai+theta)
            # OJO si se llama desde fuera de motor_calculo, yC,zC estaran
            # en ejes dados, y sumar Gvec sobra (hay que restarlo luego)
        return(np.array([y,z]))
    
    # calculo de vector tangente es[], y normal a la dcha en[], dado s
    def es_s(self,s):
        theta=s*self.signo/self.R
        es=np.array([-np.sin(self.alfai+theta), np.cos(self.alfai+theta)])
        es *= self.signo
        return(es)
    def en_s(self,s):
        theta=s*self.signo/self.R
        en=np.array([np.cos(self.alfai+theta), np.sin(self.alfai+theta)])
        en *= self.signo
        return(en)


class Recta:
    def __init__(self, p=[0.,0.], es=[0.,0.], en=[0.,0.]):
    # para guardar un punto p de la recta y vector colineal y normal 
        self.p=np.array(p)
        self.es=np.array(es)
        self.en=np.array(en)
        
    def intersec(self,r2):
        # calcula la posicion s en cada recta del pto de interseccion
        
        if abs( self.es[0]*r2.es[0]+self.es[1]*r2.es[1] -1.) < 1.e-6:
            texto='// Error en intersec: las lineas no intersectan.\n'
            print(texto)
            quit() 
        
        det= -self.es[0]*r2.es[1]+ self.es[1]*r2.es[0]
        if abs(det)<1.0e-9: print('W: bad condition det in intersec.')
        
        incy= r2.p[0]- self.p[0]
        incz= r2.p[1]- self.p[1]
        s1= (-incy*r2.es[1] + incz*r2.es[0] ) /det
        s2= ( incz*self.es[0] - incy*self.es[1] ) /det
        return (s1, s2)

    def dist (self,Q):
        # calcula la distancia de la recta al punto Q. 
        r=[ Q[0]-self.p[0] , Q[1]-self.p[1] ]
        d= float( np.cross(r, self.es) )
        return d    # d positivo si esta a la dcha de la recta segun es[]



#       ##########################################
#       ####   FIN DE DEFINICION DE OBJETOS   ####
#       #### COMIENZA DEFINICION DE FUNCIONES ####
#       ##########################################

    
# trasladar filas del GUI a ptos

def filas_a_ptos():
    global ptos,tramos, n_ptos, n_tramos, filas_ptos, pto_base
    ptos={}
    for fila in filas_ptos:
        p=Pto()
        x0,x1,x2= fila[0].get(), fila[1].get(), fila[2].get()
        if x0 != '':
            try:
                ip=int(x0)
                y=float(x1)
                z=float(x2)
            except (TypeError, ValueError):
                texto= 'Hay un error en los puntos[]. Compruebe'
                texto+=' que están todos los datos, etc.'
                texto+='\n             '+ str(sys.exc_info()[0])
                messagebox.showinfo(message=texto,title='MAL',parent=v9)
                return(1)
            p.y, p.z = y,z
            ptos[ip] = p
                # y mejor reasignar n_ptos a cada paso que perderlo si la
                # lista tiene todas las hard_nfilas con datos:
            n_ptos=len(ptos) 
            try:
                kk=pto_base
            except (NameError,TypeError, ValueError):
                pto_base=ip   # queda el primero o bien el que manda el usuario
        
        else:   # campo vacio => se termino la lista
            return(0)


# trasladar filas del GUI a tramos

def filas_a_tr():
    global ptos, tramos, n_ptos, n_tramos, filas_tramos
    tramos={}
    for fila in filas_tramos:
        tr=Tramo()
        x=[]
        for i in range(6): x.append(fila[i].get())
        if x[0] != '':
            try:
                it=int(x[0])
                ini, fin, e = int(x[1]), int(x[2]), float(x[3])
            except (TypeError, ValueError):
                texto= 'Hay un error en los tramos[]. Por favor compruebe '
                texto+='si debe o no haber puntos decimales, etc'
                texto+='\n             '+ str(sys.exc_info()[0])
                messagebox.showinfo(message=texto,title='MAL',parent=v9)
                return(1)
            
            tipoT=0         # si hay R valido se reasigna =1 mas abajo
            if x[4] != '':
                try:
                    R = float(x[4])
                except(TypeError, ValueError):
                    texto= 'Hay un radio R en los tramos[] que no se entiende.'
                    texto+='Compruebe los puntos decimales, espacios, etc.'
                    texto+='\n             '+ str(sys.exc_info()[0])
                    messagebox.showinfo(message=texto,title='MAL',parent=v9)
                    return(1)
                signo= -1.0 if x[5]=='-' else 1.0
                tipoT=1

            tr.tipoT=tipoT
            tr.inciT=[ini,fin]
            tr.e=e
            if tipoT ==1:
                tr.R=R
                tr.signo=signo

            tramos[it] = tr
            n_tramos=len(tramos) 
            
        else:   # campo vacio => se termino la lista
            return(0)


# Ventana de opciones 

def entraOpc():
    global uV_despl, uT_despl, pto_base, ptos,etiq2D,etiq3D,peque3D
    
    def hechoOpc():
        global uV_despl, uT_despl, pto_base, ptos,etiq2D,etiq3D,peque3D
        uV_despl= int(uV_di.get())
        uT_despl= int(uT_di.get())
        pto_base= int(pto_di.get())
        etiq2D=int(etiq2D_di.get())
        etiq3D=int(etiq3D_di.get())
        peque3D=int(peque3D_di.get())

        v4.destroy()
        v9.focus_force()

    v4=Toplevel(v9)
    v4.title('Opciones')

    uV_di= StringVar()
    uT_di= StringVar()
    pto_di= StringVar()
    etiq2D_di= StringVar()
    etiq3D_di= StringVar()
    peque3D_di=StringVar()
    
    try:
        a,b= uV_despl,uT_despl
        kk1,kk2= int(a),int(b)
        uV_di.set(a)
        uT_di.set(b)
    except (ValueError,NameError):
        uV_di.set('1')
        uT_di.set('2')
        
    try:
        kk1=int(pto_base)
    except (ValueError,NameError):
        kk1=1
    pto_di.set(str(kk1))
    
    try:
        kk1, kk2= int(etiq2D), int(etiq3D)
        etiq2D_di.set(kk1)
        etiq3D_di.set(kk2)
    except (ValueError,NameError):
        etiq2D_di.set(1)
        etiq3D_di.set(1)
    
    try:
        kk1= int(peque3D)
        peque3D_di.set(kk1)
    except (ValueError,NameError):
        peque3D_di.set(0)

    ttk.Separator(v4, orient=VERTICAL).grid(row=0, column=0, sticky='nse')

    a= '\nEs posible elegir qué punto se tomará como inmóvil  \n'
    a+='(ux=0) en los cálculos. La elección es arbitraria.\n' 
    a+='Esta opción requiere recalcular la sección (hágalo\n'
    a+='si la usa). Las demás opciones no lo requieren,\n'
    a+='considere si pueden servirle mejor.'
    ttk.Label(v4, text=a, foreground='green').grid(row=0,
        column=2, sticky='w')

    e1=ttk.Entry(v4, textvariable=pto_di, width=5)
    e1.grid(row=1,column=0, columnspan=2)
    
    ttk.Label(v4, text='Punto de movimiento ux nulo').grid(row=1,
        column=2, sticky='w')

    ttk.Separator(v4, orient=HORIZONTAL).grid(row=2,
        column=2,sticky='ew')

    a = 'Las gráficas 3D de desplazamientos de alabeo\n'
    a += 'pueden mostrarse desplazadas para su mejor\n'
    a += 'visualización. Las opciones son:'
    ttk.Label(v4, text=a, foreground='green').grid(row=3,
        column=2, sticky='w')

    ttk.Label(v4, text='ux_V').grid(row=4,column=0)
    ttk.Label(v4, text='ux_T').grid(row=4,column=1)
    
    
    ttk.Radiobutton(v4, variable=uV_di, value='0').grid(row=5,
        column=0)
    ttk.Radiobutton(v4, variable=uT_di, value='0').grid(row=5,
        column=1)
    ttk.Label(v4, text='Dibujar tal como se calcula').grid(row=5,
        column=2, sticky='w')
    
    
    ttk.Radiobutton(v4, variable=uV_di, value='1').grid(row=6,
        column=0)
    ttk.Radiobutton(v4, variable=uT_di, value='1').grid(row=6,
        column=1)
    ttk.Label(v4, text='Desplazar afuera de la sección').grid(row=6,
        column=2, sticky='w')

    ttk.Radiobutton(v4, variable=uV_di, value='2').grid(row=7,
        column=0)
    ttk.Radiobutton(v4, variable=uT_di, value='2').grid(row=7,
        column=1)
    ttk.Label(v4,text='Centrar al promedio el trazado').grid(row=7,
        column=2, sticky='w')
    
    ttk.Separator(v4, orient=HORIZONTAL).grid(row=8,
        column=2, sticky='ew')
        
    ttk.Separator(v4, orient=VERTICAL).grid(row=3,
        column=0, sticky='nse')
    
    a= 'Las etiquetas numéricas de las gráficas a veces\n'
    a+='molestan. Puede deshabilitar su trazado aquí.\n'
    a+='Estas opciones no afectan a la ventana de\n'
    a+='comprobación de datos ni a la ficha de la sección.'
    ttk.Label(v4, text=a, foreground='green').grid(row=9,
        column=2, sticky='w')
    ttk.Separator(v4, orient=VERTICAL).grid(row=9,
        column=0, sticky='nse')
    
    ttk.Label(v4,text='Etiquetar gráficas 2D').grid(row=10,
        column=2, sticky='w')
    ttk.Checkbutton(v4, variable=etiq2D_di).grid(
        row=10, column=0, sticky='ns', columnspan=2)
    
    ttk.Label(v4,text='Etiquetar gráficas 3D').grid(
        row=11, column=2, sticky='w')
    ttk.Checkbutton(v4, variable=etiq3D_di).grid(
        row=11, column=0, sticky='ns', columnspan=2)
    
    ttk.Separator(v4, orient=HORIZONTAL).grid(row=12,
        column=2, sticky='ew')
    ttk.Separator(v4, orient=VERTICAL).grid(row=13,
        column=0, sticky='nse')
        
    a= 'Los detalles de las graficas 3D de alabeo se\n'
    a+='aprecian mejor dibujadas a escala grande, pero\n'
    a+='resultan más intuitivas a una escala pequeña.'
    ttk.Label(v4, text=a, foreground='green').grid(row=13,
        column=2, sticky='w')
    
    ttk.Label(v4,text='Moderar los trazados de alabeo').grid(row=14,
        column=2, sticky='w')
    ttk.Checkbutton(v4, variable=peque3D_di).grid(
        row=14, column=0, sticky='ns', columnspan=2)
    
    
    ttk.Button(v4, text='Hecho', command=hechoOpc).grid(
        row=15, column=2)

    for hijo in v4.winfo_children(): hijo.grid_configure(padx=6,
        pady=6)

    v4.geometry('+410-70')
    v4.focus()
    e1.focus()
    v4.mainloop()




# GUI -  primera  ventana (de presentacion y elecciones) 

def presenta_elige():
    global elige,n_f, nfcompleto


    ##### Ventana hija para licencia #####

    def licencia():
        v1=Toplevel(v0)

        texto='Este programa es Software Libre (Free Software). Como '
        texto+='tal se le aplican los términos de la "GNU General Public '
        texto+='License", en su versión 2 o bien (como usted prefiera) '
        texto+='en una versión posterior. Básicamente, usted puede: \n\n'
        texto+='* Usar libremente el programa \n'
        texto+='* Realizar copias del mismo y distribuirlas libremente \n'
        texto+='* Estudiar su código para aprender cómo funciona \n'
        texto+='* Realizar modificaciones/mejoras del programa \n\n'
        texto+='Bajo las siguientes condiciones:  \n\n'
        texto+='* Las modificaciones realizadas al programa deben hacerse '
        texto+='públicas bajo una licencia como la presente (así el '
        texto+='software puede mejorar con aportaciones realizadas '
        texto+='sobre el trabajo de otros, como se hace en la ciencia).\n'
        texto+='* Las modificiones y trabajos derivados en general deben '
        texto+='incluir el reconocimiento al autor original (no puede '
        texto+='decir que es usted quien ha escrito este programa).\n'
        texto+='* En este caso, debe mencionar al autor original como: \n'
        texto+='Juan Carlos del Caño, profesor en la Escuela de \n'
        texto+='Ingenierías Industriales de la Universidad de  \n'
        texto+='Valladolid (Spain)  \n\n'
        texto+='Este programa se distribuye con la esperanza de que sea '
        texto+='útil, pero SIN NINGUNA GARANTIA, ni siquiera la garantía '
        texto+='de comerciabilidad o de adecuación para un propósito '
        texto+='particular. Lea la GNU General Public License para más '
        texto+='detalles.\n\n'
        texto+='Usted debería haber recibido una copia de la GNU General '
        texto+='Public License junto con este programa. Si no es así,  '
        texto+='escriba a la Free Software Foundation: \n'
        texto+='Inc., 51 Franklin Street, Fifth Floor, \n'
        texto+='Boston, MA 02110-1301, USA.'

        ttk.Label(v1, text='Resumen de términos de licencia\n', 
            font=('', 16)).grid(row=0, column=0, columnspan=2, pady=5)

        tcaja = Text(v1, width=45, height=30,wrap='word', font=('Sans',9),
            background='#D9D9D9', foreground='green', border=None, padx=20, pady=12)
        tcaja.grid(column=0, row=1, padx=8, sticky=(N,W,E,S))
        tcaja.insert('1.0',texto)

        scb = ttk.Scrollbar(v1, orient=VERTICAL, command=tcaja.yview)
        scb.grid(column=1, row=1, sticky='ns')
        
        tcaja['yscrollcommand'] = scb.set

        ttk.Button(v1, text='Entendido', width=9, command=v1.destroy).grid(
            column=0, row=2, pady=4, columnspan=2)

        tcaja['state']='disabled'

        v1.grid_columnconfigure(0, weight=1)
        v1.grid_rowconfigure(0, weight=1)
        v1.grid_rowconfigure(1, weight=4)
        v1.grid_rowconfigure(2, weight=1)
        v1.geometry('+240+60')

        v1.focus()
        v1.mainloop()


    ##### Ventana hija para breviario #####

    def breviario():

        v2=Toplevel(v0)
            
        texto='- Para qué sirve este programa:\n\n'
        texto+='Analiza secciones de pared delgada de barras resistentes '
        texto+='(vigas etc), típicamente metálicas, en régimen lineal. A '
        texto+='partir de los datos geométricos de la sección calcula el '
        texto+='centro de áreas, los ejes principales de inercia, las inercias y '
        texto+='módulos resistentes, el núcleo central, la constante de '
        texto+='torsión (J), el módulo de alabeo (Ia), el centro de '
        texto+='esfuerzos cortantes (E), el tensor de giros a cortante '
        texto+='(ver documentación), y las áreas a cortante.\n'
        texto+='Si se proporcionan valores para algunos o todos los '
        texto+='esfuerzos (cortante, axil, flector, torsor), se dibujan '
        texto+='además varias gráficas de alto valor ilustrativo, que  '
        texto+='incluyen gráficas 3D practicables del alabeo de la sección.'
        texto+='\nLa salida de texto ofrece información adicional.'
        
        texto+='\n\n- Cómo lo hace:\n\n'
        texto+='Se asume de forma consistente el modelo de pared delgada, '
        texto+='tanto en lo relativo a geometría como a tensiones. El usuario '
        texto+='debe juzgar en primer lugar si dicho modelo es adecuado para '
        texto+='analizar su problema.\n'
        texto+='La geometría de una sección de pared delgada consta prácticamente '
        texto+='siempre de tramos rectos y/o tramos circulares. Estos últimos '
        texto+='corresponden con frecuencia a radios de acuerdo, aunque también '
        texto+='pueden existir arcos de circunferencia \'grandes\' en la sección '
        texto+='(secciones tipo tubo modificadas, etc).\n'
        texto+='El programa proporciona tramos rectos y tramos circulares, con lo que '
        texto+='la geometría de la sección será modelada exactamente en lo que '
        texto+='al modelo de pared delgada se refiere. Además, el cálculo de '
        texto+='desplazamientos y flujos de cortante y/o de torsión se realiza '
        texto+='también de forma \'exacta\' utilizando expresiones analíticas. Los '
        texto+='cálculos derivados como la constante de torsión, el módulo de '
        texto+='alabeo, y el centro de esfuerzos cortantes se calculan también de '
        texto+='forma exacta usando expresiones analíticas (gracias al módulo SymPy '
        texto+='de cálculo simbólico que fue utilizado para su obtención). La '
        texto+='presencia o no de bucles cerrados en la sección es indiferente.\n'

        texto+='\n- Pautas para un uso correcto:\n\n'
        texto+=' * La geometría de la sección se modela con TRAMOS rectos '
        texto+='y circulares de espesor constante que unen sucesivos '
        texto+='PUNTOS de la linea media del perfil. Se deben numerar '
        texto+='los puntos y proporcionar sus coordenadas (y,z). '
        texto+='Se deben numerar los tramos y proporcionar la  '
        texto+='numeración de sus puntos extremos y su espesor. En el '
        texto+='caso de tramo circular, se proporcionará también su '
        texto+='radio y el signo del ángulo (positivo antihotario) al '
        texto+='recorrer el arco desde su primer punto hasta el segundo.\n'
        texto+=' * No tiene sentido \'discretizar\' los tramos rectos o '
        texto+='circulares de la sección en subtramos, ya que  un '
        texto+='único tramo proporciona ya la solución exacta. En los '
        texto+='raros casos de tramos de espesor continuamente variable '
        texto+='sí estaría indicada la aproximación por varios tramos, '
        texto+='ya que el \'tramo\' del programa es de espesor constante.\n'
        texto+=' * La numeración de puntos y tramos puede hacerse libremente '
        texto+='usando números enteros positivos. Dicha numeración puede '
        texto+='no ser correlativa ni compacta (puede no incluir a todos  '
        texto+='los números menores que el mayor proporcionado).\n'

        texto+=' * El programa entiende siempre que el arco entre dos  '
        texto+='puntos se recorre por el camino angular mas corto. '
        texto+='Esto simplifica la entrada de datos, si bien implica '
        texto+='que los tramos circulares serán menores de 180º. Si su problema '
        texto+='tiene un arco de 180º o mayor, simplemente utilice dos '
        texto+='tramos circulares para modelarlo.\n'
        texto+=' * Navegación: <Tab> & <Shift><Tab> para desplazarse '
        texto+='entre items de una ventana. <Space> pulsa un botón\n'
        
        texto+='\n- Estado de desarrollo del programa:\n\n'
        texto+='Esta versión 1.0 del programa tiene toda la capacidad de\n'
        texto+='cálculo de la versión 0.8 (versión estable anterior), mas '
        texto+='la implementación de tramos circulares . La interfaz de usuario '
        texto+='es más flexible y permite modificar un problema existente tras '
        texto+='cargarlo, asi como realizar análisis de sucesivos casos de carga.\n'
        texto+='Aunque no es probable que llegue a necesitarse, el fichero de datos '
        texto+='del problema sigue siendo accesible y fácilmente editable.\n'
        texto+='Espero que ThinSecBeam le sea útil.'
        
        ttk.Label(v2, text='Notas breves.\n', 
            font=('', 16)).grid(row=0, column=0, columnspan=2, pady=5)

        tcaja = Text(v2, width=45, height=30,wrap='word', font=('Sans',9),
            background='#D9D9D9', foreground='green', border=None, padx=20, pady=12)
        tcaja.grid(column=0, row=1, padx=8, sticky=(N,W,E,S))
        tcaja.insert('1.0',texto)

        scb = ttk.Scrollbar(v2, orient=VERTICAL, command=tcaja.yview)
        scb.grid(column=1, row=1, sticky='ns')
        
        tcaja['yscrollcommand'] = scb.set

        ttk.Button(v2, text='Entendido', width=9, command=v2.destroy).grid(
            column=0, row=2, pady=4, columnspan=2)

        tcaja['state']='disabled'

        v2.grid_columnconfigure(0, weight=1)
        v2.grid_rowconfigure(0, weight=1)
        v2.grid_rowconfigure(1, weight=4)
        v2.grid_rowconfigure(2, weight=1)
        v2.geometry('+250+70')

        v2.focus()
        v2.mainloop()



    ##### Para el prb por defecto #####

    def default_prb():
        global elige
        elige='default' # sin mas ventanas
        v0.destroy()



    ##### La ventana de inicio (por fin) #####

    v0 = Tk()
    v0.title("ThinSecBeam 1.0")

    cuadro = ttk.Frame(v0, padding='9 3 3 3') 
    cuadro.grid(column=0, row=0, sticky=(N, W, E, S))

    ttk.Label(cuadro, text='ThinSecBeam', font=('', 44)).grid(row=0,
        column=0, columnspan=4)
    ttk.Label(cuadro, text='Thin-walled Sections of Beams', 
        font=('Courier', 16)).grid(row=1, column=0, columnspan=4)
    ttk.Label(cuadro, text='by:   Juan Carlos del Caño\n').grid(row=2,
        column=0, columnspan=4)

    # hago la parte izda de la ventana

    ttk.Separator(cuadro, orient=HORIZONTAL).grid(column=0, row=3,
        columnspan=4, sticky='ew')
        
    texto= 'Esto es Software Libre (Free Software), lo cual\n'
    texto +='le otorga a usted algunos derechos pero también\n'
    texto +='conlleva algunas obligaciones. \nPor favor lea la licencia.'
    ttk.Label(cuadro, text=texto, foreground='green').grid(row=4, 
        column=0, columnspan=3, sticky='w')

    ttk.Separator(cuadro, orient=HORIZONTAL).grid(row=5, column=0,
        columnspan=3, sticky='ew')
        
    texto=  'Si usted no conoce para qué sirve este programa,\n'
    texto +='o tiene una duda básica respecto del mismo,\n'
    texto +='consulte estas breves notas.\n'
    texto +='O bien consulte la documentación.'
    ttk.Label(cuadro, text=texto, foreground='green').grid(row=6, 
        column=0, columnspan=3, sticky='w')

    ttk.Separator(cuadro, orient=HORIZONTAL).grid(row=7, column=0,
        columnspan=3, sticky='ew')
        
    texto = 'Puede empezar un problema nuevo desde cero.\n'
    texto +='Deberá especificar la geometría de la sección\n'
    texto +='(puntos y tramos) y algunos otros datos generales\n'
    texto +='del problema. También puede cargar desde archivo\n'
    texto +='un problema que haya sido guardado previamente.'
    ttk.Label(cuadro, text=texto, foreground='green').grid(row=8,
        column=0, columnspan=3, sticky='w')

    ttk.Separator(cuadro,orient=HORIZONTAL).grid(row=9,column=0,sticky='ew')
    ttk.Separator(cuadro,orient=HORIZONTAL).grid(row=9,column=2,sticky='ew')
    ttk.Label(cuadro, text='o bien:').grid(row=9,column=1)
        
    texto = 'Puede cargar un ejemplo preprogramado internamente\n'
    texto +='que ilustra las capacidades básicas del programa,\n'
    texto +='y con el que puede comprobarse si la instalación de\n'
    texto +='ThinSecBeam fue correcta.'

    ttk.Label(cuadro, text=texto, foreground='green').grid(row=10,
        column=0, columnspan=3, sticky='w')
        
    ttk.Separator(cuadro, orient=HORIZONTAL).grid(row=11, column=0,
        columnspan=4, sticky='ew')


    # ahora hago la parte derecha

    ttk.Separator(cuadro,orient=VERTICAL).grid(row=3,column=3,
        rowspan=9, sticky='ns')
    ttk.Button(cuadro, text='Licencia',
        command=licencia).grid(row=4, column=3)

    ttk.Button(cuadro, text='Breviario', command=breviario).grid(
        row=6, column=3)

    ttk.Button(cuadro, text='Resolver', command=v0.destroy).grid(
        row=8, column=3)

    ttk.Button(cuadro, text='Ejemplo', command=default_prb).grid(
        row=10, column=3)

    for hijo in cuadro.winfo_children():
        hijo.grid_configure(padx=12, pady=8)

    v0.geometry('+70-70')
    v0.focus()
    v0.mainloop()
    
    return()



#### volcado (para probar) ####

def salida_terminal():
    global n_f, nfcompleto
    global ptos,tramos,pto_base,Gvec,A_secc,eps_L,Iy,Iz,Iyz,det,esferico
    global thetaPR,psi,eta,Wpsi,Weta,ppsi,peta,E,J_secc,I_a,cerrado,theta_secc
    global uT_medio,nucleo,C_tensor
    global n_ptos,n_tramos,yzsonPR,Ipsi,Ieta,v_psi,v_eta

    if esferico:
        print ('\nEl tensor de inercia parece esferico. Se toma como tal')
        print ('(todas las direcciones son principales, con igual inercia).')

    if (not esferico):
        if yzsonPR :
            print ('\nParece que yz son ejes principales. Se toman como tales')


    print('\n  debido a Vy=1:')
    print('ipto          ux            deseq')
    for ip, pto in ptos.items():
        a=0.
        for it in pto.inciP:
            tr=tramos[it]
            a += tr.qVy[tr.inciT.index(ip)]
        print('{:3d} {:16.8e} {:16.8e} '.format(ip, pto.uVy, a))

    print('\ni_tr      flujo_i         flujo_j ')
    for it,tr in tramos.items():
        print('{:3d} {:16.8e} {:16.8e} '.format(it, tr.qVy[0], tr.qVy[1]))

    print('\n  debido a Vz=1:')
    print('ipto          ux            deseq')
    for ip, pto in ptos.items():
        a=0.
        for it in pto.inciP:
            tr=tramos[it]
            a += tr.qVz[tr.inciT.index(ip)]
        print('{:3d} {:16.8e} {:16.8e} '.format(ip, pto.uVz, a))

    print('\ni_tr      flujo_i         flujo_j ')
    for it,tr in tramos.items():
        print('{:3d} {:16.8e} {:16.8e} '.format(it, tr.qVz[0], tr.qVz[1]))

    print('\n\n###  C.E.C:  Ey, Ez  ###')
    print (E[0],' , ',E[1])
    print('en ejes dados es:')
    print(E[0]+Gvec[0], ' , ',  E[1]+Gvec[1])

    print('\n\nAngulo unitario de torsion = ', theta_secc)
    print('Modulo de torsion J = ', J_secc)

    print('\nDebido a T=1:')
    print('pto            uT             deseq')
    for ip, pto in ptos.items():
        a=0.
        for it in pto.inciP:
            tr=tramos[it]
            a+=tr.qT[tr.inciT.index(ip)] 
        print('{:3d} {:16.8e} {:16.8e} '.format(ip, ptos[ip].uT, a))

    print('\nEl promedio de uT vale: ', uT_medio)

    print('\nEl modulo de alabeo Ia= ', I_a)

    print('\n Tensor de giros a cortante C_ij: ')
    print(np.array(C_tensor))

    print ('\n\n .../ fin del resumen por terminal /...\n\n')
    return()
    
 

# funcion para poner props basicas de tramo a partir de los datos iniciales
def basico_tr():
    global ptos, tramos

    for it,tr in tramos.items():
        ini,fin=tr.inciT
        medy,medz= (ptos[fin].y+ptos[ini].y)/2,(ptos[fin].z+ptos[ini].z)/2
        incy,incz= (ptos[fin].y-ptos[ini].y) , (ptos[fin].z - ptos[ini].z)
        d=np.sqrt(incz*incz+incy*incy)

        if tr.tipoT==0:
            L=d
            alfa= np.arctan2 (incz, incy)
            if alfa <0.: alfa += 2*np.pi
            tr.L=L
            tr.alfa=alfa
            tr.GT=[medy,medz]
        else:
            R, e, signo = tr.R, tr.e, tr.signo
            a=2*np.arcsin(d/(2*R))
            alfa = signo*a
            if abs(alfa)<1.e-5:
                texto= 'El arco entre los puntos {:3d},{:3d} es muy '
                texto+='pequeño (tramo {:3d}). Revise los datos.'
                texto+='\n             '+ str(sys.exc_info()[0])
                messagebox.showinfo(message=texto.format(ini,fin,it),
                    title='MAL',parent=v9)
                return(1)
            L=R*a
            yC= medy - np.cos(alfa/2)*signo*R*incz/d
            zC= medz + np.cos(alfa/2)*signo*R*incy/d
            alfai=np.arctan2(ptos[ini].z-zC, ptos[ini].y-yC)
            if alfai <0.: alfai += 2*np.pi
            
            tr.L=L
            tr.alfa=alfa
            tr.alfai=alfai
            tr.yC, tr.zC = yC, zC

            a = np.sin(alfai+alfa)-np.sin(alfai)
            b = np.cos(alfai+alfa)-np.cos(alfai)
            tr.GT=[yC + a*R/alfa, zC - b*R/alfa]


# funcion que da sxx en un punto de la seccion
def sxx_yz(Nx,My,Mz,y,z):
    global Iy, Iz, Iyz, det, A_secc
    a= (Mz*Iy-My*Iyz)*y
    b= (My*Iz-Mz*Iyz)*z
    sxx= Nx/A_secc + (a+b)/det
    return(sxx)
    



#       ################################
#       ####                        ####
#       ####    MOTOR DE CALCULO    ####
#       ####                        ####
#       ################################

def motor_calculo():
    global n_f, nfcompleto, anchoy, anchoz
    global ptos,tramos,pto_base,Gvec,A_secc,eps_L,Iy,Iz,Iyz,det,esferico
    global thetaPR,psi,eta,Wpsi,Weta,ppsi,peta,E,J_secc,I_a,cerrado,theta_secc
    global uT_medio,nucleo,C_tensor, C0y,C0z, uVy_medio, uVz_medio
    global n_ptos,n_tramos,yzsonPR,Ipsi,Ieta,v_psi,v_eta


    def haz_nucleo():
        global anchoy,anchoz,eps_L ,yzsonPR,nucleo
        
        candidatas =[]  # lista de lineas neutras candidatas
        tics = []       # [yG, zG, i] de los ptos & de algunos sobre tramos curvos.
                        # i es un indice =-3333 si el tic proviene de un pto[],
                        # i= itramo si el tic proviene de un tramo. Permitira no
                        # formar LNs entre tics de la misma curva
                        # NO ES GLOBAL porque no es equiespaciado (otros muestreos
                        # como el de C_tensor lo requieren)

        def valida(linea):
            # devuelve True si todos los tics estan al mismo lado de la linea
            # tics[] esta en el scope de esta funcion, con ello se cuenta
            signo=0
            for a in tics:
                p=[a[0],a[1]]
                d=linea.dist(p)
                if abs(d) > eps_L:
                    nuevo=np.sign(d)
                    if (signo==1 and nuevo==-1) or (signo==-1 and nuevo==1):
                        return (False)
                    elif signo==0:
                        signo = nuevo # cuando caze un signo ya no lo cambiara
            return(True)
                        

        # tics en todos los ptos[]
        for pt in ptos.values():  tics.append([pt.yG, pt.zG, -3333])
        
        # & tic cada 5º aprox en tramos curvos
        for it, tr in tramos.items():  
            if tr.tipoT == 1:
                n = abs( int (tr.alfa/0.07) )
                if n != 0:
                    for a in np.linspace (0,tr.alfa,n): # tics intermedios & extremos
                        c1, s1 = np.cos(tr.alfai+ a) , np.sin(tr.alfai+ a)
                        y= tr.yC + tr.R * c1
                        z= tr.zC + tr.R * s1
                        tics.append([y,z,it])
                        
                        linea = Recta() # aprovechamos para crear LN candidadas de tics
                        linea.p=[y,z]
                        linea.en=np.array([c1, s1])*np.sign(tr.alfa)  # hacia la dcha del tr
                        linea.es=np.array([-s1,c1])*np.sign(tr.alfa)  # hacia el avance de alfa
                        candidatas.append(linea)
        
        # asignar ancho y,z &  eps_L. Son globales, para luego dibujar etc
        anchoy, anchoz = 0., 0.
        for i in range (len(tics)-1):
            for j in range (i+1, len(tics)):
                d= abs( tics[i][0]-tics[j][0] )
                if d > anchoy: anchoy=d
                d= abs( tics[i][1]-tics[j][1] )
                if d > anchoz: anchoz=d
        eps_L=(anchoy + anchoz)/1.e5
        
        # limpiar los tics repetidos (los habra)
        i, n = 0, len(tics)
        while i<n:
            a=[]
            for j in range(i+1,n):
                d= np.sqrt((tics[j][0]-tics[i][0])**2+(tics[j][1]-tics[i][1])**2)
                if d < eps_L:
                    a.append(j)
            for j in sorted(a,reverse=True): del tics[j]
            n=len(tics)
            i += 1

        # validar las LN de tics que hemos propuesto
        a=[]
        for i in range(len(candidatas)): 
            if not valida(candidatas[i]): a.append(i)
        
        for i in sorted(a,reverse=True): del candidatas[i]

        # creamos LN candidata con cada par de tics que no sean del mismo tr curvo
        for i in range(len(tics)-1):
            for j in range(i+1,len(tics)):
                if (tics[i][2] != tics[j][2]) or ((tics[i][2]+tics[j][2]))<0 :
                    r=Recta()
                    r.p=[tics[i][0],tics[i][1]]
                    incy=tics[j][0] - tics[i][0]
                    incz=tics[j][1] - tics[i][1]
                    L=np.sqrt((incy)**2+(incz)**2)
                    r.es=[incy/L, incz/L]   # en[] no hara falta 
                    if valida(r): candidatas.append(r)
        
        # Las LN candidatas ya son firmes. Puede haber repeticiones (a borrar),  
        # lo hago controlando (evitando) repeticiones de puntos de presion:
            
        a=[]    # sera el nucleo sin ordenar
        for r in candidatas:
            c= ( -r.p[0]*r.es[1] + r.p[1]*r.es[0] ) * A_secc
            ppy = ( Iz*r.es[1] -Iyz*r.es[0]) / c
            ppz = (-Iy*r.es[0] +Iyz*r.es[1]) / c
            poner=True
            for p in a:
                d=np.sqrt( (p[0]-ppy)**2 + (p[1]-ppz)**2)
                if d<eps_L: # hay ya un p.pres. con la misma posicion
                    poner=False
                    break
            if poner: a.append([ppy,ppz])
        
        # Los puntos del nucleo estan en a[]. Los ordeno en antihorario desde G
        b=[]        # tendra angulo 0 -> 360 del pto presion 
        for pp in a:
            theta= np.arctan2(pp[1],pp[0])
            if theta < 0.: theta += 2*np.pi
            b.append(theta)

        nucleo=[x for _,x in sorted(zip(b,a))]      # toma comando!
        nucleo.append(nucleo[0])            # para cerrar el nucleo
        
        return(nucleo)






    ############################
    ### comienzo de calculos ###
    ############################
    
    # rellenamos propiedades basicas de tramos a partir de tipoT,i,j,e,R,signo
    if basico_tr():return() # si hay error lo ha dicho ya, que salga

    # calcular el centro de areas de la seccion Gvec[]
    sumy, sumz, A_secc, Gvec = 0., 0. , 0. , [0.,0.]
    for tr in tramos.values():
        a = tr.L * tr.e
        A_secc += a
        sumy += a * tr.GT[0] 
        sumz += a * tr.GT[1]
    Gvec[0]=sumy/A_secc
    Gvec[1]=sumz/A_secc

    # cambiar a origen Gvec coordenadas de ptos, centros, GT. Aprovecho para 
    # asignar eps_L provisionalmente (se reasigna mejor en haz_nucleo)

    ymax=1.e-55
    for pto in ptos.values():
        pto.yG = pto.y-Gvec[0] # OJO OJO: CENTROS DESDE G SOLO SE MANTENDRA
        pto.zG = pto.z-Gvec[1] # EN EL MOTOR DE CALCULO. LUEGO SE LLAMA A
                               # haz_ficha QUE VUELVE A LLAMAR A basico_tr
                               # Y VUELVE A PONER LOS CENTROS DESDE O 
        if ymax < abs(pto.yG): ymax=abs(pto.yG)   # hale, sin complicarlo mas
    eps_L = ymax /4.e4
    #print('\n\n Primera asignacion de eps_L= ', eps_L)
             
    for tr in tramos.values():
        tr.GT[0] -= Gvec[0]
        tr.GT[1] -= Gvec[1]
        if tr.tipoT == 1:
            tr.yC -= Gvec[0]
            tr.zC -= Gvec[1]


    for p in ptos.values(): p.inciP=[] # no se pq hace falta esto.
                                         # en la 0.95 no hacia falta
    # construir incidencias de ptos inciP[]
    for it, tr in tramos.items():
        for j in tr.inciT:
            ptos[j].inciP.append(it)

    # calcular los momentos de inercia de la seccion - ejes por G
    Iy, Iz, Iyz = 0., 0., 0.
    for tr in tramos.values():
        i,j= tr.inciT
        yi=ptos[i].yG
        zi=ptos[i].zG
        L=tr.L
        e=tr.e
        if tr.tipoT== 0:
            salfa= np.sin(tr.alfa)
            calfa= np.cos(tr.alfa)
            Iz  += (yi*yi+ ((calfa*L)**2)/3.+ yi*calfa*L)*e*L
            Iy  += (zi*zi+ ((salfa*L)**2)/3.+ zi*salfa*L)*e*L
            Iyz += (yi*zi+ (yi*salfa + zi*calfa)*L/2.+ 
                (salfa*calfa*L*L)/3.)*e*L
        elif tr.tipoT==1:
            alfa, alfai, R, yC, zC= tr.alfa, tr.alfai, tr.R, tr.yC, tr.zC
            # las coordenadas del centro yC,zC ya estan en ejes por G
            salfai, salfai2= np.sin(alfai), np.sin(2*alfai)
            calfai, calfai2= np.cos(alfai), np.cos(2*alfai)
            s1y1, s2y2= np.sin(alfai+alfa), np.sin(2*alfai +2*alfa)
            c1y1, c2y2= np.cos(alfai+alfa), np.cos(2*alfai +2*alfa)
            # si alfa es negativo, ds=-R d(alfa), no + como en la formula
            # => Iy, Iz saldran negativos. Hay que cambiar signos:
            
            a= zC*zC*alfa+ R*R*(alfa- s2y2/2+ salfai2/2)/2- 2*R*zC*(c1y1-calfai)
            Iy += a*R*e *np.sign(alfa)
            
            a= yC*yC*alfa+ R*R*(alfa+ s2y2/2- salfai2/2)/2+ 2*R*yC*(s1y1-salfai)
            Iz += a*R*e *np.sign(alfa)
            
            a= yC*zC*alfa- R*yC*(c1y1-calfai)+ R*zC*(s1y1-salfai)
            a   += R*R*(s1y1*s1y1- salfai*salfai)/2.
            Iyz += a*R*e *np.sign(alfa)

    det= Iy*Iz-Iyz**2 
     

    # Calcular las inercias principales y direcciones principales.
    a= np.sqrt(0.25*(Iy-Iz)**2+Iyz**2 )  
    b=Iy+Iz+abs(Iyz) 
    esferico=False
    if a < b*1.e-5:
    #    El tensor de inercia parece esferico. Se toma como tal
        Iy=Iz
        Iyz=0.
        Ipsi,Ieta = Iz,Iz
        v_eta , v_psi = [0,1], [1,0]
        thetaPR=0.
        esferico=True
    else:
        # Notese que el tensor de inercia J tiene Jyy=Izz y viceversa,
        # por ello los autovectores salen intercambiados para usar I.
        # Llamaremos eta al eje pr fuerte (psi al debil)
        a=[[Iz,Iyz],[Iyz,Iy]]
        [Ipsi, Ieta], [v_eta, v_psi]= np.linalg.eigh(a)
        # thetaPR <90º es desde z hacia el ejePR fuerte eta. Positivo antihorario. 
        if v_eta[1]<0. : v_eta= -v_eta    # que v_eta mire p'arriba (como z)
        if v_psi[0]<0. : v_psi= -v_psi    # que v_psi mire pa la derecha (como y)
        
        if   abs(v_eta[0]) < 1.e-6:
            # yz parecen principales. Se toman como tales
            yzsonPR=True
            thetaPR= 0.
        elif abs(v_eta[1]) < 1.e-6:
            # yz parecen principales. Se toman como tales
            yzsonPR=True
            thetaPR=np.pi /2
        else:
            yzsonPR=False
            thetaPR=np.arctan(-v_eta[0]/v_eta[1])

    psi= Recta( [0.,0.], v_psi, [0.,0.] )
    eta= Recta( [0.,0.], v_eta, [0.,0.] )


    # calcular los modulos resistentes (ejes por G)
    deta, dpsi = 0., 0.
    peta, ppsi = [0.,0.] , [0.,0.]
    for ip,pt in ptos.items():
        y,z= pt.yG, pt.zG
        d= abs( eta.dist([y,z]) )
        if d > deta:
            deta = d
            peta=[y,z]
        d= abs( psi.dist([y,z]) )
        if d > dpsi:
            dpsi = d
            ppsi=[y,z]

    Weta = Ieta / deta
    Wpsi = Ipsi / dpsi


    # calcular los momentos estaticos de los tramos
    for iT, tr in tramos.items():
        a= tr.L * tr.e
        tr.Qy = a*tr.GT[1]
        tr.Qz = a*tr.GT[0] # mejor que andar con funciones
        
        tr.iQy= tr.iQy_s(tr.L)
        tr.iQz= tr.iQz_s(tr.L) # se discrimina alli segun el tipoT 
    
    print('\n iTramo       Qy         Qz        iQy       iQz')
    for iT, tr in tramos.items():
        print('{:3d} {:10.5f} {:10.5f} {:10.5f} {:10.5f} '.format(
                                    iT, tr.Qy, tr.Qz, tr.iQy, tr.iQz))

    # correspondencias de nudos quitando el de ux=0 (para el solver)
    corr=[]
    for i in ptos.keys():
        if i !=  pto_base:
            corr.append(i)
    n=len(corr)

    # inicializamos flujos en extremos de tramo como si fuese ux=0 en los ptos
    for tr in tramos.values():
        # para Vy=1
        a= -(tr.iQz*Iy-tr.iQy*Iyz)/(det*tr.L)
        b= -a -(tr.Qz*Iy-tr.Qy*Iyz)/det
        tr.qVy=[a,b]    
        # para Vz=1
        a= -(tr.iQy*Iz-tr.iQz*Iyz)/(det*tr.L)
        b= -a - (tr.Qy*Iz-tr.Qz*Iyz)/det
        tr.qVz=[a,b]


    ##############################
    ##  resolver caso de Vy=1   ##
    ##############################

    a=np.zeros((n,n))
    b=np.zeros((n))
    ptos[pto_base].uVy=0.

    for tr in tramos.values():
        ip,jp=tr.inciT
        c= tr.e/tr.L
        if ip != pto_base:  # ecuacion del pto i
            b[corr.index(ip)] -= tr.qVy[0]
            a[corr.index(ip),corr.index(ip)] += c
            if jp != pto_base:
                a[corr.index(ip),corr.index(jp)] -= c
        if jp != pto_base:  # ecuacion del pto j
            b[corr.index(jp)] -= tr.qVy[1]
            a[corr.index(jp),corr.index(jp)] += c
            if ip != pto_base:
                a[corr.index(jp),corr.index(ip)] -= c
        
    x=np.linalg.solve(a,b)

    # trasladar la solucion ux obtenida a los puntos
    for ip in corr:
        ptos[ip].uVy=x[corr.index(ip)]

    # corregir los flujos de extremo de tramo con los terminos en ux
    for tr in tramos.values():
        ip,jp= tr.inciT
        c= (tr.e / tr.L)*(ptos[ip].uVy-ptos[jp].uVy)
        tr.qVy[0] += c
        tr.qVy[1] -= c


    ##############################
    ##  resolver caso de Vz=1   ##
    ##############################

    # a es la misma, no se recalcula
    b=np.zeros((n))
    ptos[pto_base].uVz=0.

    for tr in tramos.values():
        ip,jp=tr.inciT
        if ip != pto_base:  # ecuacion del pto i
            b[corr.index(ip)] -= tr.qVz[0]
        if jp != pto_base:  # ecuacion del pto j
            b[corr.index(jp)] -= tr.qVz[1]
    # seria facil resolver ambos de una pasada?  Quiza mas tarde...

    x=np.linalg.solve(a,b)

    # trasladar la solucion ux obtenida a los puntos
    for ip in corr:
        ptos[ip].uVz=x[corr.index(ip)]

    # corregir los flujos de extremo de tramo con los terminos en ux
    for tr in tramos.values():
        ip,jp= tr.inciT
        c= (tr.e / tr.L)*(ptos[ip].uVz-ptos[jp].uVz)
        tr.qVz[0] += c
        tr.qVz[1] -= c


    ###########################################
    ##  calcular centro de esf cortantes E   ##
    ###########################################

    E=[0.,0.]

    # para Vy=1 encontraremos zE=E[1] con origen G (con Vz=1 encontramos yE=E[0])
    MG_Vy, MG_Vz = 0., 0.
    for it,tr in tramos.items():
        # calcular el Mto respecto de G de las tensiones del tramo
        # positivo si sale de la figura (sentido +x)
        a,b = tr.MG()
        MG_Vy += a
        MG_Vz += b
        
    E[1]= -MG_Vy # es zE= -M/Vy con Vy=1
    E[0]=  MG_Vz # es yE= +M/Vz con Vz=1


    #############################
    #### Calculos de torsion ####
    #############################


    J_abierto=0.
    for tr in tramos.values():
        J_abierto += tr.L*(tr.e**3)
    J_abierto /= 3.

    for tr in tramos.values():  tr.Asect= tr.Asect_s(tr.L)

    ptos[pto_base].uT=0.   # el ux del pto_base se toma nulo

    # lo hago diferente para cerrados y abiertos porque no quiero dejar
    # que se resuelva con una matriz mal condicionada

    if len(tramos)+1 > len(ptos):
        print('\n\nT: se resuelve como cerrado')
        cerrado=True
        
        n= len(corr)+ 1     # esa ultima sera la ecuacion de Suma de Mtos
        a= np.zeros((n,n))
        b= np.zeros((n))
        n -= 1              # apunta a la ultima ec. (suma mtos) en a[]

        for tr in tramos.values():  # montar el sistema de ecuaciones
            i,j= tr.inciT
            if i != pto_base: i_ = corr.index(i)
            if j != pto_base: j_ = corr.index(j)
            c= tr.e/tr.L
            
            if j != pto_base:
                # a la ecuacion deseq pto j:
                a[j_,j_] += c
                a[j_,n]  += c*tr.Asect
                if i != pto_base:
                    a[j_,i_] -= c
                
            if i != pto_base:
                # a la ecuacion deseq pto i:
                a[i_,i_] += c
                a[i_,n]  -= c*tr.Asect
                if j != pto_base:
                    a[i_,j_] -= c

                # a la ec. de Suma Mtos:
            a[n,n]  += c*tr.Asect**2
            if j != pto_base:   a[n,j_] += c*tr.Asect   
            if i != pto_base:   a[n,i_]  -= c*tr.Asect

        # Las rigideces de abierto se incluyen porque es lo mas correcto. La opcion
        # de no incluirlas generaba confusion sin un beneficio claro:
        a[n,n] += J_abierto

        b[n] = 1.    # El Mx aplicado (llamado T en el programa)
        
        x=np.linalg.solve(a,b)

        # trasladar la solucion ux obtenida a los puntos, y ang. giro. unitario
        theta_secc= x[n] 
        for i in corr:
            i_ = corr.index(i)
            ptos[i].uT=x[i_]

        # asignar flujos de extremo de tramo & los no-flujos

        for tr in tramos.values():
            i,j= tr.inciT
            c= tr.e / tr.L
            b = ( ptos[j].uT -ptos[i].uT + tr.Asect*theta_secc )*c
            tr.qT = [-b,b]
            tr.tau_bi= theta_secc * tr.e

        # calculo J de la seccion
        J_secc = 1./theta_secc

    elif len(tramos)+1 == len(ptos):
        print('\n\n T:se resuelve como abierto')
        cerrado=False   # (es abierto)
        
        J_secc =J_abierto
        theta_secc= 1. /J_abierto
        n= len(corr)
        a= np.zeros((n,n))
        b= np.zeros((n))
        
        # en este caso tengo una ec por tramo (su incr. ux es proporcional
        # a su area sectorial). Como es n_ptos-1=tramos, es lo justo para 
        # resolver los ux. Y con matriz bien condicionada!:
        
        corr_tr=[]
        for itr in tramos.keys():
            corr_tr.append(itr) # correlacion entre tramo y fila de a[]
        
        for itr,tr in tramos.items():
            itr_ = corr_tr.index(itr)
            tr.qT= [0., 0.]
            tr.tau_bi= - theta_secc*tr.e
            i,j= tr.inciT
            if i != pto_base: 
                i_ = corr.index(i)
                a[ itr_, i_]= 1
            if j != pto_base:
                j_ = corr.index(j)
                a[ itr_, j_]= -1
                
            b[itr_]= tr.Asect * theta_secc
        
        x=np.linalg.solve(a,b)
        
        # trasladar los ux de torsion (uT) obtenidos a los puntos
        for i in corr:
            i_ = corr.index(i)
            ptos[i].uT=x[i_]


    else:
        print('El numero de ptos no puede ser mayor que el de tramos +1')
        print('Abortando')
        quit()


    ###########################################
    #####  calculo del modulo de alabeo  ######
    ###########################################

    # lo primero necesito el uT_medio, y para ello integr u(s)*dA, 
    # y para ella int 2w(s) ds , que en el tramo llamo c:

    uT_medio= 0.
        
    for tr in tramos.values():
        i,j = tr.inciT
        e, L = tr.e, tr.L
        if tr.tipoT == 0:
            yi, zi = ptos[i].yG, ptos[i].zG
            yj, zj = ptos[j].yG, ptos[j].zG
            es_y, es_z = (yj-yi)/L, (zj-zi)/L
            c= (yi-E[0])*es_z - (zi-E[1])*es_y
            c *= L*L/2
        else:
            alfa, alfai, R, yC, zC = tr.alfa, tr.alfai, tr.R, tr.yC, tr.zC
            signo = np.sign(alfa)
            c  = (yC-E[0])*(-np.cos(alfai+alfa)+np.cos(alfai)-alfa*np.sin(alfai))
            c += (zC-E[1])*(-np.sin(alfai+alfa)+np.sin(alfai)+alfa*np.cos(alfai))
            c += R*alfa**2 /2
            c *= R*R*signo
        # c tiene integr 2*area_sectorial del tramo. Completo int u(s)*dA:
        c = e*L*ptos[i].uT + tr.qT[1]*L*L/2 - c*e*theta_secc
        uT_medio += c

    uT_medio /= A_secc


    # ya esta uT_medio. Ahora I_a= int_seccion (u_centrado/theta_secc)**2 *e*ds
    I_a = 0.
    for tr in tramos.values():
        I_a += tr.e * tr.int2()  # con int2() esto es asi de facil

    I_a /= theta_secc**2
        
    #########################################
    #####  calculo del nucleo central  ######
    #########################################

    nucleo = haz_nucleo()   # todo se hace alli. Y se captan dimensiones & eps_L


    #########################################
    #####  tensor de giros a cortante  ######
    #########################################

    # El C_tensor requiere el muestreo en una nube de puntos equiespaciados en la
    # seccion. Chapuzas como dejar que un punto figure dos veces (los extremos
    # de tramo) son comodas de programar pero hacen que un termino del tensor
    # que en casos de simetria clara deberia ser =0, y esperamos e-15, quede 
    # tipo e-3 o e-4 (orden de 1/n_tics en la seccion). Por eso se adopta la 
    # estrategia de que los tics sean interiores al tramo, entendiendo que hay que
    # añadir manualmente o como sea tics en cada pto[].

    # uVy & uVz_medios se calculan tambien numericamente y aprovecharan este array 
    # de Tramo(). 
    # Va a usarse tambien para dibujar, y aqui viene bien que haya repeticiones
    # tipo conexion de tramos. Los arrays de dibujar seran uVy_plot, uVz_plot, 
    # qVy_plot, qVz_plot, uT_plot. Los qT, de torsion son ctes y se haran
    # al vuelo, asi como los calculos de sxx.

    # n_tics sera el num aprox de muestras no-nodales en la seccion. Asi queda
    # facil de cambiar. Incluso se puede hacer accesible al usuario (**pensar**)

    # asignacion de n_tics de forma que el tramo menor tenga 4+ tics interiores
    Ltot, Lmin= 0., 1.e90
    for tr in tramos.values():
        Ltot += tr.L
        if tr.L<Lmin: Lmin=tr.L
    n_tics=max(65, int(4.1*Ltot/Lmin))
    print('\n n_tics nominal= ',n_tics)

    uVy_medio, uVz_medio = 0., 0.
    a,b,c = [], [], []
    i=0     # contador de tics de verdad

    print()
    for itr,tr in tramos.items():
        tr.tics, tr.uVy_plot, tr.uVz_plot = [],[],[]
        tr.qVy_plot, tr.qVz_plot, tr.uT_plot = [],[],[]
        n= int(n_tics*tr.L/Ltot)
        '''if n<3: # esto viene mal para Cij, estropea la uniform. -> comprobar JC
            n=3
            print('\nEl tramo ',itr,' es muy pequeño. \nSe le dan 3 tics internos.')
        '''
        print('tramo ',itr,' con ',n,' tics interiores.')
        s_tic=tr.L/(n+1.)
        for s in np.linspace (s_tic,  tr.L - s_tic,  n): 
            i += 1
            y,z= tr.yz_s(s)
            uVy, uVz = tr.uV_s(s, 1., 0.), tr.uV_s(s, 0., 1.)
            uT=tr.uT_s(s)
            qVy, qVz = tr.qV_s(s, 1., 0.), tr.qV_s(s, 0., 1.)
            
            tr.tics.append([y,z,s])
            tr.uVy_plot.append(uVy)
            tr.uVz_plot.append(uVz)
            tr.uT_plot.append(uT)
            tr.qVy_plot.append(qVy)
            tr.qVz_plot.append(qVz)
            a.append([1.,y,z])
            b.append(uVy)
            c.append(uVz)
            uVy_medio += uVy
            uVz_medio += uVz

    for ip, pt in ptos.items(): 
        # quedan por sumar valores en ptos (no-tics de momento)
        # solo para C_tensor & para los uV medios
        i += 1
        y,z = pt.y, pt.z
        a.append([1.,y,z])
        b.append(pt.uVy)
        c.append(pt.uVz)
        uVy_medio += pt.uVy
        uVz_medio += pt.uVz

    print('\n\n n_tics real para C[] & uV medios = ',i)
    uVy_medio /= float(i)
    uVz_medio /= float(i)

    C_tensor=np.zeros((2,2))
    C0y, C_tensor[0][0],C_tensor[0][1] = np.linalg.lstsq(a,b,rcond=-1)[0]
    C0z, C_tensor[1][0],C_tensor[1][1] = np.linalg.lstsq(a,c,rcond=-1)[0]

    a,b,c = None, None, None    # liberamos ese espacio
        
        
    for it,tr in tramos.items():
        # para los plots añadimos i & j de tramo a sus tics[], que viene bien
        i,j=tr.inciT
        y,z = ptos[i].y , ptos[i].z
        uVy, uVz, uT = ptos[i].uVy, ptos[i].uVz, ptos[i].uT
        qVy, qVz = tr.qVy[0], tr.qVz[0]

        tr.tics.insert( 0, [y,z, 0.] )
        tr.uVy_plot.insert(0, uVy)
        tr.uVz_plot.insert(0, uVz)
        tr.uT_plot.insert (0, uT)
        tr.qVy_plot.insert(0, -qVy)  # recordar, q del solver es + hacia el pto
        tr.qVz_plot.insert(0, -qVz)
        
        y,z = ptos[j].y , ptos[j].z
        uVy, uVz, uT = ptos[j].uVy, ptos[j].uVz, ptos[j].uT
        qVy, qVz = tr.qVy[1], tr.qVz[1]

        tr.tics.append( [y,z, tr.L] )
        tr.uVy_plot.append(uVy)
        tr.uVz_plot.append(uVz)
        tr.uT_plot.append(uT)
        tr.qVy_plot.append(qVy)
        tr.qVz_plot.append(qVz)

    # quiero usar numpy arrays para graficos (funciona el slice etc):
    for tr in tramos.values():
        tr.tics=     np.array(tr.tics)
        tr.uVy_plot= np.array(tr.uVy_plot)
        tr.uVz_plot= np.array(tr.uVz_plot)
        tr.uT_plot=  np.array(tr.uT_plot)
        tr.qVy_plot= np.array(tr.qVy_plot)
        tr.qVz_plot= np.array(tr.qVz_plot)

    # Ahora esta almacenada la info de dibujar uV & qV en los objetos Tramo & Pto 

    return()



###################################
####  Cargar un prb existente  ####
###################################

def a_cargar():
    global elige, n_f, nfcompleto, ptos, tramos, n_ptos, n_tramos, pto_base

    if elige=='default':
        n_ptos, n_tramos, pto_base =7,8,1
        ptos, tramos = {}, {}
        for ip in range(1,8): ptos[ip]=Pto()
        ptos[1].y, ptos[1].z =  -5. , 0.
        ptos[2].y, ptos[2].z =   5. , 0.
        ptos[3].y, ptos[3].z = -10. , 10.
        ptos[4].y, ptos[4].z =  -5. , 10.
        ptos[5].y, ptos[5].z =   5. , 10.
        ptos[6].y, ptos[6].z =  10. , 10.
        ptos[7].y, ptos[7].z =   0. , 14.142

        for it in (1,2,4,5,6,7,8,9): tramos[it]=Tramo()
        tramos[1].tipoT, tramos[1].inciT= 0, [1,4]
        tramos[2].tipoT, tramos[2].inciT= 0, [1,2]
        tramos[4].tipoT, tramos[4].inciT= 0, [2,5]
        tramos[5].tipoT, tramos[5].inciT= 0, [3,4]
        tramos[6].tipoT, tramos[6].inciT= 0, [4,5]
        tramos[7].tipoT, tramos[7].inciT= 0, [5,6]
        tramos[8].tipoT, tramos[8].inciT= 1, [3,7]
        tramos[9].tipoT, tramos[9].inciT= 1, [6,7]
        
        for it in (1,2,4,5,6,7,8,9): tramos[it].e= 0.5
        tramos[8].R, tramos[8].signo= 14.142, -1.
        tramos[9].R, tramos[9].signo= 14.142,  1.
        
        entry_My.delete(0,'end')
        entry_Mz.delete(0,'end')
        entry_Nx.delete(0,'end')
        entry_Vy.delete(0,'end')
        entry_Vz.delete(0,'end')
        entry_T.delete(0,'end')

        entry_My.insert(0,'22.')
        entry_Mz.insert(0,'33.')
        entry_Nx.insert(0,'4.')
        entry_Vy.insert(0,'-2.')
        entry_Vz.insert(0,'7.')
        entry_T.insert(0,'10.')
        
    else:
        nfcompleto=filedialog.askopenfilename(parent=v9, initialdir=dir_home,
            title='Abrir archivo')
        if not nfcompleto: return()
        
        n_f= path.basename(nfcompleto)
        v9.title('ThinSecBeam 1.0 - '+n_f)
        
        ptos, tramos = {}, {}

        f=open(nfcompleto, 'r')
        
        a=f.readline()
        a=f.readline()
        a=f.readline()
        a=f.readline()
        b=a.split()
        n_ptos, n_tramos, pto_base = int(b[0]), int(b[1]), int(b[2])
        
        a=f.readline()
        a=f.readline()
        for i in range(n_ptos):
            a=f.readline()
            b=a.split()
            ip,y,z= int(b[0]), float(b[1]), float(b[2])
            ptos[ip]=Pto()
            ptos[ip].y, ptos[ip].z = y, z

        a=f.readline()
        a=f.readline()
        for i in range(n_tramos):
            a=f.readline()
            b=a.split()
            it,tipoT,ip,jp,e,=int(b[0]),int(b[1]),int(b[2]),int(b[3]),float(b[4])
            tramos[it]=Tramo()
            tr=tramos[it]
            tr.tipoT, tr.inciT, tr.e = tipoT, [ip,jp], e
            if tipoT == 1:
                R, s = float(b[5]), b[6]
                signo= -1. if s == '-' else 1.
                tr.R=R
                tr.signo=signo

        a=str(f.readline)
        f.close()
        
    # rellenamos el GUI con los datos leidos:
    
    for i in range (hard_nfilas):
        for j in range (3):
            filas_ptos[i][j].delete(0,'end')
        for j in range (6):
            filas_tramos[i][j].delete(0,'end')
    
    for i, ip in enumerate(ptos):
        filas_ptos[i][0].insert(0,str(ip))
        filas_ptos[i][1].insert(0,str(ptos[ip].y))
        filas_ptos[i][2].insert(0,str(ptos[ip].z))
    for i, it in enumerate(tramos):
        filas_tramos[i][0].insert(0,str(it))
        filas_tramos[i][1].insert(0,str(tramos[it].inciT[0]))
        filas_tramos[i][2].insert(0,str(tramos[it].inciT[1]))
        filas_tramos[i][3].insert(0,str(tramos[it].e))
        if tramos[it].tipoT == 1:
            filas_tramos[i][4].insert(0,str(tramos[it].R))
            a='-' if tramos[it].signo <0 else '+'
            filas_tramos[i][5].insert(0,a)
    
    boton_txt.state(['disabled'])
    boton_guardar.state(['disabled'])
    boton_goM.state(['disabled'])
    boton_goV.state(['disabled'])
    boton_goT.state(['disabled'])
    boton_yaTr.focus()
    
def a_cargar_mas():
    global elige, n_f,nfcompleto, ptos,tramos, n_ptos,n_tramos, pto_base

    texto ='Esta función permite resolver problemas con mas '
    texto+='tramos de los previstos en la interfaz. Recuerde que '
    texto+='los tramos rectos y circulares se resuelven usando expresiones '
    texto+='analíticas, por lo que la precisión no aumentará poniendo mas.\n\n'
    texto+='Si aun así desea resolver una sección con muchos tramos, '
    texto+='proporcione a continuación el fichero de datos que usted habrá '
    texto+='editado manualmente. El programa calculará con esos datos y '
    texto+='la interfaz (simplificada) le dejará en el punto al que hubiese '
    texto+='llegado tras pulsar el botón ''"Calcula"'' después de introducir los '
    texto+='puntos y tramos.'

    a=messagebox.askokcancel(title='Pregunta', message='¿Está seguro?',
        detail=texto,icon='question', parent=v9)
    print(a)
    if not a: return()
    
    nfcompleto=filedialog.askopenfilename(parent=v9, initialdir=dir_home,
        title='Abrir archivo')
    n_f= path.basename(nfcompleto)
    v9.title('ThinSecBeam 1.0 - '+n_f)
    
    for hijo in frame_ptos.winfo_children(): hijo.destroy()
    for hijo in frame_tramos.winfo_children(): hijo.destroy()
    boton_cargar.destroy()
    boton_guardar.destroy()
    boton_mas.destroy()
    
    ptos,tramos= {}, {}
    
    f=open(nfcompleto,'r')
    
    a=f.readline()
    a=f.readline()
    a=f.readline()
    a=f.readline()
    b=a.split()
    n_ptos, n_tramos, pto_base = int(b[0]), int(b[1]), int(b[2])
    
    a=f.readline()
    a=f.readline()
    for i in range(n_ptos):
        a=f.readline()
        b=a.split()
        ip,y,z= int(b[0]), float(b[1]), float(b[2])
        ptos[ip]=Pto()
        ptos[ip].y, ptos[ip].z = y, z

    a=f.readline()
    a=f.readline()
    for i in range(n_tramos):
        a=f.readline()
        b=a.split()
        it,tipoT,ip,jp,e,=int(b[0]),int(b[1]),int(b[2]),int(b[3]),float(b[4])
        tramos[it]=Tramo()
        tr=tramos[it]
        tr.tipoT, tr.inciT, tr.e = tipoT, [ip,jp], e
        if tipoT == 1:
            R, s = float(b[5]), b[6]
            signo= -1. if s == '-' else 1.
            tr.R=R
    
    motor_calculo()
    haz_ficha(recarga=False)
    salida_terminal()
    return()

    

######################################
####  Informe de texto a fichero  ####
######################################

def haz_informe():
    global n_f, nfcompleto
    global ptos,tramos,pto_base,Gvec,A_secc,eps_L,Iy,Iz,Iyz,det,esferico
    global thetaPR,psi,eta,Wpsi,Weta,ppsi,peta,E,J_secc,I_a,cerrado,theta_secc
    global uT_medio,nucleo,C_tensor
    global n_ptos,n_tramos,yzsonPR,Ipsi,Ieta,v_psi,v_eta

    a='Filename para el informe:'
    nfinforme=filedialog.asksaveasfilename(parent=v9, title=a,
        initialdir= dir_home)
    
    f= open(nfinforme, 'w')


    print('#####  Informe generado por ThinSecBeam 1.0  #####', file=f)
    print('-'*60, file=f)
    
    print('\nIdentificación: ', nfinforme, file=f)
    print('-'*60, file=f)
    print('\nPuntos:', file=f)
    print('Punto base = ', pto_base, file=f)
    print('       _____ejes dados_____    _____ejes por G_____', file=f)
    print('{:>3s} {:>11s} {:>11s} {:>11s} {:>11s}  {:s}'.format(
         'ip','y','z','yG','zG','inciP'), file=f)
    for ip,p in ptos.items():
        print('{:>3d} {:> 11.5g} {:> 11.5g} {:> 11.5g} {:> 11.5g} {:}'.format(
            ip, p.y, p.z, p.yG, p.zG, p.inciP), file=f)
    print('-'*60, file=f)
    
    print('\nTramos - parámetros básicos:', file=f)
    print('{:>4s} {:>6s} {:>6s} {:>6s} {:>11s} {:>11s}'.format(
         'it', 'tipoT','pto_i','pto_j','e','L'), file=f)
    
    for it,tr in tramos.items():
        texto='{:>4d} {:>6d} {:>6d} {:>6d} {:> 11.5g} {:> 11.5g}'
        print (texto.format(it, tr.tipoT, tr.inciT[0], tr.inciT[1],tr.e, tr.L), file=f)
    
    print('\nTramos - parámetros geométricos:', file=f)
    print('{:^4s} {:>11s} {:>11s} {:>11s} {:>11s} {:>11s} '.format(
        'it','alfa','alfai','R','yC','zC'), file=f)
    
    for it,tr in tramos.items():
        if tr.tipoT == 1:
            texto='{:^4d} {:> 11.5g} {:> 11.5g} {:> 11.5g} {:> 11.5g} {:> 11.5g} '
            print (texto.format(it, tr.alfa, tr.alfai, tr.R, tr.yC, tr.zC), file=f)
        else:
            print('{:^4d} {:> 11.5g}'.format(it, tr.alfa), file=f)
    
    print('\nTramos - parámetros de cálculo:', file=f)
    print('{:>4s} {:>11s} {:>11s} {:>11s} {:>11s} {:>11s} '.format(
         'it', 'Qy', 'Qz', 'iQy', 'iQz', 'Asect'),file=f)
    for it,tr in tramos.items():
        texto='{:^4d} {:> 11.5g} {:> 11.5g} {:> 11.5g} {:> 11.5g} {:> 11.5g}'
        print (texto.format(it, tr.Qy, tr.Qz, tr.iQy, tr.iQz, tr.Asect), file=f)
    print('-'*60, file=f)
    
    print('\nSección - parámetros globales:\n', file=f)

    print('{:>11s} {:>11s} {:>11s} {:>11s} {:>11s}'.format(
        'yG','zG','A','ACy','ACz'),file=f)
    print('{:> 11.5g} {:> 11.5g} {:> 11.5g} {:> 11.5g} {:> 11.5g}'.format(
        Gvec[0],Gvec[1],A_secc,1/C_tensor[0][0], 1/C_tensor[1][1]),file=f)
    if (  abs(C_tensor[0][0])+abs(C_tensor[1][1])<
         (abs(C_tensor[1][0])+abs(C_tensor[0][1]))*1.e3  ):
        texto= '     Aviso: Las areas a cortante no parecen tener sentido,\n'
        texto+='     al menos en ejes yz.'
        print(texto,file=f)
    else:
        print('Las areas a cortante parecen tener sentido en ejes yz.',file=f)

    print('\n{:>11s} {:>11s} {:>11s} {:>11s} {:>11s} {:>8s}'.format(
        'Iy','Iz','Iyz','Ipsi','Ieta', 'thetaPR'),file=f)
    print('{:> 11.5g} {:> 11.5g} {:> 11.5g} {:> 11.5g} {:> 11.5g} {:> 8.4f} '.format(
        Iy,Iz,Iyz,Ipsi,Ieta, thetaPR),file=f)

    a='si' if esferico else 'no'
    b='si' if yzsonPR else 'no'
    print('\n{:>11s} {:>11s} {:^18s} {:^18s}'.format(
        'Wpsi','Weta','I esférico','yz principales'),file=f)
    print('{:> 11.5g} {:> 11.5g} {:^18s} {:^18s}'.format(
        Wpsi, Weta, a, b),file=f)
    print('     (W_psi está condicionado por yG,zG= {:10.4f} {:10.4f})'.format(
        ppsi[0], ppsi[1]), file=f)
    print('     (W_eta está condicionado por yG,zG= {:10.4f} {:10.4f})'.format(
        peta[0], peta[1]), file=f)
        
    print('\n{:>11s} {:>11s} {:>11s} {:>11s}'.format(
        'yE','zE','yE(G)','zE(G)'),file=f)
    print('{:> 11.5g} {:> 11.5g} {:> 11.5g} {:> 11.5g}'.format(
         E[0]+Gvec[0], E[1]+Gvec[1], E[0], E[1]),file=f)
    
    print('\n{:^11s} {:^11s}'.format('J','I_a'),file=f)
    print('{:> 11.5g} {:> 11.5g}'.format(J_secc, I_a),file=f)
    
    print('\nPara un momento Mx=1, y para el pto base elegido, el', file=f)
    print('promedio exacto de los desplazamientos ux de torsión es:', file=f)
    print('            G*ux_medio = {:^14.7e}'.format(uT_medio), file=f)

    print('\nTensor de Giros a Cortante C_ij, con i,j en y,z:', file=f)
    print(np.array(C_tensor), file=f)
    
    print('\nEn ejes principales de inercia, C_psi,eta es:', file=f)
    a=np.array([v_psi,v_eta])
    C_psieta= np.matmul(np.matmul(a,C_tensor),np.transpose(a))
    print(np.array(C_psieta), file=f)
    
    print('\nEn ejes principales propios, C_diag es:', file=f)
    a=np.array([v_psi,v_eta])
    C_diag,v= np.linalg.eigh(C_tensor)
    print('C_I={:11.5g} ;  C_II={:11.5g}'.format(C_diag[0],C_diag[1]), file=f)
    print('con vectores propios:', file=f)
    print('n_I = {:6.3f},  {:6.3f}'.format(v[0][0],v[0][1]),file=f)
    print('n_II= {:6.3f},  {:6.3f}'.format(v[1][0],v[1][1]),file=f)
    a=np.arctan(v[0][1]/v[0][0])*180/np.pi if v[0][0]!=0 else 90
    print('(ángulo respecto a yz= {:6.3f} grados)'.format(a), file=f)
    
    print('-'*60, file=f) 
    
    print('\nInformación orientativa sobre pandeo lateral', file=f)
    print( '--- (por favor consulte la documentación) ---\n',file=f)

    print('Se asume flexión según el eje fuerte (eta),', file=f)
    print('y que el material es un acero de construcción.\n', file=f)
    
    Llim= 16.02*np.sqrt(I_a/J_secc)
    texto='Orientativamente:\nPara longitud de la barra = {:11.5g} y superior, es\n'
    texto+='razonable despreciar el efecto del módulo de alabeo\n'
    print(texto.format(Llim),file=f)
    
    print('Valores orientativos de Mcr/E para algunas longitudes de barra,',file=f)
    print('siendo E~ 2.1e5 MPa (conviértase a las uds del problema):',file=f)
    
    Mcr= lambda L : Ieta*(np.pi/L)**2 *np.sqrt(I_a/Ieta +
        (L/np.pi)**2 *(J_secc/Ieta)/2.6)
     
    tabla=[]
    for L in np.linspace(5*anchoy,30*anchoy,6):
        tabla.append([L, Mcr(L)])
    
    print('\n    L    ', end='',file=f)
    for i in tabla:
        print('{:>11.4g}'.format(i[0]), end='', file=f)

    print('\n  Mcr/E  ', end='', file=f)
    for i in tabla:
        print('{:>11.4g}'.format(i[1]), end='', file=f)
    
    
    print('\n'+'-'*60, file=f)    
    
    print('\nNúcleo central de la sección: (coordenadas y,z en', file=f)
    print('ejes por G, ordenadas en sentido antihorario):', file=f)
    print('    y          z', file=f)
    for a in nucleo: print('{:> 11.5g} {:> 11.5g}'.format(a[0],a[1]), file=f)
    
    print('-'*60, file=f)

    print('\nSolucion en extremos de tramos para esfuerzos unitarios:\n', file=f)
    print('{:^4s} {:^11s} {:^11s} {:^11s}'.format(
        'i pto','ux(Vy=1)','ux(Vz=1)','ux(Mx=1)'), file=f)
    for ip,p in ptos.items():
        print('{:^4d} {:> 11.5g} {:> 11.5g} {:> 11.5g}'.format(
            ip, p.uVy, p.uVz, p.uT), file=f)
    
    print('\n{:^4s} {:^22s} {:^22s} {:>11s}'.format(
        '', 'q(Vy=1)', 'q(Vz=1)','q(Mx=1)'), file=f)
    print('{:^4s} {:^11s} {:^11s} {:^11s} {:^11s} {:^11s}'.format(
        'i_tr', 'i', 'j', 'i', 'j', 'j'), file=f)
    for it,tr in tramos.items():
        texto='{:^4d} {:> 11.5g} {:> 11.5g} {:> 11.5g} {:> 11.5g} {:> 11.5g}'
        print(texto.format(it,tr.qVy[0],tr.qVy[1],tr.qVz[0],tr.qVz[1],tr.qT[1]), file=f)

    print('-'*60, file=f)    
    
    print('\nInformación detallada para esfuerzos unitarios:\n',file=f)
    texto='{:^4s} {:^11s} {:^11s} {:^11s} {:^11s} {:^11s} {:^11s} {:^11s} {:^11s}'
    print(texto.format('i_tr','y','z','s','ux(Vy)', 'ux(Vz)','q(Vy)', 'q(Vz)','ux(Mx)'),file=f)
    texto='{:^4d} {:> 11.5g} {:> 11.5g} {:> 11.5g} {:> 11.5g} {:> 11.5g} {:> 11.5g} {:> 11.5g} {:> 11.5g}'
    for it,tr in tramos.items():
        for i in range(len(tr.tics)):
            print(texto.format(it, tr.tics[i][0], tr.tics[i][1], tr.tics[i][2],
                tr.uVy_plot[i], tr.uVz_plot[i], tr.qVy_plot[i], tr.qVz_plot[i], 
                tr.uT_plot[i]),file=f)

    texto='\n---------------------  fin del informe  ----------------------'
    print (texto, file=f)
    
    
    f.close()    
    
    texto='El informe se ha guardado en la ruta:\n'
    texto+= nfinforme
    messagebox.showinfo(message=texto, title='Confirmación', parent=v9)
    return()



#### Para despues de asignar ptos[] & tramos[]:
#### correr el motor de calculo y sacar ficha.

def haz_ficha(recarga=True, para=True):
    global n_f, nfcompleto, anchoy, anchoz
    global ptos,tramos,pto_base,Gvec,A_secc,eps_L,Iy,Iz,Iyz,det,esferico
    global thetaPR,psi,eta,Wpsi,Weta,ppsi,peta,E,J_secc,I_a,cerrado,theta_secc
    global uT_medio,nucleo,C_tensor
    global n_ptos,n_tramos,yzsonPR,Ipsi,Ieta,v_psi,v_eta
    
    if recarga:
        filas_a_ptos()
        filas_a_tr()
    motor_calculo()
    salida_terminal()
    
    boton_goM.state(['!disabled'])
    boton_goV.state(['!disabled'])
    boton_goT.state(['!disabled'])
    boton_txt.state(['!disabled'])
    try:  # guardar puede no existir si se viene desde +...
        boton_guardar.state(['!disabled'])
    except (TclError):
        pass
    
    # cerrar ventanas graficas y sacar la ventana de datos
    cierraplots()
    plt.figure('Datos')
    pintaPtosTramos(0,1,1,1)
    
    # ventana de ficha de la seccion / propiedades estaticas - (dibujo y texto)
    ficha = plt.figure('Ficha', figsize=(10., 6.9))
    gs = gridspec.GridSpec(1, 2, width_ratios=[5, 3]) 
    ax= plt.subplot(gs[0])

    # props estaticas- subplot de dibujo
    pintaPtosTramos(0,0,1,1)
    eps=(anchoy+anchoz)/58

    plt.arrow(Gvec[0],Gvec[1], anchoy*0.39, 0, 
        head_width=eps, linestyle='solid', linewidth=0.8,
        overhang=0.2, color='g')
    plt.text(Gvec[0]+anchoy*0.40,Gvec[1]+eps, 'y', fontsize=12, color='g')

    plt.arrow(Gvec[0],Gvec[1],0, anchoz*0.39,
        head_width=eps,  linestyle='solid', linewidth=0.8,
        overhang=0.2, color='g')
    plt.text(Gvec[0]+eps,Gvec[1]+anchoz*0.40, 'z', fontsize=12, color='g')

    plt.plot(Gvec[0],Gvec[1], marker='o', markersize=4, color='g')
    plt.text(Gvec[0]+eps,Gvec[1]+eps/2,'G',fontsize=12, fontweight='bold',
        color='g', horizontalalignment='center',verticalalignment='center')

    # dibujar el nucleo central
    y,z= [], []
    for p in nucleo:
        y.append(p[0]+Gvec[0])
        z.append(p[1]+Gvec[1])
    plt.plot(y,z, 'm--', linewidth=1)#, alpha=0.5)

    # dibujar el c.e.c.
    a,b=E[0]+Gvec[0],E[1]+Gvec[1]
    plt.plot(a,b, marker='o', markersize=4, color='grey')
    plt.text(a+eps,b-eps/2,'E',fontsize=12, fontweight='bold',
        color='grey', horizontalalignment='center',verticalalignment='center')
    
    # dibujar los ejes prles de inercia & el angulo
    rel=(Ieta/Ipsi)**0.25

    c=20*eps 	# el eje psi mas largo
    a, b= Gvec+ c*np.array(psi.es) ,  Gvec- c*np.array(psi.es)
    plt.plot([a[0],b[0]],[a[1],b[1]], color='b', linewidth=1)
    plt.text(a[0],a[1],r'$\xi$',color='b', fontsize=11)

    c=20*eps/rel	# eje eta mas corto
    a, b = Gvec + c*np.array(eta.es), Gvec- c*np.array(eta.es)
    plt.plot([a[0],b[0]],[a[1],b[1]], color='b', linewidth=1)
    plt.text(a[0],a[1]+eps,r'$\eta$',color='b', fontsize=11)

    if 0.1< abs(np.sin(thetaPR)) :   # arco solo si no enguarra
        arcJC(Gvec[0],Gvec[1], eps*7, np.pi/2, thetaPR, col='b',lw=1.)
        b=(np.pi+thetaPR)/2
        a=Gvec+ np.array([np.cos(b),np.sin(b)])*eps*8.
        plt.text(a[0],a[1], r'$\theta$',color='b', fontsize=11,
             horizontalalignment='center',verticalalignment='center')

    plt.title ('{:}'.format(n_f))
    plt.margins(0.06)
    plt.grid(True)
    plt.axis('equal')
    
    
    # props estaticas- subplot de texto
    ax= plt.subplot(gs[1])
    
    texto  = 'Baricentro G (ejes dados):\n'
    texto += '  Gy={:10.5f}\n  Gz={:10.5f}\n\n'.format(Gvec[0],Gvec[1])
    texto += 'Inercias (ejes por G):\n'
    texto += '  Iy={:9.4e}\n  Iz={:9.4e}\n  Iyz={:9.4e}\n\n'.format(
                                                                Iy,Iz,Iyz)

    texto += 'Inercias Principales (ejes por G):\n'
    texto += r'  I$\xi$'+'= {:9.4e}\n'.format(Ipsi)
    texto += r'  I$\eta$'+'= {:9.4e}\n'.format(Ieta)
    texto += r'  $\theta$'+'= {:10.5f}º\n'.format(thetaPR*180/np.pi)
    texto += r'  W$\xi$'+'= {:9.4e}\n'.format(Wpsi)
    texto += r'  W$\eta$'+'= {:9.4e}\n\n'.format(Weta)

    texto += 'Cte de Torsion '+'$J_T$={:10.5e}\n'.format(J_secc)
    texto += 'Mod de Alabeo  '+'$I_A$'+'={:10.5e}\n\n'.format(I_a)


    texto += 'C.E.Cortantes - ejes G (ejes O):\n'  
    texto += '  Ey={:10.5f} ({:10.5f})\n'.format(E[0],E[0]+Gvec[0])
    texto += '  Ez={:10.5f} ({:10.5f})\n\n'.format(E[1],E[1]+Gvec[1])

    texto += 'Area fisica= {:10.5f}\n'.format(A_secc)
    texto += r'$AC_Z$'+'= {:10.5f}'.format(1./ C_tensor[1][1])
    texto += ' (si procede)\n'
    texto += r'$AC_Y$'+'= {:10.5f}'.format(1./ C_tensor[0][0])
    texto += ' (si procede)\n'

    texto += '\nNucleo Central :  --------'

    plt.text(0.0, 0.02, texto, fontsize=12,bbox=dict(facecolor='#EEFFEE', 
        edgecolor='#EEF0FF', lw=8, boxstyle='round, pad=0.8'))
    
    plt.axis('off')
    #plt.tight_layout()
    if para: plt.show()



def a_guardar():
    global n_f, nfcompleto, ptos, tramos, pto_base
    
    nfcompleto = filedialog.asksaveasfilename(parent=v9,initialdir=dir_home,
        title='Guardar problema')
    if bool(nfcompleto):
        n_f = path.basename(nfcompleto)
        v9.title('ThinSecBeam 1.0 - '+n_f)
        
        if not pto_base: 
            for i in enumerate(ptos):
                pto_base=i[1]
                break  # toma el codigo del 1er pto que encuentra
        
        f=open(nfcompleto,'w')
        
        print('Identificación: ', nfcompleto,file=f)
        print('\nN_ptos   N_tramos   pto_base',file=f)
        print(' {:3d}  {:7d}  {:7d}'.format(len(ptos),len(tramos),pto_base),file=f)
        
        print('\nipto       y        z',file=f)
        for ip,p in ptos.items():
            print('{:3d}  {:8.3f}  {:8.3f}'.format(ip, p.y, p.z),file=f)
        
        print('\nitr   tipoT   i     j    e        R     +/-',file=f)
        
        for it,tr in tramos.items():
            if tr.tipoT == 0:
                print('{:3d} {:5d} {:5d} {:5d} {:7.3f}'.format(it, tr.tipoT,
                    tr.inciT[0], tr.inciT[1], tr.e),file=f)
            else:
                a='-' if tr.signo == -1 else '+' 
                print('{:3d} {:5d} {:5d} {:5d} {:7.3f}  {:7.3f}  {:}'.format(it,
                    tr.tipoT, tr.inciT[0], tr.inciT[1], tr.e, tr.R, a),file=f)
        print('###  Fin de datos suministrados  ###', file=f)
        f.close()
        a='El fichero de datos se ha guardado en:\n{:}'.format(nfcompleto)
        messagebox.showinfo(message=a, title='Confirmación',
             icon='info', parent=v9)
    f.close()
    
    
def a_salir():
    if messagebox.askokcancel(message='¿Quiere cerrar el programa? ',
                detail='Los datos no guardados se perderán.', default='cancel',
                icon='question', title='Confirmacion:',parent=v9) :
        for i in plt.get_fignums(): plt.close(i)
        v9.destroy()
        try: # si es linux que cierre el terminal tambien
            kill(getppid(), signal.SIGHUP)        
        except:
            exit()
        # con extension .pyw no saca terminal.

def cierraplots():
    for i in plt.get_fignums(): plt.close(i)


def masfilas():
    frame_ptos.destroy()
    frame_tramos.destroy()



################################
####  Funciones de dibujar  ####
################################

def arcJC(yC,zC,R,alfai,alfa,col='k',lw=2):
    step=abs(int(alfa/0.05))    # => step de unos 3 grados = 0.05rad
    y,z = [],[]
    for theta in np.linspace(alfai, alfai+alfa, step):
        y.append(yC+R*np.cos(theta))
        z.append(zC+R*np.sin(theta))
    plt.plot(y,z, color=col, linewidth=lw)


def yz3D_iguales(ax):
    ''' Finalmente no uso esta funcion porque va mal. Deja el trazado
    aun mas deformado. La solucion mejor (mientras lo arreglen) es
    redimensionar la ventana para que se vea proporcionado. 
    No obstante la idea era:
    
    Actualmente Matplotlib ya no deja opcion ax.set_aspect('equal')
    ni ax.axis('equal') en trazados 3D, por lo que la seccion sale deformada
    en su plano. Esta funcion es para llamarla justo antes de plt.show() para 
    poner ejes iguales en el plano de la seccion (es xy en matplotlib
    aunque nuestra nomenclatura sea yz). 
    
    Entrada= ax, unos ejes de matplotlib, por ej sacados de plt.gca().
    '''

    x_minmax = ax.get_xlim3d()
    x_ancho = abs(x_minmax[1] - x_minmax[0])
    x_medio = (x_minmax[1] + x_minmax[0])/2
    y_minmax = ax.get_ylim3d()
    y_ancho = abs(y_minmax[1] - y_minmax[0])
    y_medio = (y_minmax[1] + y_minmax[0])/2

    plot_ancho = max([x_ancho, y_ancho])/2
    #print('xy medios & semiancho=',x_medio,y_medio,plot_ancho)
    ax.set_xlim3d([x_medio - plot_ancho, x_medio + plot_ancho])
    ax.set_ylim3d([y_medio - plot_ancho, y_medio + plot_ancho])



def pintaPtosTramos(flagCl,flagN,flagP,flagT):
    # Flags: 
    # Cl= borrar la ventana. Implica que la ventana es ='Datos' !!
    #    y ademas que se llamara a filas_a_ptos & filas_a_tramos,
    #    que borran ptos{} & tramos{}, por lo que el usuario tiene
    #    que volver a correr motor_calculo (dar a Calcula) para 
    #    sacar resultados nuevos (por ej un Informe).
    #    A pensar: ¿correr motor_calculo cada vez que se vaya a hacer algo?
    #              ¿deshabilitar botones de dibujo/comprobacion tras Calcula?
    # N= que ponga etiqueta (numeros) a los ptos & tramos, y flecha a los tr
    # P= pintar los puntos
    # T= pintar los tramos
    global ptos,tramos, anchoy, anchoz

    if flagCl : 
        plt.figure('Datos')
        plt.clf()
    
    if flagP:
        if flagCl: filas_a_ptos()
        # necesito anchos yz provisionales:
        anchoy, anchoz = 0., 0.
        for pi in ptos.values():
            for pj in ptos.values():
                d= abs(pi.y - pj.y)
                if d > anchoy: anchoy=d
                d= abs(pi.z - pj.z)
                if d > anchoz: anchoz=d
        dy, dz= anchoy/99, anchoz/99
        eps=(anchoy+anchoz)/58
        
        for ip,p in ptos.items():
            plt.plot(p.y, p.z, marker='o', markersize=5, color='green')
            if flagN:
                plt.text(p.y , p.z+dz, '{:2d}'.format(ip),
                color='g', fontweight='bold', fontsize=10)
    if flagT:
        if flagCl: 
            filas_a_tr()
            boton_guardar.state(['disabled'])
            boton_txt.state(['disabled'])
            boton_goM.state(['disabled'])
            boton_goV.state(['disabled'])
            boton_goT.state(['disabled'])

        # necesito ya centros y alfas (aunque se recalculen en motor):
        basico_tr()
        
        for it,tr in tramos.items():
            i,j= tr.inciT[0], tr.inciT[1]
            yi,zi= ptos[i].y, ptos[i].z
            yj,zj= ptos[j].y, ptos[j].z
            if tr.tipoT == 0:  # es recto
                plt.plot([yi,yj],[zi,zj], color='k', linewidth=2)
                es=np.array([yj-yi, zj-zi])  # (no es unitario)
                if flagN:
                    p_medio=np.array([yi,zi]) + es/2.1
                    alfa=np.arctan2((zj-zi),(yj-yi))
                    alfatext=alfa
                    if  alfa < -1.396 :
                        alfatext += np.pi
                    elif alfa > 1.571:
                        alfatext -= np.pi
                    plt.text(p_medio[0], p_medio[1],'{:2d}'.format(it), 
                        fontsize=10, rotation=alfatext*180/np.pi,
                        horizontalalignment='right',
                        verticalalignment='bottom',
                        bbox=dict(facecolor='g', edgecolor='white', 
                        alpha=0.2, boxstyle='round'))
                    # si se numera, habra flecha tambien:
                    p_medio=np.array([yi,zi]) + es/1.9
                    plt.arrow(p_medio[0],p_medio[1],
                        dy*np.cos(tr.alfa),dy*np.sin(tr.alfa),
                        head_width=2.5*dy,linestyle='solid', 
                        linewidth=1.2, overhang=0.6, color='g')
                    # & habra espesores tambien:
                    p_medio=np.array([yi,zi]) + es/2.4
                    plt.text(p_medio[0]-eps, p_medio[1]-eps,
                        '{:5.3g}'.format(tr.e), fontsize=10,style='italic',
                        rotation=alfatext, horizontalalignment='center',
                        verticalalignment='center')
                
                
            else:  # es circular
                #ax=plt.gca()
                #a,b=(tr.alfai*180/np.pi,(tr.alfai+tr.alfa )*180/np.pi)
                #if tr.alfa<0.: a,b=b,a
                arcJC(tr.yC, tr.zC, tr.R, tr.alfai, tr.alfa)

                if flagN:
                    a= tr.alfai + tr.alfa*0.48
                    p_medio=[tr.yC+tr.R*np.cos(a), tr.zC+tr.R*np.sin(a)]
                    plt.text(p_medio[0], p_medio[1],'{:2d}'.format(it), 
                        fontsize=10, 
                        horizontalalignment='right',
                        verticalalignment='bottom',
                        bbox=dict(facecolor='g', edgecolor='white', 
                        alpha=0.2, boxstyle='round'))
                    # idem, si se numera, habra flecha tambien:
                    a= tr.alfai + tr.alfa*0.59
                    p_ini=np.array([tr.yC+tr.R*np.cos(a), tr.zC+tr.R*np.sin(a)])
                    a= tr.alfai + tr.alfa*0.591
                    p_fin=np.array([tr.yC+tr.R*np.cos(a), tr.zC+tr.R*np.sin(a)])
                    es=p_fin-p_ini
                    plt.arrow(p_ini[0],p_ini[1], es[0], es[1],
                        head_width=2.5*dy,linestyle='solid', 
                        length_includes_head=True,
                        linewidth=1.2, overhang=0.6, color='g')
                    # y habra espesor:
                    a= tr.alfai + tr.alfa*0.67
                    p_medio=[tr.yC+(tr.R-eps)*np.cos(a), tr.zC+(tr.R-eps)*np.sin(a)]
                    a *= 180/np.pi
                    alfatext= a-90 if a<180 else a+90 
                    plt.text(p_medio[0], p_medio[1],'{:5.3g}'.format(tr.e), 
                        fontsize=10,style='italic',
                        rotation= alfatext,
                        horizontalalignment='center',
                        verticalalignment='center')
                    

    if flagN:
        plt.xlabel('y',fontweight='bold', fontsize=10)
        plt.ylabel('z',fontweight='bold', fontsize=10)
    plt.grid(True)
    plt.axis('equal')
    plt.margins(0.1)
    if flagT: 
        try:  # guardar puede no existir si se viene desde +...
            boton_guardar.state(['!disabled'])
        except (TclError):
            pass
    if flagCl: plt.show() # solo viene con True para la fig de Datos
                          # Es un control chapuza pero suficiente



# Ventanas de sxx - flexion en 2d & 3D

def haz_M(para=True):
    global ptos, tramos, anchoy, anchoz, Gvec, etiq2D,etiq3D
    
    My_dato,Mz_dato = float(entry_My.get()),float(entry_Mz.get())
    Nx_dato = float(entry_Nx.get())
    eps=(anchoy+anchoz)/58
    
    # ventana de sxx en 2D

    plt.figure('Tensiones de flexion -2D')
    plt.title('{:}'.format(n_f))
    ax=plt.gca()    
    
    pintaPtosTramos(0,0,1,1)
    
    # pinta una flecha en sentido del flector
    alfa=np.arctan2(-Mz_dato,My_dato)
    plt.arrow(Gvec[0],Gvec[1],
            7*eps*np.cos(alfa), 7*eps*np.sin(alfa), 
            head_width= eps, linestyle='solid', linewidth=1.2,
            overhang=0.6, color='r')
    alfa -= 0.3
    if etiq2D:
        plt.text(Gvec[0]+ 8.5*eps*np.cos(alfa), 
                Gvec[1]+ 8*eps*np.sin(alfa),
                'M', color='r', fontsize=14, fontweight='normal',
                horizontalalignment='center',verticalalignment='center')    
    
    # dibujar los ejes prles de inercia
    rel=(Ieta/Ipsi)**0.25

    c=20*eps 	# el eje psi mas largo
    a, b= Gvec+ c*np.array(psi.es) ,  Gvec- c*np.array(psi.es)
    plt.plot([a[0],b[0]],[a[1],b[1]], color='m', linewidth=2,
        linestyle=':')

    c=20*eps/rel	# eje eta mas corto
    a, b = Gvec + c*np.array(eta.es), Gvec- c*np.array(eta.es)
    plt.plot([a[0],b[0]],[a[1],b[1]], color='m', linewidth=2,
        linestyle=':')
    
    # calcula y pinta la LN:

    aLN= (Mz_dato*Iy-My_dato*Iyz)/det
    bLN= (My_dato*Iz-Mz_dato*Iyz)/det
    cLN= Nx_dato/A_secc
    
    k= np.sqrt(aLN*aLN + bLN*bLN)
    es=np.array([-bLN,aLN]) /k
    en=np.array([aLN,bLN]) /k
    if   abs(bLN) < 1.e-8*abs(aLN):
        p=np.array([-cLN/aLN,0.])
    elif abs(aLN) < 1.e-8*abs(bLN):
        p=np.array([0., -cLN/bLN])
    else:
        p=np.array([-cLN/aLN,-cLN/bLN])/2
      # ese calculo es para ejes por G. Corregimos a ejes dados
    p += np.array(Gvec)

    LN=Recta(p,es,en)
    
        # saber el lado de traccion:
    k=LN.p + LN.en
    tracc_LN_en=np.sign( sxx_yz(Nx_dato,My_dato,Mz_dato, 
            k[0]-Gvec[0], k[1]-Gvec[1]))
    
    trazar_LN=True if abs(LN.dist(Gvec)) < anchoy else False
    if trazar_LN:
        r1=Recta(Gvec, LN.en, [0.,0.])  # no necesito r1.en[]
        pG= LN.p+ LN.intersec(r1)[0]*LN.es # pto de LN cercano a G
        a,b = pG + 30*eps*LN.es,   pG - 24*eps*LN.es
        # la LN:
        plt.plot([a[0],b[0]],[a[1],b[1]], color='b', linewidth=1.2)
        
        c= b + 1.3*eps* LN.en* tracc_LN_en
        if LN.es[0] != 0.:
            alfa= np.arctan(LN.es[1]/LN.es[0])
        else:
            alfa= np.sign(LN.es[1])*np.pi/2.

            # acota tr-c:
        plt.text(c[0], c[1], 'Tracc', color='b', fontsize=11,
                horizontalalignment='center',verticalalignment='center',
                rotation=alfa*180/np.pi)

        c= b - 1.3*eps*LN.en*tracc_LN_en
        plt.text(c[0], c[1], 'Compr', color='b', fontsize=11,
                horizontalalignment='center',verticalalignment='center',
                rotation=alfa*180/np.pi)
    else:
        p=[Gvec[0]-anchoy*0.6 , Gvec[1]-anchoz*0.7]
        texto='LN alejada.\nNo se traza.' 
        plt.text(p[0], p[1], texto, color='b', fontsize=10)

    # distancias menor y mayor a la LN de los ptos (con signo en)
    dmin, dmax = 1.0e99, -1.0e99
    pmin, pmax = [0.,0.],[0.,0.]
    for tr in tramos.values():
        for i in range(len(tr.tics)):
            a= np.array([tr.tics[i][0],tr.tics[i][1]])
            d=LN.dist(a)
            if d < dmin:  dmin, pmin = d, a
            if d > dmax:  dmax, pmax = d, a
    
        # lineas paralelas a LN por extremos, y la recta tope
    r1= Recta(Gvec+24*eps*LN.es, LN.en, -LN.es)    # recta tope
    rmin=Recta(pmin, LN.es, LN.en)
    rmax=Recta(pmax, LN.es, LN.en)  # // a LN por pto min & max
    qmin= pmin+ rmin.es* rmin.intersec(r1)[0]
    qmax= pmax+ rmax.es* rmax.intersec(r1)[0]  # ptos homologos de max-min sobre tope
    
        # trazado de las lineas auxiliares
    plt.plot( [[pmin[0]],[qmin[0]]] , [[pmin[1]],[qmin[1]]], '--g', lw=1)
    plt.plot( [[pmax[0]],[qmax[0]]] , [[pmax[1]],[qmax[1]]], '--g', lw=1)
    plt.plot( [[qmin[0]],[qmax[0]]] , [[qmin[1]],[qmax[1]]], 'g', lw=1) # recta tope
    
        # miro si la LN pasa entre medias (hay tr & c en la seccion)
    bmin= sxx_yz(Nx_dato,My_dato,Mz_dato,qmin[0]-Gvec[0],qmin[1]-Gvec[1])
    bmax= sxx_yz(Nx_dato,My_dato,Mz_dato,qmax[0]-Gvec[0],qmax[1]-Gvec[1])
    cambia=False if np.sign(bmin) == np.sign(bmax) else True
    
        # redefino r1 (r tope:)
    d=np.sqrt((qmax[0]-qmin[0])**2 + (qmax[1]-qmin[1])**2 ) # es L de r tope
    es= np.array([ qmax[0]-qmin[0], qmax[1]-qmin[1] ])/d
    r1=Recta([qmin[0],qmin[1]], es , LN.es) 

        # escala para trazar
    sxx_tipo= (abs(My_dato)+abs(Mz_dato))*anchoy/(4*Iz) + abs(Nx_dato)/A_secc
    sxx_scal= 9*eps/sxx_tipo 
    
        # acoto extremos
    pbmin=qmin+ abs(bmin)*sxx_scal*r1.en 
    aju=eps*(r1.en-r1.es)
    plt.text( pbmin[0]+aju[0], pbmin[1]+aju[1], 
        '{:.4}'.format(bmin),
        color='g', fontsize=11, fontweight='normal',
        horizontalalignment='center',verticalalignment='center') 

    pbmax=qmax+ abs(bmax)*sxx_scal*r1.en
    aju=eps*(r1.en+r1.es)
    plt.text( pbmax[0]+aju[0], pbmax[1]+aju[1], 
        '{:.4}'.format(bmax),
        color='g', fontsize=11, fontweight='normal',
        horizontalalignment='center',verticalalignment='center') 
    
    if cambia: 
            # trazar linea quebrada
        pcorte=qmin + r1.es * r1.intersec(LN)[0]
        plt.plot([pbmin[0], pcorte[0], pbmax[0]] ,
                 [pbmin[1], pcorte[1], pbmax[1]] , color='g', linewidth=1)
    else:   # trazar linea seguida
        plt.plot([pbmin[0], pbmax[0]] ,
                 [pbmin[1], pbmax[1]] , color='g', linewidth=1)
        
    for s in np.linspace(0.,d, 12):  # d = L de rtope (r1)
        p0=qmin+ s*r1.es
        b=sxx_yz(Nx_dato,My_dato,Mz_dato,p0[0]-Gvec[0],p0[1]-Gvec[1])
        p1=p0+abs(b)*sxx_scal *r1.en
        if b<0.: p0,p1 = p1,p0
        # que no haga flecha cerca del corte (si hay):
        if cambia:
            s0=np.sqrt( (pcorte[0]- qmin[0])**2 + (pcorte[1]- qmin[1])**2 )
        if not cambia or abs(s0-s)> d/10:
            plt.arrow(p0[0],p0[1], p1[0]-p0[0], p1[1]-p0[1],
                head_width= 0.8*eps, linestyle='solid', linewidth=0.8,
                overhang=0.2, color='g',length_includes_head=True)

    plt.margins(0.05,0.05)
    plt.grid(True)
    plt.axis('equal')


    # figura 3D de sxx

    vent_sxx3D=plt.figure('Tensiones de flexion -3D')
    ax=vent_sxx3D.gca(projection='3d')
    #ax.set_aspect('equal') sustituido por llamada a yz3D_iguales()
    plt.title('{:}'.format(n_f))

    # dibujar flecha del Momento
    alfa=np.arctan2(-Mz_dato,My_dato)
    a= Gvec+ np.array([np.cos(alfa), np.sin(alfa)])*anchoy/2.2
    ax.plot3D([Gvec[0], a[0]], [Gvec[1], a[1]], [0.,0.],
              linewidth=2, color='r')
    b=a+ np.array([np.cos(alfa+2.8), np.sin(alfa+2.8)])*anchoy/13
    ax.plot3D([a[0], b[0]], [a[1],b[1]], [0.,0.],
              linewidth=2, color='r')
    b=a+ np.array([np.cos(alfa-2.8), np.sin(alfa-2.8)])*anchoy/13
    ax.plot3D([a[0], b[0]], [a[1],b[1]], [0.,0.],
              linewidth=2, color='r')    
    
    if etiq3D:
        ax.text(a[0]+ 2*eps*np.cos(alfa+1.5), 
                 a[1]+ 2*eps*np.sin(alfa+1.5), 0.,
                'M', color='r', fontsize=13, fontweight='bold',
                horizontalalignment='center',verticalalignment='center')    
    
    # dibujar los ejes prles de inercia
    rel=(Ieta/Ipsi)**0.25

    c=16*eps 	# el eje psi mas largo
    a, b= Gvec+ c*np.array(psi.es) ,  Gvec- c*np.array(psi.es)
    ax.plot([a[0],b[0]],[a[1],b[1]], color='m', linewidth=0.8,
        linestyle='-')

    c=16*eps/rel	# eje eta mas corto
    a, b = Gvec + c*np.array(eta.es), Gvec- c*np.array(eta.es)
    ax.plot([a[0],b[0]],[a[1],b[1]], color='m', linewidth=0.8,
        linestyle='-')


        # traza la seccion & la sxx:
    for tr in tramos.values():
        ax.plot3D(tr.tics[:,0],tr.tics[:,1], 0, linewidth=2, color='k')
        a=[]
        for p in tr.tics:
            yG,zG= p[0]-Gvec[0],  p[1]-Gvec[1]
            a.append(sxx_yz(Nx_dato,My_dato,Mz_dato, yG,zG))
        ax.plot3D(tr.tics[:,0],tr.tics[:,1], a, linewidth=2, color='g')
        # y una linea aux a mitad de tr curvo (sin etiq):
        if tr.tipoT ==1:
            y,z= tr.yz_s(tr.L/2)-Gvec
            sxx=sxx_yz(Nx_dato,My_dato,Mz_dato, y-Gvec[0], z-Gvec[1])
            ax.plot3D([y, y],
                      [z, z],
                      [0. , sxx],
                        linewidth=1, linestyle='--', color='g')

    # traza lineas aux por ptos[]:
    for p in ptos.values():
        sxx=sxx_yz(Nx_dato,My_dato,Mz_dato, p.yG,p.zG)
        ax.plot3D([p.y, p.y],
                  [p.z, p.z],
                  [ 0. ,  sxx],
                    linewidth=1, linestyle='--', color='g')
        # si procede acota:
        if etiq3D:
            ax.text(p.y+ eps, p.z+eps, sxx,
                '{:9.3e}'.format(sxx), 
                color='grey', fontsize=10, fontweight='normal',
                horizontalalignment='center',verticalalignment='center'  )            

    # traza LN:
    if trazar_LN:
        a,b = pG + 24*eps*LN.es,   pG - 24*eps*LN.es
        ax.plot3D([a[0],b[0]], [a[1],b[1]], [0,0], color='b', linewidth=1.2)
        if etiq3D:
            ax.text(b[0]+ eps, b[1]+eps, 0,'LN',color='b', fontsize=11)
            am=abs(sxx_tipo)/7
            ax.text(a[0], a[1], am*0.8, ' Tracc', 'x',
                fontsize=11, color='b', fontweight='normal') 
            ax.text(a[0], a[1], -am, ' Compr', 'x', 
                fontsize=11, color='b', fontweight='normal') 
            am *= 0.8
            ax.plot( [a[0],a[0]], [a[1],a[1]], [-am, am], '--',
                color='b', linewidth=1., marker='.')

    if para: 
        #yz3D_iguales(ax)
        plt.show()
    
    
# Ventanas de uT (en 3D) y de qT (en 2D)

def haz_T(para=True):
    global ptos, tramos, anchoy, anchoz, Gvec, theta_secc
    global uT_medio, uT_despl, etiq2D,etiq3D, peque3D
    
    T_dato = float(entry_T.get())
    eps=(anchoy+anchoz)/58
    
    # mira si es abierta 
    abierto =True if len(ptos) == len(tramos)+1 else False
    
    # asigna una tension sT_tipo (en cerrado solo mira los de flujo)
    sT_tipo , L_tot = 0.,0.
    for tr in tramos.values():
        L_tot += tr.L
        if abierto:
            a= abs(theta_secc* T_dato* tr.e)
            if a > sT_tipo: sT_tipo=a
        else:
            i,j = tr.inciT
            a= abs(tr.qT[1]* T_dato / tr.e)
            if a > sT_tipo: sT_tipo=a
    scal= 3*eps/ sT_tipo


    # ventana de tensiones de torsion (tangenciales, en 2D)
    
    plt.figure('Tensiones de torsion')
    plt.title('{:}'.format(n_f))
    ax=plt.gca()
    
    pintaPtosTramos(0,0,1,1)
    
        # pinta una flecha circular en sentido del torsor
    a=r'$\circlearrowleft$' if T_dato>0. else r'$\circlearrowright$'
    ax.plot([Gvec[0]],[Gvec[1]],marker=a,markersize=35.,
            alpha=0.7, color='r')
    
    # Pinta tensiones de ambos tipos:
    for tr in tramos.values():
        s_vaiven, s_flujo = theta_secc* T_dato*tr.e ,  tr.qT[1]* T_dato/tr.e
        i,j=tr.inciT
        pi=np.array([ptos[i].y, ptos[i].z])
        pj=np.array([ptos[j].y, ptos[j].z])
        acotar_etc= True if tr.L*25 > L_tot else False 
        
        if tr.tipoT ==0:  # tramo recto
            es=np.array([np.cos(tr.alfa),np.sin(tr.alfa)])
            en=np.array([ es[1], -es[0] ])
            
            if acotar_etc:
                p_texto= pi+ 0.6*tr.L*es
                alfatext=tr.alfa*180/np.pi
                if np.cos(tr.alfa-0.17)<0: alfatext += 180 

                # pinta el texto segun sea tramo abierto o cerrado: 
                if etiq2D:
                    if abierto or abs(s_flujo/sT_tipo)< 1.e-5:
                        p_texto += eps*en 
                        plt.text (p_texto[0], p_texto[1],
                            r'$\pm$'+'{:8.3e}'.format(abs(s_vaiven)),
                            fontsize=11, fontweight='normal', color='m',
                            horizontalalignment='center',
                            verticalalignment='center',rotation= alfatext)
                    else:
                        texto='{:8.3e}\n'+r'$\pm$'+'{:8.3e}'
                        plt.text (p_texto[0], p_texto[1],
                            texto.format(s_flujo, abs(s_vaiven)),
                            fontsize=11, fontweight='normal', color='m',
                            horizontalalignment='center',
                            verticalalignment='center',rotation= alfatext)
            
            # dibuja flecha y lineas (ctes), sin enguarrar:
            if not abierto and abs(s_flujo/sT_tipo) > 1.e-5:
                if acotar_etc:
                    p_texto= pi + 0.3*tr.L*es
                    a= 0.2*eps* es* np.sign(s_flujo)
                    plt.arrow(p_texto[0], p_texto[1], a[0], a[1], 
                        head_width=eps, length_includes_head=True, 
                        overhang=0.2, head_starts_at_zero=True, color='m')

                yplot, zplot =[],[]

                pplot= pi+ s_flujo* scal* en
                yplot.append(pplot[0])
                zplot.append(pplot[1])
                
                pplot= pj+ s_flujo* scal* en
                yplot.append(pplot[0])
                zplot.append(pplot[1])

                plt.plot(yplot, zplot, '-m', linewidth=2)                

                if acotar_etc:
                    plt.plot([pi[0], yplot[0]],
                             [pi[1], zplot[0]], '--m')
                    plt.plot([pj[0], yplot[-1]],
                             [pj[1], zplot[-1]], '--m')

            if acotar_etc:   # dibuja vaiven, sea tr abierto o cerrado
                sig=np.sign(T_dato)
                b='right' if sig <0. else 'left'
                
                p_texto= pi+ eps*es+ eps*en/3.
                a= 0.2*eps* es*sig
                plt.arrow(p_texto[0], p_texto[1], a[0], a[1],
                    shape= b,
                    head_width=0.7*eps, head_length=eps, 
                    length_includes_head=True, 
                    overhang=0., head_starts_at_zero=True, color='m')
                                    
                p_texto= pi+ eps*es- eps*en/3.
                a= -0.2*eps*es*sig
                plt.arrow(p_texto[0], p_texto[1], a[0], a[1],
                    shape= b,
                    head_width=0.7*eps, head_length=eps,
                    length_includes_head=True, 
                    overhang=0., head_starts_at_zero=True, color='m')

        elif tr.tipoT==1:  # tramo curvo
            
            R,L,alfa,alfai,signo= tr.R, tr.L, tr.alfa, tr.alfai, tr.signo
            yC,zC= tr.yC, tr.zC
            
            if acotar_etc:
                p_texto= tr.yz_s(0.6*L) -Gvec
                alfatext= (alfai+0.6*alfa) -np.pi/2
                if np.cos(alfatext-0.17)<0: alfatext += np.pi
                en=np.array([np.cos(alfatext), np.sin(alfatext)])*signo
                es=np.array([-en[1],en[0]])
                alfatext *= 180/np.pi
                
                # pinta el texto segun sea tramo abierto o cerrado:
                if etiq2D:
                    if abierto or abs(s_flujo/sT_tipo)< 1.e-5:
                        p_texto += eps*en 
                        plt.text (p_texto[0], p_texto[1],
                            r'$\pm$'+'{:8.3e}'.format(abs(s_vaiven)),
                            fontsize=11, fontweight='normal', color='m',
                            horizontalalignment='center',
                            verticalalignment='center',rotation= alfatext)
                    else:
                        texto='{:8.3e}\n'+r'$\pm$'+'{:8.3e}'
                        plt.text (p_texto[0], p_texto[1],
                        texto.format(s_flujo, abs(s_vaiven)),
                            fontsize=11, fontweight='normal', color='m',
                            horizontalalignment='center',
                            verticalalignment='center',rotation= alfatext)

            # dibuja flecha y lineas (ctes), sin enguarrar:
            if not abierto and abs(s_flujo/sT_tipo)> 1.e-5:
                if acotar_etc:
                    p_texto=tr.yz_s(0.3*L)-Gvec
                    es=tr.es_s(0.3*L) 
                    a=0.2*eps*np.sign(s_flujo)*es
                    plt.arrow(p_texto[0], p_texto[1], a[0], a[1], 
                        head_width=eps, length_includes_head=True, 
                        overhang=0.2, head_starts_at_zero=True, color='m')
                        
                Rpinta=R + signo* s_flujo* scal
                arcJC(yC,zC,Rpinta,alfai,alfa,col='m',lw=2)
                if acotar_etc:  # lineas guia en los extremos
                    a=[yC,zC]+Rpinta*np.array([np.cos(alfai),np.sin(alfai)])
                    plt.plot([pi[0],a[0]] , [pi[1],a[1]] ,'--m')
                    a=[yC,zC]+Rpinta*np.array([np.cos(alfai+alfa),np.sin(alfai+alfa)])
                    plt.plot([pj[0],a[0]] , [pj[1],a[1]] ,'--m')
                    
            if acotar_etc:   # dibuja vaiven, sea tr abierto o no
                sig=np.sign(T_dato)
                b='right' if sig <0. else 'left'
                en, es = tr.en_s(1.5*eps), tr.es_s(1.5*eps)                

                p_texto=tr.yz_s(1.5*eps)+ eps*en/3.- Gvec
                a= +0.2*eps* es *sig
                plt.arrow(p_texto[0], p_texto[1], a[0], a[1],
                    shape= b,
                    head_width=0.7*eps, head_length=eps, 
                    length_includes_head=True, 
                    overhang=0., head_starts_at_zero=True, color='m')
                
                p_texto=tr.yz_s(1.5*eps)- eps*en/3.- Gvec
                a= -a
                plt.arrow(p_texto[0], p_texto[1], a[0], a[1],
                    shape= b,
                    head_width=0.7*eps, head_length=eps,
                    length_includes_head=True, 
                    overhang=0., head_starts_at_zero=True, color='m')





    # ventana de desplazamientos de alabeo por torsion (en 3D)
    
    vent_uV=plt.figure('Alabeo por Torsión')
    ax=vent_uV.gca(projection='3d')
    #ax.set_aspect('equal') sustituido por llamada a yz3D_iguales()
    plt.title('{:}'.format(n_f))
    
    if uT_despl=='': uT_despl=2  # se inicializa 'centrado al promedio'

    uTmax=0.
    for p in ptos.values():
        if abs(p.uT) > abs(uTmax) : uTmax= p.uT  #  max de uT
    uTmax *= abs(T_dato)
    
        # avisar si el alabeo parece nulo
    if abs(uTmax*T_dato)/(sT_tipo*25*eps) < 5.e-5:
        ax.text(Gvec[0], Gvec[1]+12*eps, uTmax/3,
            '¿Alabeo nulo?', 
            color='b', fontsize=34, fontweight='bold', alpha=0.5,
            horizontalalignment='center',verticalalignment='center'  )
    
        # pintar flecha en direccion de Mx
    sig= np.sign(T_dato)
    a= sig*abs(uTmax)/2.   # a= L del vector a escala de uT
    if peque3D: a*= 4
    b= a/9.          # L cabeza de la flecha
    c= eps/3.        # semiancho cabeza de flecha 
    ax.plot([Gvec[0], Gvec[0]], 
            [Gvec[1], Gvec[1]], 
            [0.,      a],
            linewidth=2, color='r')
    ax.plot([Gvec[0], Gvec[0]], 
            [Gvec[1], Gvec[1]-c], 
            [a,       a-b],
            linewidth=3, color='r')
    ax.plot([Gvec[0], Gvec[0]], 
            [Gvec[1], Gvec[1]+c], 
            [a,       a-b],
            linewidth=3, color='r')
    if etiq3D:
        ax.text(Gvec[0]+ eps, Gvec[1]+eps, a+b/2,
            'T', color='r', fontsize=14, fontweight='normal',
            horizontalalignment='center',verticalalignment='center'  )

        # asignacion del corrimiento para trazar uT:
    if uT_despl == 2:  # centrado al promedio
        umed= T_dato*uT_medio
        d=-umed
    elif uT_despl == 1:  # afuera de la seccion
        d=-uTmax
    else:   # sera uV_despl=0, tal como se calcula
        d=0.    

        # dibuja ejes principales
    rel=(Ieta/Ipsi)**0.25

    c=12*eps 	# el eje psi mas largo
    a, b= Gvec+ c*np.array(psi.es) ,  Gvec- c*np.array(psi.es)
    ax.plot([a[0],b[0]],[a[1],b[1]], color='b', linewidth=1)

    c=12*eps/rel	# eje eta mas corto
    a, b = Gvec + c*np.array(eta.es), Gvec- c*np.array(eta.es)
    ax.plot([a[0],b[0]],[a[1],b[1]], color='b', linewidth=1)
    
        # trazado basico de seccion & los uV
    for it, tr in tramos.items():
        ax.plot(tr.tics[:,0],tr.tics[:,1], color='k')
        uT_plot= T_dato*np.array(tr.uT_plot)
        ax.plot3D(tr.tics[:,0],tr.tics[:,1],uT_plot[:]+d,
            linewidth=2, color='m')

    if peque3D:
        # un puntin perpendicular, que quede ux mas pegado a la seccion
        a=[4*uTmax, 5*uTmax, 4*uTmax]
        ax.plot3D([Gvec[0]], [Gvec[1]], [a[uT_despl]] )
        ax.plot3D([Gvec[0]], [Gvec[1]], [-a[uT_despl]] )
        
        # en ptos: lineas guia y etiquetar ptos y despl
    for ip,p in ptos.items():
        a= T_dato*p.uT
        ax.plot([p.y, p.y],[p.z, p.z],[0.,a+d],'--',lw=0.8,color ='m')
        if etiq3D:
            ax.text(p.y, p.z, a+d , '{:.4}'.format(a), 'x', 
                fontsize=11, color='m', fontweight='normal')
            ax.text(p.y, p.z, 0, ip, 'x', fontsize=10, color='g', 
                fontweight='bold') 

    if para: 
        #yz3D_iguales(ax)
        plt.show()




# Ventanas de qV (en 2D) y de uV (en 3D)

def haz_V(para=True):
    global ptos, tramos, anchoy, anchoz, Gvec, C_tensor, C0y, C0z, E
    global uVy_medio,uVz_medio, uV_despl, etiq2D,etiq3D,peque3D
    
    Vy_dato, Vz_dato = float(entry_Vy.get()), float(entry_Vz.get())


    # Ventana de flujos de cortante (2D)
    
    plt.figure('Flujos de tension (qV)')
    plt.title('{:}'.format(n_f))

    pintaPtosTramos(0,0,1,1)
    eps=(anchoy+anchoz)/58

        # dibuja E (para que se sepa donde estaria V):
    plt.plot(Gvec[0]+E[0],Gvec[1]+E[1], marker='o', markersize=4, color='grey')
    if etiq2D:
        plt.text(Gvec[0]+E[0]+eps,Gvec[1]+E[1]+eps/2,'E',fontsize=11, 
            fontweight='bold', color='grey', 
            horizontalalignment='center',verticalalignment='center')        
    
        # dibuja vector V (por E o mas cerca de G)
    alfa=np.arctan2(Vz_dato,Vy_dato)
    d=np.sqrt(E[0]**2+E[1]**2)
    (a,b)= (E[0],E[1]) if d<anchoy/2 else (E[0]*0.6,E[1]*0.6)
    a +=Gvec[0]
    b +=Gvec[1] 
    plt.arrow(a,b,
            7*eps*np.cos(alfa), 7*eps*np.sin(alfa), 
            head_width= eps, linestyle='solid', linewidth=1.2,
            overhang=0.6, color='r')
    alfa -= 0.2
    if etiq2D:
        plt.text(a+ 8*eps*np.cos(alfa), 
                b+ 8*eps*np.sin(alfa),
                'V', color='r', fontsize=14, fontweight='normal',
                horizontalalignment='center',verticalalignment='center')


        # Escala para dibujar= q_scal. Long total= L_tot:
    qmax, L_tot = 0., 0.
    for it, tr in tramos.items(): 
        L_tot += tr.L
        k=abs(max([tr.qVy[0],tr.qVy[1],tr.qVz[0],tr.qVz[1]], key=abs))
        if qmax < k : qmax = k  # lo reutilizo luego a valor de tr, no perdura
    q_scal=(anchoy+anchoz)/(11*qmax*(abs(Vy_dato)+abs(Vz_dato)))
    
    for it,tr in tramos.items():
        acotar_etc= True if tr.L*28 > L_tot else False 
        if tr.tipoT == 0:
            i,j= tr.inciT[0], tr.inciT[1]
            yi,zi=ptos[i].y , ptos[i].z
            yj,zj=ptos[j].y , ptos[j].z
            incy,incz= yj-yi , zj-zi
            es=np.array([incy, incz]) /tr.L
            en=np.array([ es[1] , -es[0] ]) # sera + a la dcha del tr
            y,z = [], []
            
            valoresq=[]
            for i in range(len(tr.tics)):  # traza la funcion
                v= np.array([tr.tics[i][0], tr.tics[i][1]])
                a=tr.qVy_plot[i]*Vy_dato + tr.qVz_plot[i]*Vz_dato
                valoresq.append(a)
                w = v+ en*a*q_scal
                y.append(w[0])
                z.append(w[1])
                if acotar_etc: # lineas punteadas
                    if i==0 and abs(valoresq[0]) > abs(qmax)/1.e7:
                        # chapucilla: ya sea qmax de otro tr, o absoluto
                        # se trata solo de evitar trazar en extremos de tr
                        plt.plot([v[0],w[0]] , [v[1],w[1]], color='g',
                            linestyle='--', linewidth=1)
                    if i== len(tr.tics)-1 and abs(valoresq[i]) > abs(qmax)/1.e7:
                        plt.plot([v[0],w[0]] , [v[1],w[1]], color='g',
                            linestyle='--', linewidth=1)
            plt.plot(y,z, color='g', linewidth=1)
            
            haymax=False
            for i in range(1,len(tr.tics)-1):  # busca maximo local
                s0= np.sign(valoresq[i] - valoresq[i-1])
                s2= np.sign(valoresq[i] - valoresq[i+1])
                if s0+s2 !=0 :   # cubre caso de incremento=0
                    haymax=True
                    iqmax, qmax= i,valoresq[i]
                    break   # hay 0 ó 1 maximos en un tr recto
            
            alfatext=np.arctan(np.tan(tr.alfa))*180/np.pi
            if  alfatext < -80 : alfatext += 180 # angulo para etiquetas
            
            if haymax:  # interpolar parabolicamente (en s) el max & acotar
                c=np.zeros((3,3))
                c[0][0]=1.
                c[0][1]=tr.tics[iqmax-1][2]
                c[0][2]=tr.tics[iqmax-1][2]**2

                c[1][0]=1.
                c[1][1]=tr.tics[iqmax][2]
                c[1][2]=tr.tics[iqmax][2]**2

                c[2][0]=1.
                c[2][1]=tr.tics[iqmax+1][2]
                c[2][2]=tr.tics[iqmax+1][2]**2
                
                d=np.array([valoresq[iqmax-1],valoresq[iqmax],valoresq[iqmax+1]])
                a=np.linalg.solve(c,d)
                smax=-a[1]/(2*a[2]) if abs(a[2])>1.e-12 else 0.
                qmax= a[0] + a[1]*smax + a[2]*smax*smax
                
                p0= np.array([yi,zi]) + es*smax
                p1= p0 + qmax*q_scal*en
                
                y=[p0[0],p1[0]] # linea perpend, al tr en el maximo
                z=[p0[1],p1[1]]
                plt.plot(y,z, color='m', linestyle='--', linewidth=0.5)
                
                p1 += 1.5*eps*en*np.sign(qmax)  # texto del maximo
                if etiq2D:
                    texto='{:10.3e}\ns={:9.3f}'.format(abs(qmax),smax)
                    plt.text(p1[0],p1[1],texto,color='m',fontsize=10,
                    horizontalalignment='center',verticalalignment='center',
                    rotation=alfatext)
                
                d=0.1*eps*es*np.sign(qmax)  # flecha en el maximo
                plt.arrow(p0[0], p0[1], d[0],d[1], 
                head_width=0.8*eps, length_includes_head=True, 
                overhang=0.4, head_starts_at_zero=True, color='g')

            if acotar_etc: 
                # flechas al ppio+1 & final-1
                d=0.1*eps*es*np.sign(valoresq[1]) # flecha inicial
                plt.arrow(tr.tics[1][0], tr.tics[1][1], d[0],d[1], 
                head_width=0.8*eps, length_includes_head=True, 
                overhang=0.4, head_starts_at_zero=True, color='g')
                
                i=len(tr.tics)-2
                d=0.1*eps*es*np.sign(valoresq[i]) # flecha final
                plt.arrow(tr.tics[i][0], tr.tics[i][1], d[0],d[1], 
                head_width=0.8*eps, length_includes_head=True, 
                overhang=0.4, head_starts_at_zero=True, color='g')
            
            # acotacion de texto en los extremos (se acota aunq sea pequeño)
            if len(ptos[tr.inciT[0]].inciP) >1:
                p=ptos[tr.inciT[0]]
                a=p.inciP
                it_elegido=a[1] if len(a)==2 else False
                if not it_elegido or it==it_elegido:
                    p0=[yi,zi]-eps*es+ en*(valoresq[0]*q_scal+eps*np.sign(valoresq[0]))
                    if etiq2D:
                        texto='{:10.3e}'.format(abs(valoresq[0]))
                        plt.text(p0[0],p0[1],texto,color='g',fontsize=10,
                        horizontalalignment='center',verticalalignment='center',
                        rotation=alfatext)
            j=len(tr.tics)-1
            if len(ptos[tr.inciT[-1]].inciP) >1:
                p=ptos[tr.inciT[-1]]
                a=p.inciP
                it_elegido=a[1] if len(a)==2 else False
                if not it_elegido or it==it_elegido:
                    p0=[yj,zj]+eps*es+ en*(valoresq[j]*q_scal+eps*np.sign(valoresq[j]))
                    texto='{:10.3e}'.format(abs(valoresq[j]))
                    if etiq2D:
                        plt.text(p0[0],p0[1],texto,color='g',fontsize=10,
                        horizontalalignment='center',verticalalignment='center',
                        rotation=alfatext)

        elif tr.tipoT==1:
            i,j= tr.inciT[0], tr.inciT[1]
            yi,zi=ptos[i].y , ptos[i].z
            yj,zj=ptos[j].y , ptos[j].z
            incy,incz= yj-yi , zj-zi
            R,L,alfa,alfai,signo= tr.R, tr.L, tr.alfa, tr.alfai, tr.signo
            yC,zC= tr.yC, tr.zC

            y,z = [], []
            valoresq=[]
            for i in range(len(tr.tics)):  # traza la funcion
                v= np.array([tr.tics[i][0], tr.tics[i][1]])
                a=tr.qVy_plot[i]*Vy_dato + tr.qVz_plot[i]*Vz_dato
                valoresq.append(a)
                w = v+ tr.en_s(tr.tics[i][2])*a*q_scal
                y.append(w[0])
                z.append(w[1])
                if acotar_etc: # lineas punteadas
                    if i==0 and abs(valoresq[0]) > abs(qmax)/1.e7: 
                        # misma chapucilla con qmax
                        plt.plot([v[0],w[0]] , [v[1],w[1]], color='g',
                            linestyle='--', linewidth=1)
                    if i== len(tr.tics)-1 and abs(valoresq[i]) > abs(qmax)/1.e7:
                        plt.plot([v[0],w[0]] , [v[1],w[1]], color='g',
                            linestyle='--', linewidth=1)
            # si quiero suavizar yz (si hay pocos ptos), seria aqui.
            # la alternativa de aumentar n de tics es mas simple y 
            # elegante, aunque sea fuerza bruta
            plt.plot(y,z, color='g', linewidth=1)

            haymax, qmax = False, 0.
            for i in range(1,len(tr.tics)-1):  # busca maximo local
                s0= np.sign(valoresq[i] - valoresq[i-1])
                s2= np.sign(valoresq[i] - valoresq[i+1])
                if s0+s2 !=0 :   # cubre incr=0
                    haymax=True
                    if abs(valoresq[i])>abs(qmax):
                        iqmax, qmax= i,valoresq[i]
                    # si hubiese mas de un max, quedara el mayor
            
            if haymax:  # interpolar parabolicamente (en s) el max & acotar
                c=np.zeros((3,3))
                c[0][0]=1.
                c[0][1]=tr.tics[iqmax-1][2]
                c[0][2]=tr.tics[iqmax-1][2]**2

                c[1][0]=1.
                c[1][1]=tr.tics[iqmax][2]
                c[1][2]=tr.tics[iqmax][2]**2

                c[2][0]=1.
                c[2][1]=tr.tics[iqmax+1][2]
                c[2][2]=tr.tics[iqmax+1][2]**2
                
                d=np.array([valoresq[iqmax-1],valoresq[iqmax],valoresq[iqmax+1]])
                a=np.linalg.solve(c,d)
                smax=-a[1]/(2*a[2]) if abs(a[2])>1.e-12 else 0.
                qmax= a[0] + a[1]*smax + a[2]*smax*smax # calculado esta
                
                p0=tr.yz_s(smax) -Gvec  # Se ha llamado mas veces a basico_tr despues de
                                        # motor_calculo, con lo que los centros estan
                                        # otra vez desde O y la funcion tr.yz_s tiene 
                                        # un Gvec de mas que hay que corregir  :-S
                en, es = tr.en_s(smax), tr.es_s(smax)
                p1= p0 + qmax*q_scal*en


                y=[p0[0],p1[0]] # linea perpend, al tr en el maximo
                z=[p0[1],p1[1]]
                plt.plot(y,z, color='m', linestyle='--', linewidth=0.5)
                
                alfatext=np.arctan2(es[1], es[0])*180/np.pi
                if alfatext > 101 : alfatext -= 180 
                if alfatext < -80 : alfatext += 180 # angulo para etiqueta
                
                p1 += 1.5*eps*en*np.sign(qmax)  # texto del maximo
                if etiq2D:
                    texto='{:10.3e}\ns={:9.3f}'.format(abs(qmax),smax)
                    plt.text(p1[0],p1[1],texto,color='m',fontsize=10,
                    horizontalalignment='center',verticalalignment='center',
                    rotation=alfatext)
                
                d=0.1*eps*es*np.sign(qmax)  # flecha en el maximo
                plt.arrow(p0[0], p0[1], d[0],d[1], 
                head_width=0.8*eps, length_includes_head=True, 
                overhang=0.4, head_starts_at_zero=True, color='g')

            if acotar_etc: 
                # flechas al ppio+1 & final-1
                es=tr.es_s(tr.tics[1][2])
                d=0.1*eps*es*np.sign(valoresq[1]) # flecha inicial
                plt.arrow(tr.tics[1][0], tr.tics[1][1], d[0],d[1], 
                head_width=0.8*eps, length_includes_head=True, 
                overhang=0.4, head_starts_at_zero=True, color='g')
                
                i=len(tr.tics)-2
                es=tr.es_s(tr.tics[i][2])
                d=0.1*eps*es*np.sign(valoresq[i]) # flecha final
                plt.arrow(tr.tics[i][0], tr.tics[i][1], d[0],d[1], 
                head_width=0.8*eps, length_includes_head=True, 
                overhang=0.4, head_starts_at_zero=True, color='g')
            
            # Acotacion de texto en los extremos. No conviene cuando
            # el pto tiene exactamente 2 tr. De ahi el lio adicional
            # Pero se acota aunq sea pequeño (no suele molestar)

            if len(ptos[tr.inciT[0]].inciP) >1:
                p=ptos[tr.inciT[0]]
                a=p.inciP
                it_elegido=a[1] if len(a)==2 else False
                if not it_elegido or it==it_elegido:
                    es,en=tr.es_s(0.),tr.en_s(0.)
                    alfatext=np.arctan2(es[1], es[0])*180/np.pi
                    if alfatext > 101 : alfatext -= 180 
                    if alfatext < -80 : alfatext += 180
                    p0=[yi,zi]-eps*es+ en*(valoresq[0]*q_scal+
                        eps*np.sign(valoresq[0]))
                    if etiq2D:
                        texto='{:10.3e}'.format(abs(valoresq[0]))
                        plt.text(p0[0],p0[1],texto,color='g',fontsize=10,
                        horizontalalignment='center',verticalalignment='center',
                        rotation=alfatext)
            j=len(tr.tics)-1
            if len(ptos[tr.inciT[-1]].inciP) >1:
                p=ptos[tr.inciT[-1]]
                a=p.inciP
                it_elegido=a[1] if len(a)==2 else False
                if not it_elegido or it==it_elegido:
                    es,en=tr.es_s(tr.L),tr.en_s(tr.L)
                    alfatext=np.arctan2(es[1], es[0])*180/np.pi
                    if alfatext > 101 : alfatext -= 180 
                    if alfatext < -80 : alfatext += 180 # angulo para etiqueta
                    p0=[yj,zj]+eps*es+ en*(valoresq[j]*q_scal+eps*np.sign(valoresq[j]))
                    if etiq2D:
                        texto='{:10.3e}'.format(abs(valoresq[j]))
                        plt.text(p0[0],p0[1],texto,color='g',fontsize=10,
                        horizontalalignment='center',verticalalignment='center',
                        rotation=alfatext)

        else:
            print('tipo de tramo no reconocido: tr=',it,' tipo=',tipoT)
            quit()
    
    
    
    # ventana de desplazamientos por cortante uV (en 3D)
    
    vent_uV=plt.figure('Alabeo por Cortante')
    ax=vent_uV.gca(projection='3d')
    #ax.set_aspect('equal') sustituido por llamada a yz3D_iguales()
    
    if uV_despl=='': uV_despl=1  # se inicializa 'afuera de la seccion'

        # dibuja vector V (la flecha a mano porque arrow no va en 3D)
    alfa=np.arctan2(Vz_dato,Vy_dato)
    a= Gvec[0]+ 10*eps*np.cos(alfa)
    b= Gvec[1]+ 10*eps*np.sin(alfa)
    ax.plot( [Gvec[0],a], [Gvec[1],b], [0.,0.],  lw=2, color ='r' )

    beta = alfa+ 160.*np.pi/180. 	# semiangulo flecha 20º
    aa= a+ 1.7*eps*np.cos(beta)
    bb= b+ 1.7*eps*np.sin(beta)
    ax.plot( [a,aa], [b,bb], [0.,0.],  lw=3, color ='r' )

    beta = alfa- 160.*np.pi/180. 	# semiangulo flecha 20º
    aa= a+ 1.7*eps*np.cos(beta)
    bb= b+ 1.7*eps*np.sin(beta)
    ax.plot( [a,aa], [b,bb], [0.,0.],  lw=3, color ='r' )
    
    if etiq3D:
        alfa += 0.27
        ax.text(Gvec[0]+ 10*eps*np.cos(alfa), Gvec[1]+ 10*eps*np.sin(alfa),
                0., 'V', color='r', fontsize=14, fontweight='normal',
                horizontalalignment='center',verticalalignment='center'  )
        
        # dibuja ejes principales
    rel=(Ieta/Ipsi)**0.25

    c=12*eps 	# el eje psi mas largo
    a, b= Gvec+ c*np.array(psi.es) ,  Gvec- c*np.array(psi.es)
    ax.plot([a[0],b[0]],[a[1],b[1]], color='b', linewidth=1)

    c=12*eps/rel	# eje eta mas corto
    a, b = Gvec + c*np.array(eta.es), Gvec- c*np.array(eta.es)
    ax.plot([a[0],b[0]],[a[1],b[1]], color='b', linewidth=1)

        # busca umax (no se puede combinar como u medio):
    umax=0.
    for tr in tramos.values():
        for i in range(len(tr.tics)):
            u=Vy_dato*tr.uVy_plot[i]+Vz_dato*tr.uVz_plot[i]
            if abs(umax)<abs(u):umax=u

        # asignacion del corrimiento para trazar ux:
    if uV_despl == 2:  # centrado al promedio
        uV_medio= Vy_dato*uVy_medio + Vz_dato*uVz_medio
        d=-uV_medio
        print('uV_medio es=', uV_medio)
    elif uV_despl == 1:  # afuera de la seccion
        d=-umax
    else:   # sera uV_despl=0, tal como se calcula
        d=0.
    
        # trazado basico de seccion & los uV
    for it, tr in tramos.items():
        ax.plot(tr.tics[:,0],tr.tics[:,1], color='k')
        uV_plot= Vy_dato*np.array(tr.uVy_plot)+Vz_dato*np.array(tr.uVz_plot)
        ax.plot3D(tr.tics[:,0],tr.tics[:,1],uV_plot[:]+d,
            linewidth=2, color='grey')
    
    if peque3D:
        # un puntin perpendicular, que quede ux mas pegado a la seccion
        a=[4*umax, 5*umax, 4*umax]
        ax.plot3D([Gvec[0]], [Gvec[1]], [a[uT_despl]] )
        a=[3*umax, 2*umax, 4*umax]
        ax.plot3D([Gvec[0]], [Gvec[1]], [-a[uT_despl]] )
      
        # en ptos: lineas guia y etiquetar ptos y despl
    for ip,p in ptos.items():
        a= Vy_dato*p.uVy+Vz_dato*p.uVz
        ax.plot([p.y, p.y],[p.z, p.z],[0.,a+d],'--',lw=0.8,color ='grey')
        if etiq3D:
            ax.text(p.y, p.z, a+d , '{:.4}'.format(a), 'x', 
                fontsize=11, color='grey', fontweight='normal')
            ax.text(p.y, p.z, 0, ip, 'x', fontsize=10, color='g', 
                fontweight='bold') 
    
        # trazado del plano de ajuste
    uf= lambda y,z: (Vy_dato*(C_tensor[0][0]*y+ C_tensor[0][1]*z+ C0y)+
                     Vz_dato*(C_tensor[1][0]*y+ C_tensor[1][1]*z+ C0z)+d)
    for tr in tramos.values():
        a=[]
        for i in range(len(tr.tics)):
            a.append(uf(tr.tics[i][0], tr.tics[i][1]))
        ax.plot3D(tr.tics[:,0],tr.tics[:,1], a, linewidth=1,color='m')
        
        # dibuja la linea de ux=0  (ay+bz+c=0)
    a=Vy_dato*C_tensor[0][0]+Vz_dato*C_tensor[1][0]
    b=Vy_dato*C_tensor[0][1]+Vz_dato*C_tensor[1][1]
    c=Vy_dato*C0y+Vz_dato*C0z + d
    if a==0:
        p, e = np.array([0, -c/b]), np.array([1., 0.])
    else:
        k=np.array([-b,a])/np.sqrt(a*a+b*b)
        p, e = np.array([-c/a, 0]), k

    LNuV=Recta(p, e, np.array([-e[1],e[0]]))
    e=np.array([Vy_dato,Vz_dato])/np.sqrt(Vy_dato**2+Vz_dato**2)
    LV=Recta(Gvec, e, np.array([-e[1],e[0]]))
    
    s1,s2=LV.intersec(LNuV)
    p= LV.p + s1*LV.es
    p1=p+ 20*eps*LNuV.es
    p2=p- 20*eps*LNuV.es
    ax.plot(  [ p1[0],p2[0] ],  [ p1[1],p2[1] ], [ 0.,0. ], '--',
            color='r',  linewidth=0.5 )
    if etiq3D:
        ax.text(p1[0], p1[1], 0., '$u_{x}=0$',  
                fontsize=11, color='r', fontweight='normal')

    if para: 
        #yz3D_iguales(ax)
        plt.show()




    


yzsonPR,esferico = '',''
anchoy, anchoz, eps_L = 0., 0., 0.
theta_secc = 0.
uV_despl, uT_despl, etiq2D, etiq3D, peque3D = '', '', 1, 1, 0

presenta_elige()

        
# ---------------------------------------------------
# ----    Ventana GUI de interaccion principal   ----
#----------------------------------------------------


v9 = Tk()
v9.geometry('740x600+6+6')
v9.title('ThinSecBeam 1.0 ')

estilo = ttk.Style()
estilo.configure('jc.TButton', foreground='green')
estilo.configure('jc.TLabelframe.Label', foreground ='green')
estilo.configure('jc_blue.TButton',backgroung='#e6e6ff')
estilo.configure('jc_red.TButton',backgroung='#ffe6e6', foreground='ffffe6')
#estilo.configure('jc.TFrame', background='#e6e6ff')
#fondoazul = ttk.Frame(v9, padding='4', style='jc.TFrame')

hard_nfilas=24
filas_ptos, filas_tramos= [], []
nfcompleto=''
dir_home=path.expanduser('~')

fondoazul = ttk.Frame(v9, padding='4')
fondoazul.grid(column=0, row=0, rowspan=4, columnspan=3, sticky=(N, W, E, S))

frame_ptos= ttk.Labelframe(fondoazul, text='Puntos',style='jc.TLabelframe',
    width=200,height=500)
frame_tramos= ttk.Labelframe(fondoazul, text='Tramos',style='jc.TLabelframe', 
    width=450,height=500)
frame_ptos.grid(column=0,row=0,rowspan=6,ipadx=3)
frame_tramos.grid(column=1,row=0,rowspan=6,ipadx=3)

frame_M= ttk.Labelframe(fondoazul, text='Flector',style='jc.TLabelframe', 
    width=170,height=170)
frame_M.grid(column=2,row=0,padx=3, pady=3, sticky='s')

frame_V= ttk.Labelframe(fondoazul, text='Cortante',style='jc.TLabelframe', 
    width=170,height=130)
frame_V.grid(column=2,row=1,padx=3, pady=3)

frame_T= ttk.Labelframe(fondoazul, text='Torsor',style='jc.TLabelframe', 
    width=170,height=80)
frame_T.grid(column=2,row=2,padx=3, pady=5)

texto=' Cerrar\ngráficas'
boton_cierraplots=ttk.Button(fondoazul, text=texto, width=9,
    command=cierraplots).grid(column=2, row=3, padx=3, pady=3)

frame_botones=ttk.Frame(fondoazul, width=170, height=110)
frame_botones.grid(column=2,row=5, padx=3, pady=3)

boton_dibuPt=ttk.Button(frame_ptos, text='Dibuja',width=6, 
    command=lambda: pintaPtosTramos(1,1,1,0))
boton_dibuPt.grid(column=0, row=hard_nfilas+1, columnspan=2, sticky='sw',
    padx=3, pady=3)

boton_dibuTr=ttk.Button(frame_tramos, text='Dibuja',width=6,
    command= lambda: pintaPtosTramos(1,1,1,1))
boton_dibuTr.grid(column=0, row=hard_nfilas+1,columnspan=3, sticky='sw', 
    padx=3, pady=3)
    
boton_opc=ttk.Button(frame_tramos, text='Opciones',width=8,
    command=entraOpc).grid(column=3, row=hard_nfilas+1, sticky='sw', padx=3, pady=3)

boton_yaTr=ttk.Button(frame_tramos,text='Calcula',width=8, command=haz_ficha)
boton_yaTr.grid(column=4, row=hard_nfilas+1, columnspan=2, sticky='se', padx=3, pady=3)

# en mas,cargar &guardar el grid aparte, para poder destruirlos etc:

boton_mas =ttk.Button(fondoazul,text='...+',width=5,command=a_cargar_mas)
boton_mas.grid(column=2,row=4, padx=5, pady=(45,3))

boton_txt=ttk.Button(frame_botones,text='Informe',width=8,
    command=haz_informe, state='disabled')
boton_txt.grid(column=0, row=0, sticky='sw', columnspan=2, padx=3, pady=2)

boton_cargar=ttk.Button(frame_botones, text='Cargar', width=6,command=a_cargar)
boton_cargar.grid(column=1, row=0, sticky='se', padx=3, pady=2)

boton_guardar=ttk.Button(frame_botones, text='Guardar',width=8,
    command=a_guardar, state='disabled')
boton_guardar.grid(column=0, row=1,sticky='sw', padx=3, pady=2)

boton_salir=ttk.Button(frame_botones, text='Salir',width=6, style='jc.TButton',
    command=a_salir).grid(column=1, row=1,sticky='se', padx=3, pady=3)

label_ptos1=ttk.Label(frame_ptos,text='pto').grid(column=0, row=0, padx=3, pady=2)
label_ptos2=ttk.Label(frame_ptos,text='y').grid(column=1, row=0, padx=3, pady=2)
label_ptos3=ttk.Label(frame_ptos,text='z').grid(column=2, row=0, padx=3, pady=2)

label_tr1=ttk.Label(frame_tramos,text='tr').grid(column=0, row=0, padx=3, pady=2)
label_tr2=ttk.Label(frame_tramos,text='ini').grid(column=1, row=0, padx=3, pady=2)
label_tr3=ttk.Label(frame_tramos,text='fin').grid(column=2, row=0, padx=3, pady=2)
label_tr4=ttk.Label(frame_tramos,text='e').grid(column=3, row=0, padx=3, pady=2)
label_tr5=ttk.Label(frame_tramos,text='R').grid(column=4, row=0, padx=3, pady=2)
label_tr6=ttk.Label(frame_tramos,text=u'\u00B1',foreground='green').grid(
    column=5, row=0, padx=3, pady=2)


# creamos los Entry para los ptos

for ifilas in range(hard_nfilas):
    fila=[]
    for j in range(3):
        ancho=3 if j==0 else 10
        e=ttk.Entry(frame_ptos, width=ancho)
        e.grid(row=ifilas+1, column=j)
        fila.append(e)
    filas_ptos.append(fila)


fondoazul.columnconfigure(0, weight=2)
fondoazul.columnconfigure(1, weight=3)
fondoazul.columnconfigure(2, weight=1)
fondoazul.rowconfigure(0, weight=1)

# y los Entry para los tramos
for ifilas in range(hard_nfilas):
    fila=[]
    a=[3,3,3,10,10,2]
    for j in range(6):
        ancho=a[j]
        e=ttk.Entry(frame_tramos, width=ancho)
        e.grid(row=ifilas+1, column=j)
        fila.append(e)
    filas_tramos.append(fila)


# los entry para las cargas Vy Vz T My Mz

label_Vy=ttk.Label(frame_V,text='Vy').grid(
    column=0, row=0, padx=3, pady=2)
label_Vz=ttk.Label(frame_V,text='Vz').grid(
    column=0, row=1, padx=3, pady=2)
entry_Vy=ttk.Entry(frame_V, width=10)
entry_Vy.grid(column=1, row=0, padx=3, pady=2)
entry_Vz=ttk.Entry(frame_V, width=10)
entry_Vz.grid(column=1, row=1, padx=3, pady=2)

boton_goV=ttk.Button(frame_V, text='Go!',width=4,# style='jc.TButton', 
    command=haz_V, state='disabled')
boton_goV.grid(column=0, row=2, columnspan=2, padx=3, pady=5)

label_T=ttk.Label(frame_T,text='Mx').grid(column=0, row=0, padx=3, pady=2)
entry_T=ttk.Entry(frame_T, width=10)
entry_T.grid(column=1, row=0, padx=3, pady=2)
boton_goT=ttk.Button(frame_T, text='Go!',width=4,command=haz_T, state='disabled')
boton_goT.grid(column=0, row=1, columnspan=2, padx=3, pady=5)

label_Nx=ttk.Label(frame_M,text='Nx').grid(column=0, row=0, padx=3, pady=2)
label_My=ttk.Label(frame_M,text='My').grid(column=0, row=1, padx=3, pady=2)
label_Mz=ttk.Label(frame_M,text='Mz').grid(column=0, row=2, padx=3, pady=2)
entry_Nx=ttk.Entry(frame_M, width=10)
entry_Nx.grid(column=1, row=0, padx=3, pady=2)
entry_My=ttk.Entry(frame_M, width=10)
entry_My.grid(column=1, row=1, padx=3, pady=2)
entry_Mz=ttk.Entry(frame_M, width=10)
entry_Mz.grid(column=1, row=2, padx=3, pady=2)
boton_goM=ttk.Button(frame_M, text='Go!',width=4,command=haz_M, state='disabled')
boton_goM.grid(column=0, row=3, columnspan=2, padx=3, pady=5)

v9.columnconfigure(0, weight=1)
v9.rowconfigure(0, weight=1)

if elige == 'default':
    a_cargar()
    elige='' # por si quiero cargar otro prb tras el default
    haz_ficha(recarga=False, para=False)
    haz_M(para=False)
    haz_V(para=False)
    haz_T(para=True)

boton_cargar.focus()

v9.protocol('WM_DELETE_WINDOW', a_salir)
v9.mainloop()










