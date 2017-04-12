#!/usr/bin/env python
#-*- coding : utf8 -*-

"""
Contact: 
Date: 06/03/2017
Description: 
- A function parsing a PDB file into a dictionary
"""
from math import sqrt
import numpy

def ParsingPDB (pdbFile):#fonction qui parse un fichier pdb
   
    infile = open(pdbFile, "r")
    lines = infile.readlines() 

    dico_molecule = {}									
    chainList = []
    rList = []

    print("ok")
    for line in lines:
        if line[:4:] == 'ATOM':						
        
            chaine = line[21]
            if chaine not in chainList:
                chainList.append(chaine)
                dico_molecule[chaine] = {}
                resList=[]
        
            res = line[23:26]
            if res not in resList :
                resList.append(res)
                dico_molecule[chaine][res] = {}
                atomList=[]
                rList=[]
            
            atom = line[13:16]
            if atom not in atomList:
                atomList.append(atom)
                dico_molecule[chaine][res][atom] = {}
            
            if line[17:20] not in rList:
                rList.append(line[17:20])
    
            dico_molecule[chaine][res][atom]['x'] = line[31:38]
            dico_molecule[chaine][res][atom]['y'] = line[39:46]
            dico_molecule[chaine][res][atom]['z'] = line[47:54]
            dico_molecule[chaine][res][atom]['id'] = line[7:11]
        
            dico_molecule[chaine][res]['resname'] = rList
            dico_molecule[chaine]['reslist'] = resList
            dico_molecule[chaine][res]['atomlist'] = atomList
        
    infile.close()
    return(dico_molecule)


#Calcule de la distance entre 2 residus selon 2 methodes de calcule:
    #1. Distance la plus courte entre les 2 residus
    #2. Distance entre les centres de masse

def Distance_res_min(dico1,dico2):   #Mode de calcul 1 : distance la plus courte
#Dico 1 et Dico2 sont des dico de residus
    min=200 #On fixe une distance min
    
    for atom1 in dico1['atomlist']:#Calcul de la distance entre les 2 residus
        coor1=[dico1[atom1]["x"],dico1[atom1]["y"],dico1[atom1]["z"]]
        for atom2 in dico2["atomlist"] :
            coor2=[dico2[atom2]["x"],dico2[atom2]["y"],dico2[atom2]["z"]]
            distance=Distance(float(coor1[0]),float(coor1[1]),float(coor1[2]),float(coor2[0]),float(coor2[1]),float(coor2[2]))
            
            #print(distance)
            if distance>0:   #A mettre ? car si=0, les residus sont a la meme position...
                if distance<min: #On en deduit la distance min
                    min=distance
    print(min)
    return(min)


def Distance_res_masse(dico1,dico2):
     #Mode de calcul 2 : entre les centres de masse
     #x,y, et z sont les coordoonnees du centre
     #Dico1 et 2 sont des dicos de residus
            
    
    #1.definir les coordonnes du centre de masse X,Y,Z
    x=0
    y=0
    z=0
    
    X=0
    Y=0
    Y=0
    cpt=0
    distance=0
    distance2=0
    #Les coordonnees se definissent comme une moyenne des coordonnees des 2 residus
    for atom1 in dico1['atomlist']:
        coor1=[dico1[atom1]["x"],dico1[atom1]["y"],dico1[atom1]["z"]]
        for atom2 in dico2["atomlist"] :
            coor2=[dico2[atom2]["x"],dico2[atom2]["y"],dico2[atom2]["z"]]
            #On fait la somme des coordonnees de tous les atomes
            x=float(coor1[0])+float(coor2[0])
            y=float(coor1[1])+float(coor2[1])
            z=float(coor1[2])+float(coor2[2])
            cpt=cpt+1 
           
            X=x/cpt #On divise par le nombre d'atomes
            Y=y/cpt
            Z=z/cpt
            
            #3. Calculer la distance entre les residus et le centre de masse
            distance=Distance(float(coor1[0]),float(coor1[1]),float(coor1[2]),X,Y,Z)
            distance2=Distance(float(coor2[0]),float(coor2[1]),float(coor2[2]),X,Y,Z)
        print("d1=")
        print(distance)
            
        print("d2=")
        print(distance2)
    return(distance)
    
def Distance(x1,y1,z1,x2,y2,z2): #Calcule la distance entre deux points
    return(sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2))


if __name__ == '__main__':
    res1=ParsingPDB("/Users/mathildebertrand/Desktop/M1-S2/python/tp3/Documents/arginine.pdb")
    res2=ParsingPDB("/Users/mathildebertrand/Desktop/M1-S2/python/tp3/Documents/arginine.pdb")
    #Distance_res_min(res1['A'][' -3'],res2['A'][' -3'])
    Distance_res_masse(res1['A'][' -3'],res2['A'][' -3'])
