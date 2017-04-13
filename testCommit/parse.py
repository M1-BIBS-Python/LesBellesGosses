#!/usr/bin/env python
#-*- coding : utf8 -*-

"""
Author : LesBellesGosses
Date : 12/04/2017
Description : Projet Barstar
"""

from math import sqrt
import numpy

def ParsingPDB (pdbFile):#fonction qui parse un fichier pdb
   
    infile = open(pdbFile, "r")
    lines = infile.readlines()

    dico_models={} # une cle associee a chaque conformation de la molecule

    for line in lines:
        if line[:5:] == "MODEL":
            dico_molecule=line[13:].strip()
            dico_models[dico_molecule] = {} 
            dico_models[dico_molecule]["chains"] = []
            
        if line[:4:] == 'ATOM':                     
        
            chain = line[21]
            if chain not in dico_models[dico_molecule].keys():
                dico_models[dico_molecule][chain] = {}
                dico_models[dico_molecule]["chains"].append(chain)
                dico_models[dico_molecule][chain]["reslist"]=[]
        
            res = line[23:26]
            if res not in dico_models[dico_molecule][chain].keys() :
                dico_models[dico_molecule][chain]["reslist"].append(res)
                dico_models[dico_molecule][chain][res] = {}
                dico_models[dico_molecule][chain][res]["atomlist"]=[]
            
            atom = line[13:16]
            dico_models[dico_molecule][chain][res]["atomlist"].append(atom)
            dico_models[dico_molecule][chain][res][atom] = {}
    
            dico_models[dico_molecule][chain][res][atom]['x'] = line[31:38]
            dico_models[dico_molecule][chain][res][atom]['y'] = line[39:46]
            dico_models[dico_molecule][chain][res][atom]['z'] = line[47:54]
            dico_models[dico_molecule][chain][res][atom]['id'] = line[7:11]
        
            dico_models[dico_molecule][chain][res]['resname'] = line[17:20]
        
    infile.close()
    return(dico_models)


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
            x=[float(coor1[0]),float(coor2[0])]
            y=[float(coor1[1]),float(coor2[1])]
            z=[float(coor1[2]),float(coor2[2])]
            
            centre=CM(x,y,z)
            #3. Calculer la distance entre les residus et le centre de masse
            distance=Distance(float(coor1[0]),float(coor1[1]),float(coor1[2]),centre[0],centre[1],centre[2])
            distance2=Distance(float(coor2[0]),float(coor2[1]),float(coor2[2]),centre[0],centre[1],centre[2])
        print("d1=")
        print(distance)
        print("d2=")
        print(distance2)
    return(distance)

# je propose qu'on fasse une fonction de calcul du centre de masse comme ci-dessous, je trouve que c'est plus simple et plus general
def CM(listx,listy,listz):
    x=sum(listx)/float(len(listx))
    y=sum(listy)/float(len(listy))
    z=sum(listz)/float(len(listy))
    coords=[x,y,z]
    return coords
    
def Distance(x1,y1,z1,x2,y2,z2): #Calcule la distance entre deux points
    return(sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2))

def RMSD(list_delta): #calcul de RMSD
    distcarre=[]
    for delta in list_delta:
        distcarre.append(delta**2)
    return (sqrt((sum(distcarre))/float(len(list_delta))))

def giration(listx,listy,listz,dico1):
    #Dapres ce que jai compris : cest la distance moyenne du centre dun nuage de point avec lensemble des points du nuages
    
    N =0#Nombre de residus totaux
    somme=0
    #1. Calcul de la distance entre le residus et le centre de masse
    for atom in dico1['atomlist']:
        coor=[float(dico1[atom]["x"]),float(dico1[atom]["y"]),float(dico1[atom]["z"])]
        centre=CM(listx,listy,listz)
        distance=Distance(coor[0],coor[1],coor[2],centre[0],centre[1],centre[2])
        somme=somme+distance #On fait la somme de toutes les distances
        N=N+1
    #2. On retoune la racine de la moyenne de toutes ces distances au carres
    return (sqrt((1/N)*sum**2))


if __name__ == '__main__':
    
    import argparse,os,glob,shutil,sys

    ########################################################
    #                   Arguments
    ########################################################

    parser = argparse.ArgumentParser()
    parser.add_argument("-p", help="Path to directory where pdb files are stored or Path to PDB file") 
    parser.add_argument("-a", help="Type of analysis: global or local or both")
    args = parser.parse_args()
    if args.p == None or args.a == None: # si l'un des arguments est vide
        parser.print_help()
        parser.exit()
    path=args.p
    analyse=args.a

    ########################################################
    #               Etapes preliminaires
    ########################################################

    #on supprime le dossier contenant les resultats de l'execution precedente (s'il existe) afin d'eviter les chevauchements
    try:
        shutil.rmtree("%s/PythonProgResults"%os.path.dirname(path)) 
    except OSError:
        pass
    
    try:
        os.chdir(path)
        fichiers=glob.iglob("*.pdb") # si le path est vers un dossier, on lira tous les fichiers pdb
    except OSError:
        fichiers=glob.iglob("%s"%path) # si le path est vers un fichier, on lira ce ficher
        pass

    #########################################################
    #                       Main
    #########################################################

    if (analyse == "global"):
        for fichier in fichiers:
            dico=ParsingPDB(fichier)
            
            #calcul RMSD de chaque conformation par rapport a la structure d'origine
            dico_RMSD={}
            for key in dico:
                list_delta=[]
                for chain in dico[key]["chains"]:
                    for res in dico[key][chain]["reslist"]:
                        list_delta.append(Distance(float(dico[key][chain][res]['CA ']['x']),float(dico[key][chain][res]['CA ']['y']),float(dico[key][chain][res]['CA ']['z']),float(dico['0'][chain][res]['CA ']['x']),float(dico['0'][chain][res]['CA ']['y']),float(dico['0'][chain][res]['CA ']['z']))) # compare dist entre les Ca
                coords=RMSD(list_delta)
                dico_RMSD[key]=coords
            
    elif (analyse == "local"):
        for fichier in fichiers:
            do sth
    else :
        for fichier in fichiers:
            do sth
    #########################################################
    #            Ecrire les fichiers de sortie
    #########################################################
    
    os.mkdir("%s/PythonProgResults"%os.path.dirname(path)) #ce dossier sert a stocker les fichiers sortis
    out=open("%s/PythonProgResults/output_analyse_global"%os.path.dirname(path),"w")
    with open(fichier, "r") as filin:
        line=filin.readline()
        while line != "ENDMDL\n":
            line=filin.next()
            out.write(line)
    filin.close()
    out.write("\nRMSD results\n")
    
    for i in range(len(dico)):
        print i
        out.write("\nconformation %s: %s"%(i,dico_RMSD["%s"%i]))
        
    out.close()


     
