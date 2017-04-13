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

def giration(dico):#c'est le dico[key] qu'on passe ici
    listCMres=[] #list permet de stocker le centre de masse de chaque residu
    globx=[] #list permet de stocker le x de tous les atomes d'une prot
    globy=[]
    globz=[]
    for chain in dico["chains"]:
        for res in dico[chain]["reslist"]:
            listx=[]
            listy=[]
            listz=[]
            for atom in dico[chain][res]["atomlist"]:
                listx.append(float(dico[chain][res][atom]['x']))
                listy.append(float(dico[chain][res][atom]['y']))
                listz.append(float(dico[chain][res][atom]['z'])) 
            globx.extend(listx)
            globy.extend(listy)
            globz.extend(listz)
            listCMres.append(CM(listx,listy,listz))
    CMprot=CM(globx,globy,globz)
    list_dist=[]
    for xyz in listCMres:
        list_dist.append(Distance(CMprot[0],CMprot[1],CMprot[2],xyz[0],xyz[1],xyz[2]))
    return max(list_dist)   

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

    #########################################################
    #                       Main
    #########################################################

    if (analyse == "global"):
        for fichier in fichiers:
            dico=ParsingPDB(fichier)
            
            #calcul RMSD de chaque conformation par rapport a la structure d'origine
            dico_RMSD={}
            dico_Giration={}
            for key in dico:
                list_delta=[]
                dico_Giration[key]=giration(dico[key])
                for chain in dico[key]["chains"]:
                    for res in dico[key][chain]["reslist"]:
                        list_delta.append(Distance(float(dico[key][chain][res]['CA ']['x']),float(dico[key][chain][res]['CA ']['y']),float(dico[key][chain][res]['CA ']['z']),float(dico['0'][chain][res]['CA ']['x']),float(dico['0'][chain][res]['CA ']['y']),float(dico['0'][chain][res]['CA ']['z']))) # compare dist entre les Ca         
                dico_RMSD[key]=RMSD(list_delta)
            
    elif (analyse == "local"):
        for fichier in fichiers:
            #do sth
    else :
        for fichier in fichiers:
            #do sth
    #########################################################
    #            Ecrire les fichiers de sortie
    #########################################################
    
    os.mkdir("%s/PythonProgResults"%os.path.dirname(path)) #ce dossier sert a stocker les fichiers sortis
    out=open("%s/PythonProgResults/output_analyse_global"%os.path.dirname(path),"w")
    
    #ecrire la structure d'origine
    with open(fichier, "r") as filin:
        line=filin.readline()
        while line != "ENDMDL\n":
            line=filin.next()
            out.write(line)
    filin.close()
    
    #ecrire dans l'ordre les resultats du RMSD de la comparaison des conformations avec la structure d'origine
    out.write("\nConformation \t Giration \t\t RMSD results\n")
    for i in range(len(dico)):
        out.write("\n%s: \t\t %.12f \t\t %.12f"%(i,dico_RMSD["%s"%i],dico_Giration["%s"%i]))
        
    out.close()


     
