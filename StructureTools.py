#!/usr/bin/env python
#-*- coding : utf8 -*-

"""
Author : LesBellesGosses
Description : Projet Barstar
"""
from math import sqrt

#parser un fichier pdb
def ParsingPDB (pdbFile):
    """but : creer un dictionnaire a partir d'un fichier pdb
    input : le nom d'un fichier pdb
    output : un dictionnaire
    """
    infile = open(pdbFile, "r")
    lines = infile.readlines()

    dico_models={} # une cle associee a chaque conformation de la molecule

    for line in lines:

        if line[:5:] == "MODEL": #Si la ligne commence par MODEL,On rajoute le numero de conformation

            dico_molecule=line[9:14].strip()
            dico_models[dico_molecule] = {}
            dico_models[dico_molecule]["chains"] = []

        if line[:4:] == 'ATOM':                     #Si la ligne commence par ATOM

            chain = line[21]
            if chain not in dico_models[dico_molecule].keys():
                dico_models[dico_molecule][chain] = {} # je me demande si on garde la cle chain car il n'y pas de chain dans fichier pdb
                dico_models[dico_molecule]["chains"].append(chain)
                dico_models[dico_molecule][chain]["reslist"]=[]

            res = line[23:26].strip()
            if res not in dico_models[dico_molecule][chain].keys() :
                dico_models[dico_molecule][chain]["reslist"].append(res)
                dico_models[dico_molecule][chain][res] = {}
                dico_models[dico_molecule][chain][res]["atomlist"]=[]

            atom = line[13:16].strip()
            dico_models[dico_molecule][chain][res]["atomlist"].append(atom)
            dico_models[dico_molecule][chain][res][atom] = {}

            dico_models[dico_molecule][chain][res][atom]['x'] = line[31:38]
            dico_models[dico_molecule][chain][res][atom]['y'] = line[39:46]
            dico_models[dico_molecule][chain][res][atom]['z'] = line[47:54]
            dico_models[dico_molecule][chain][res][atom]['id'] = line[7:11]

            dico_models[dico_molecule][chain][res]['resname'] = line[17:20]

    infile.close()
    return(dico_models)

#calcul le centre de masse (en negligeant la masse atomique)
def CM(listx,listy,listz):
    """but : calculer le centre de masse
    input : l'ensemble des abscisses, des ordonnees et des cotes en liste
    output : une liste contenant l'abscisse, l'ordonnee et la cotes du centre de masse
    """
    x=sum(listx)/float(len(listx))
    y=sum(listy)/float(len(listy))
    z=sum(listz)/float(len(listy))
    coords=[x,y,z]
    return coords


#Calcul de distance entre deux points
def Distance(x1,y1,z1,x2,y2,z2):
    """but : calculer la distance dans l'espace tridimensionnelle
    input : les coordonnees de deux points
    output : la distance entre ces deux points
    """
    return(sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2))


#calcul de RMSD
def RMSD(list_delta):
    """but : calculer le RMSD
    input : une liste de distances
    output : la valeur de RMSD
    """
    distcarre=[]
    for delta in list_delta:
        distcarre.append(delta**2)
    return (sqrt((sum(distcarre))/float(len(list_delta))))


#creer un dictionnaire de centre de masse pour une proteine
def CMglob(dico):
    """but : calculer le centre de masse de chaque residu ainsi que le centre de masse de la proteine
    input : dico de proteine pour une conformation donnee (dico[conformation])
    output : un dictionnaire contenant le centre de masse de chaque residu et le centre de masse de la proteine
    """
    globx=[] #list permet de stocker le x de tous les atomes d'une prot
    globy=[]
    globz=[]
    glob={}
    glob["residulist"]=[]
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
            glob[res]=CM(listx,listy,listz)
            glob["residulist"].append(dico[chain][res]["resname"])
    glob["prot"]=CM(globx,globy,globz)
    
    
#creer des classes
def createClass(dico, bestscore, nbcl) :
    """but : un dictionnaire permettant de classer les valeurs d'un dictionnaire en nombre de classes que les utilisateurs souhaitaient
    input : dictionnaire, maximum de la liste, nombre de classes
    output : dictionnaire contenant chaque element de la liste comme cle, et sa classe comme valeur
    """
    classe={} #dictionnaire de la classe aux elements
    compteur=0
    seuil=0
    while len(classe) != nbcl:
        compteur=compteur+1
        classe[compteur]=[]
        seuil=(nbcl-compteur)*(bestscore/nbcl)


        for cle,element in dico.items(): #on parcours le dico
            if element >= seuil:
                classe[compteur].append(cle) #on range le numero du residu
                           
    dico_etoclass={} # dictionnaire d'element a la classe
    for key in classe:
        for elem in classe[key]:
            if not elem in dico_etoclass:
                dico_etoclass[elem]=key

    return dico_etoclass
