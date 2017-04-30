#!/usr/bin/env python
#-*- coding : utf8 -*-

"""
Author : LesBellesGosses
Description : Projet Barstar
"""
from math import sqrt
import numpy as np
import matplotlib.pyplot as plt

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


#Fonction qui lit les premieres lignes du fichier pdb et mets dans une liste les valeurs du temps
def Temps (pdbFile):
    """but : creer un liste contenant les valeurs du temps
    input : le nom d'un fichier pdb
    output : une liste
    """
    temps = []
    infile = open (pdbFile, "r")
    lines = infile.readlines()
    for line in lines :
        if line[:5:] == "TITLE":
            timet=line[65:80].strip()
            temps.append(timet)
    infile.close()
    return(temps)
    

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
    return glob # dico contient le CM de chaque residu et le CM de la prot
 
    
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

    
def graph(ordonnee,abscisse,ordonne2,titre,titrey,titrex,type_graph):
    """but : Representer les resultats sous forme de graphique pour les interpreter
    input :les coordonnees x,y et y2 : y2 permet de superposer 2 graphs si on le souhaite(si on veut faire un graph unique, y2 sera vide) et le type de graph (point ou line)
    titrey et titrex sont les legendes des coordonees des axes y et x
    output : un graphique
    """
    absc=[] #Liste qui va contenir les coordonnees de labscisse
    ordo=[] #Liste qui va contenir les coordonnnees de lordonees
    ordo2=[] #Pour le cas ou on veut superposer des graphs

    if type(ordonnee) is list : #si cest une liste
        ordo=ordonnee

    else: #sinon cest un dico
        #j'ai trouve pk cette fct bloque pour local: les dicos ne contiennent pas la cle "0" (diff que global"), donc je propose une solution ci-dessous et au moins ca marche pour un test, apres il faut verifier si c'est correcte
        if "0" in ordonnee:
            for i in range(len(ordonnee)): #on parcours le dictionnaire
                ordo.append(ordonnee["%s"%i]) #Et on range dans une liste les differentes valeurs contenues dans le dico, dans lordre dans lesquelles on les a trouve
        else:
            for i in range(1,len(ordonnee)):
                ordo.append(ordonnee["%s"%i])

    if type(ordonne2) is list: #Si cest une liste
        ordo2=ordonne2
    else:
        if "0" in ordonnee:
            for i in range(len(ordonne2)):
                ordo2.append(ordonne2["%s"%i])
        else:
            for i in range(1,len(ordonne2)):
                ordo2.append(ordonne2["%s"%i])


    if type(abscisse) is list:
        absc=abscisse
    else:
        for i in range(len(abscisse)):
            absc.append((abscisse["%s"%i]))


    #Si la liste ordo2 est vide, on fait une simple representation
    if not ordo2:
        y=np.array(ordo)
        x=np.array(absc)
        plt.title(titre)
        plt.xlabel(titrex)
        plt.ylabel(titrey)
        if type_graph=="point": #Si on veut representer des points
            plt.scatter(x,y)

        if type_graph=="line": #Si on veut representer des lignes
            plt.plot(x,y)


    #Sinon on fait des graphs superposes
    else:
        y=np.array(ordo) #RMSD
        y2=np.array(ordo2)
        x=np.array(absc) #La conformation=numero du modele
        plt.title(titre)
        plt.xlabel(titrex)
        plt.ylabel(titrey)
        

        if type_graph =="point":
            plt.scatter(x,y,c='red')
            plt.scatter(x,y2,c='blue')

        if type_graph=="line":
            plt.plot(x,y,c='red')
            plt.plot(x,y2,c='blue')


     #On affiche le graph
   
    plt.show()
