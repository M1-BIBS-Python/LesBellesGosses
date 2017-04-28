#!/usr/bin/env python
#-*- coding : utf8 -*-
#~ python parse.py -p /Users/mathildebertrand/Desktop/M1-S2/python/Git/LesBellesGosses -a local
"""
Author : LesBellesGosses
Date : 12/04/2017
Description : Projet Barstar
"""

from math import sqrt
import matplotlib.pyplot as plt
import numpy as np
import argparse,os,glob,shutil,sys

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
    

#Enfouissement de chaque residu dans une proteine (pour l'analyse locale)
#a modifier pour tenir compte le fichier de ref
def Enfouissement(dico):#c'est le dico qu'on passe ici
    """but : etude de l'enfouissement pour tous les conformations d'une proteine 
    input : dico de proteine (avec tous ces conformations)
    output : une tuple contenant un dictionnaire de centre de masse d'une residus (pour tous conformations) ainsi qu'un dictionnaire de valeur de l'enfouissement moyenne
    """ 
    glob={}
    glob["prot"]=[]
    for conformation in dico:
        dico_CM=CMglob(dico[conformation])
        for res in dico_CM:
            if res not in glob:
                glob[res]=[]
            glob[res].append(dico_CM[res]) #centre de masse de chaque residu de chaque conformation
        glob["prot"].append(dico_CM["prot"]) #centre de masse de la proteine
    glob["residulist"]=dico_CM["residulist"] #les residus d'une proteine dans l'ordre

    for res in glob.keys():
        if res != "residulist" and res != "prot":
            for i in range (len(glob[res])):

                glob[res][i]=(Distance(glob[res][i][0],glob[res][i][1],glob[res][i][2],glob["prot"][i][0],glob["prot"][i][1],glob["prot"][i][2]))
          
    
        #2. Calcule de la distance moyenne
    moy={}
    del glob["prot"]
    for res in glob.keys():
        
        if res != "residulist":
            moy[res]= sum(glob[res])/len(glob[res])
            
    #~ print moy
    
    #~ final=[glob,moy] #Liste qui contient le dico de distance et celui de distance moyenne
    
    return (glob,moy)
    

#calcul de RMSD local
def RMSDlocal(dico1,dico2): #Calcule pour chaque residu le RMSD et renvoie le RMSD moyen pour chaque residu
    """but : calculer le RMSD moyen de chaque residu
    input : un dico de proteine et un dico de structure d'origine
    output : un dictionnaire contenant nb residu comme cle, et RMSD local comme valeur, et une cle residulist contenant en ordre les residus
    """ 
    #les dicos comme arguments
    flex_res={}
    flex_res["residulist"]=[]
    for conformation in dico1:
        for chain in dico1[conformation]["chains"]:
            for res in dico1[conformation][chain]["reslist"]:
                delta=[]
                for atom in dico1[conformation][chain][res]["atomlist"]:
                    x1=float(dico1[conformation][chain][res][atom]["x"])
                    x2=float(dico2["0"][chain][res][atom]["x"])
                    y1=float(dico1[conformation][chain][res][atom]["y"])
                    y2=float(dico2["0"][chain][res][atom]["y"])
                    z1=float(dico1[conformation][chain][res][atom]["z"])
                    z2=float(dico2["0"][chain][res][atom]["z"])
                    delta.append(Distance(x1,y1,z1,x2,y2,z2))#l'ensemble de distances des atomes d'un residu
                if res not in flex_res.keys():
                    flex_res[res]=[]
                flex_res[res].append(RMSD(delta)) # pour un residu donne, on ajoute son rmsd dans chaque conformation
                if conformation=="0":
                    flex_res["residulist"].append(dico1[conformation][chain][res]["resname"]) #les residus dans une prot en ordre

    for res in flex_res:
        if res != "residulist":
            flex_res[res]=sum(flex_res[res])/len(flex_res[res]) #pour un residu, on fait la moyenne de rmsd pour toutes les conformations
    return flex_res #dico ayant ResSeq comme cle, et RMSD local comme valeur, et une cle residulist contenant en ordre les residus


def giration(dico):#c'est le dico[key] qu'on passe ici
    dico_CM=CMglob(dico)
    list_dist=[]
    for res in dico_CM:
        if res != "residulist":
            list_dist.append(Distance(dico_CM["prot"][0],dico_CM["prot"][1],dico_CM["prot"][2],dico_CM[res][0],dico_CM[res][1],dico_CM[res][2]))
    return max(list_dist)
#################################################################
#Permet de representer des graphs
#ordonnee et abscisse sont les elements que lon veut representer : peuevnt etre des listes ou des dicos
#ordonne2 est une duexieme ordonneee (permet de superposer 2 graphs). Si on veut representer un graph simple on met une liste vide

#Rajouter les legendes !!!
#Et mettre dans un pdf

def graph(ordonnee,abscisse,ordonne2,titre,type_graph):

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

        if type_graph =="point":
            plt.scatter(x,y,c='red')
            plt.scatter(x,y2,c='blue')
            #~ axes = plt.gca()

            #~ plt.legend(['RMSD','Giration'])
        if type_graph=="line":
            plt.plot(x,y,c='red')
            plt.plot(x,y2,c='blue')


     #On affiche le graph
    plt.show()


##########################################################################
def writefile_glob(dico,dRMSD,dGiration):
    out=open("%s/PythonProgResults/GlobalAnalysis_%s"%(path,os.path.basename(fichier)),"w")

    #ecrire dans l'ordre les resultats de la comparaison des conformations avec la structure d'origine
    out.write("\nConformation \t RMSD results \t\t Giration\n")
    for i in range(len(dico)):
        out.write("\n%s: \t\t %.12f \t %.12f"%(i,dRMSD["%s"%i],dGiration["%s"%i]))

    out.close()


def writefile_local(dRMSD_moy,dEnf):
    out=open("%s/PythonProgResults/LocalAnalysis_%s"%(path,os.path.basename(fichier)),"w")

    out.write("Residue number \t\t Residue \t Mean RMSD \t\t Residue depth \n")
    for i in range(1,len(dRMSD_moy)):
        out.write("%s \t\t\t %s \t\t %.12f \t %.12f \n"%(i,dRMSD_moy["residulist"][i-1],dRMSD_moy["%s"%i],dEnf["%s"%i]))
    out.close()
#######################################################################################
def Global(fichier):

    print "Parsing:",fichier
    dico=ParsingPDB(fichier) #1/ Parse le fichier pdb
    list_temps=Temps(fichier)

    #2.calcul RMSD de chaque conformation par rapport a la structure d'origine
    #3. et aussi calcul du rayon de giration de chaque conformation

    dico_RMSD={}
    dico_Giration={}

    for key in dico:

        list_delta=[]
        dico_Giration[key]=giration(dico[key])

        for chain in dico[key]["chains"]:

            for res in dico[key][chain]["reslist"]:
                list_delta.append(Distance(float(dico[key][chain][res]['CA']['x']),float(dico[key][chain][res]['CA']['y']),float(dico[key][chain][res]['CA']['z']),float(dico_ref['0'][chain][res]['CA']['x']),float(dico_ref['0'][chain][res]['CA']['y']),float(dico_ref['0'][chain][res]['CA']['z']))) # compare dist entre les Ca

        dico_RMSD[key]=RMSD(list_delta)

    writefile_glob(dico,dico_RMSD,dico_Giration)

    list_conformation=sorted(int(i) for i in dico.keys()) #numero de conformation trie

    #Variation du RMSD en fonction du temps
    title='Evolution du RMSD en fonction du temps'
    l2=[]
    not l2

    graph(dico_RMSD,list_temps,l2,title,'line')

    #Variation du rayon de giration en fonction du temps
    title='Evolution du rayon de giration en fonction du temps'
    graph(dico_Giration,list_temps,l2,title,'line')

    #Variation du RMSD et Giration en fonction de la conformation
    graph(dico_RMSD,list_conformation,dico_Giration,title,"line")
    title='Variation RMSD/Giration en fonction de la conformation'


def Local(fichier):
    print "Parsing:",fichier
    dico=ParsingPDB(fichier) #1/ Parse le fichier pdb
    list_temps=Temps(fichier)
    dico_RMSD_moy=RMSDlocal(dico,dico_ref)      
    (dicoCM,dicoEnf)=Enfouissement(dico) #On recupere enfouissement pour chaque res selon les conformations et lenfouissement moyen
    #graph(dicoCM,list_temps,[],"essai","line")
    writefile_local(dico_RMSD_moy,dicoEnf)          
###############################################################################################

if __name__ == '__main__':

    ########################################################
    #                   Arguments
    ########################################################

    parser = argparse.ArgumentParser()
    parser.add_argument("-p", help="Path to directory where pdb files are stored")
    parser.add_argument("-a", help="Type of analysis: global or local or both")
    args = parser.parse_args()
    if args.p == None or args.a == None: # si l'un des arguments est vide
        parser.print_help()
        parser.exit()
    path=args.p
    analyse=args.a

    path=os.path.abspath(path)

    ########################################################
    #               Etapes preliminaires
    ########################################################

    #on supprime le dossier contenant les resultats de l'execution precedente (s'il existe) afin d'eviter les chevauchements
    #le dossier "PythonProgResults" sert a stocker les fichiers sortis

    try:
        shutil.rmtree("%s/PythonProgResults"%path)
    except OSError:
        pass

    os.chdir(path)
    fichiers=glob.iglob("*.pdb") # si le path est vers un dossier, on lira tous les fichiers pdb

    os.mkdir("%s/PythonProgResults"%path)

    #########################################################
    #                       Main
    #########################################################
    dico_ref=ParsingPDB("start_prot_only.pdb")#dictionnaire de reference de la structure d'origine
    if (analyse == "global"): #si vous voulez seulement une analyse globale
        for fichier in fichiers:
            if fichier!="start_prot_only.pdb":
                Global(fichier)

    elif (analyse == "local"): #si vous voulez seulement une analyse locale
        for fichier in fichiers:
            if fichier!="start_prot_only.pdb":
                Local(fichier)
    else : #si vous voulez a la fois une analyse globale et une analyse locale
        for fichier in fichiers:
            if fichier!="start_prot_only.pdb":
                Global(fichier)
                Local(fichier)

            #Analyse des resultats
                #Graph1: Distance en fonction du numero des residus
                #Idee : superposer aussi les valeurs moyennes
                #Mettre une legende pour dire a quoi correspondent les chiffres de laxe des x
            #~ y=np.array(list_dist)
            #~ x=np.array(list_res)
            #~ plt.scatter(x,y)
            #~ plt.show()

                #Graph2: RMSD moyen en fonction du temps

            #~ y2=np.array()
            #~ x=np.array(list_temps)
            #~ plt.plot(x,y2)
            #~ plt.show()

                #~ #Graph3: CM moyen en fonction du temps
            #~ y3.np.array()
            #~ x.np.array(list_temps)
            #~ plt.plot(x,y3)
            #~ plt.show()








