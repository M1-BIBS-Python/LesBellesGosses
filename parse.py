#!/usr/bin/env python
#-*- coding : utf8 -*-

"""
Author : LesBellesGosses
Date : 12/04/2017
Description : Projet Barstar
"""
from math import sqrt
#~ import matplotlib.pyplot as plt
#~ import numpy as np

def ParsingPDB (pdbFile):#fonction qui parse un fichier pdb

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

def Temps (pdbFile): #Fonction qui lit les premieres lignes du fichier pdb et mets dans une liste les valeurs du temps
    temps = []

    infile = open (pdbFile, "r")
    lines = infile.readlines()

    for line in lines :
        if line[:5:] == "TITLE":
            timet=line[65:80].strip()
            temps.append(timet)
    #print temps
    infile.close()
    return(temps)


def CM(listx,listy,listz):
    x=sum(listx)/float(len(listx))
    y=sum(listy)/float(len(listy))
    z=sum(listz)/float(len(listy))
    coords=[x,y,z]
    return coords


def Distance(x1,y1,z1,x2,y2,z2): #Calcule la distance entre deux points
    return(sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2))


def RMSD(list_delta): #calcul de RMSD
    #List_delta est une liste de distances
    distcarre=[]
    for delta in list_delta:
        distcarre.append(delta**2)
    return (sqrt((sum(distcarre))/float(len(list_delta))))


def CMglob(dico):#c'est le dico[key] qu'on passe ici
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
            glob["%s"%res]=CM(listx,listy,listz)
            glob["residulist"].append(dico[chain][res]["resname"])
    glob["prot"]=CM(globx,globy,globz)
    return glob # dico contient le CM de chaque residu et le CM de la prot


def RMSDlocal(dico1,dico2): #les dico[key] comme arguments
    flex_res={}
    for chain in dico1["chains"]:
        for res in dico1[chain]["reslist"]:
            delta=[]
            for atom in dico1[chain][res]["atomlist"]:
                x1=float(dico1[chain][res][atom]["x"])
                x2=float(dico2[chain][res][atom]["x"])
                y1=float(dico1[chain][res][atom]["y"])
                y2=float(dico2[chain][res][atom]["y"])
                z1=float(dico1[chain][res][atom]["z"])
                z2=float(dico2[chain][res][atom]["z"])
                delta.append(Distance(x1,y1,z1,x2,y2,z2))
            flex_res=RMSD(delta)
    return flex_res #dico ayant ResSeq comme cle, et RMSD local comme valeur


def giration(dico):#c'est le dico[key] qu'on passe ici
    dico_CM=CMglob(dico)
    list_dist=[]
    for res in dico_CM.keys():
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
        for i in range(len(ordonnee)): #on parcours le dictionnaire
            ordo.append(ordonnee["%s"%i]) #Et on range dans une liste les differentes valeurs contenues dans le dico, dans lordre dans lesquelles on les a trouve


    if type(ordonne2) is list: #Si cest une liste
        ordo2=ordonne2
    else:
        for i in range(len(ordonne2)):
            ordo.append(ordonne2["%s"%i])


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
        if type_graph=="point":
            plt.plot(x,y,c='red')
            plt.plot(x,y2,c='blue')


     #On affiche le graph
    plt.show()


##########################################################################"
def writefile_glob():
    out=open("%s/PythonProgResults/GlobalAnalysis_%s"%(path,os.path.basename(fichier)),"w")

    #ecrire dans l'ordre les resultats de la comparaison des conformations avec la structure d'origine
    out.write("\nConformation \t RMSD results \t\t Giration\n")
    for i in range(len(dico)):
        out.write("\n%s: \t\t %.12f \t %.12f"%(i,dico_RMSD["%s"%i],dico_Giration["%s"%i]))

    out.close()


def writefile_local():
    out=open("%s/PythonProgResults/LocalAnalysis_%s"%(path,os.path.basename(fichier)),"w")

    for i in range(len(dico_dist)):
        out.write("MODEL %s\n"%i)
        out.write("Residu sequence number \t Residu \t Distance between CM of residu and CM of protein \t RMSD \n")
        for j in range(1,len(dico_dist["%s"%i])):
            out.write("%s \t\t\t %s \t\t\t %.12f \t\t\t\t %.12f\n"%(j,dico_dist["%s"%i]["residulist"][j-1],dico_dist["%s"%i]["%s"%j],dico_RMSD["%s"%i]["%s"%j]))
    out.close()


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

    #savoir si le chemin s'agit d'un chemin vers dossier ou fichier
    est_fichier=False
    path=os.path.abspath(path)
    if os.path.isfile(path):
        fichiers=glob.iglob("%s"%path)
        path=os.path.dirname(path)
        est_fichier=True

    ########################################################
    #               Etapes preliminaires
    ########################################################

    #on supprime le dossier contenant les resultats de l'execution precedente (s'il existe) afin d'eviter les chevauchements
    #le dossier "PythonProgResults" sert a stocker les fichiers sortis

    try:
        shutil.rmtree("%s/PythonProgResults"%path)
    except OSError:
        pass

    if not est_fichier:
        os.chdir(path)
        fichiers=glob.iglob("*.pdb") # si le path est vers un dossier, on lira tous les fichiers pdb

    os.mkdir("%s/PythonProgResults"%path)

    #########################################################
    #                       Main
    #########################################################

    if (analyse == "global"): #si vous voulez seulement une analyse globale
        for fichier in fichiers:
            print "Parsing:",fichier
            dico=ParsingPDB(fichier) #1/ Parse le fichier pdb
            list_temps=Temps(fichier)

            #2.calcul RMSD de chaque conformation par rapport a la structure d'origine
            #3. et aussi calcul du rayon de giration de chaque conformation

            dico_RMSD={}
            dico_Giration={}
            list_conformation=[] #Liste qui va contenir le numero de chaque conformation

            for key in dico:

                list_delta=[]
                dico_Giration[key]=giration(dico[key])

                list_conformation.append(float(key))               #On recupere le numero de la conformation
                for chain in dico[key]["chains"]:

                     for res in dico[key][chain]["reslist"]:
                         list_delta.append(Distance(float(dico[key][chain][res]['CA']['x']),float(dico[key][chain][res]['CA']['y']),float(dico[key][chain][res]['CA']['z']),float(dico['0'][chain][res]['CA']['x']),float(dico['0'][chain][res]['CA']['y']),float(dico['0'][chain][res]['CA']['z']))) # compare dist entre les Ca

                dico_RMSD[key]=RMSD(list_delta)


            writefile_glob()

            #~ #Variation du RMSD en fonction du temps
            #~ title='Evolution du RMSD en fonction du temps'
            #~ l2=[]
            #~ not l2

            #~ graph(dico_RMSD,list_temps,l2,title,'line')

            #~ #Variation du rayon de giration en fonction du temps
            #~ title='Evolution du rayon de giration en fonction du temps'
            #~ graph(dico_Giration,list_temps,l2,title,'line')


            #~ #graph representant RMSD/Giration en fonction de la conformation
            #~ title='Variation RMSD/GIRATION en fonction de la conformation'
            #~ graph(dico_RMSD,list_conformation,dico_Giration,title,"point")


    elif (analyse == "local"): #si vous voulez seulement une analyse locale
        for fichier in fichiers:
            print "Parsing:",fichier
            dico=ParsingPDB(fichier)
            list_temps=Temps(fichier)
            
            
            #dist de chaque residu de chaque conformation par rapport au CM de la prot
            #et calcul de RMSD local
            dico_dist={}
            dico_RMSD={} #Va contenir les RMSD de chaque residus par rapport a la position dans la structure
            dico_RMSD_moyen={} #Va contenir pour chaque residu le RMSD moyen
            
            for key in dico: #pour chaque conformation
                dico_Enfouissement=CMglob(dico[key])
                dico_RMSD[key]=RMSDlocal(dico["0"],dico[key]) #comparaison de deux conformations (dans la fct, comparaison de chaque residu des deux conformations)
                RMSDlist=RMSDlocal(dico["0"],dico[key]) #comparaison de deux conformations (dans la fct, comparaison de chaque residu des deux conformations)
				
				#####################Calcul du RMSD moyen################################
                #~ for chain in dico[key]["chains"]:
                     #~ for res in dico[key][chain]["reslist"]:
						 #~ print res
						 #~ residu=dico[key][chain][res]["resname"] #on recupere le nom du residu
				
						 #On met dans un dico le nom du residu et son RMSD
						 #~ if residu not in dico_RMSD_moyen.keys(): #si le residus nest pas deja une cle
							 #~ dico_RMSD_moyen[residu]=dico_RMSD[key]
							
						 #~ else: #sinon on ajoute la valeur a la cle
							 #~ for valeur in dico_RMSD_moyen[residu]:
								#~ print valeur
								#~ dico_RMSD_moyen[residu]=(residu+dico_RMSD[key])/2 #On calcule la moyenne au fur et a mesure
                #~ print dico_RMSD_moyen
				#########################################################################		 
						
							 
                for res in dico_Enfouissement:
                    
                    if res != "prot" and res != "residulist":
                        dico_dist[key][res]=Distance(dico_Enfouissement["prot"][0],dico_Enfouissement["prot"][1],dico_Enfouissement["prot"][2],dico_Enfouissement[res][0],dico_Enfouissement[res][1],dico_Enfouissement[res][2])


                    if res == "residulist":
                        dico_dist[key]["residulist"]=dico_Enfouissement["residulist"]



                    #~ for atom in dico[chain][res]["atomlist"]:
                        #~ print atom
            #RMSD moyen pour chaque residu de chaque conformation (a ajouter)

            writefile_local()
            





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

    #~ else : #si vous voulez a la fois une analyse globale et une analyse locale
        #~ for fichier in fichiers:
            #~ #do sth






