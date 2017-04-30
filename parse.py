#!/usr/bin/env python
#-*- coding : utf8 -*-

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


#calcul de RMSD global
def RMSDglobal(dico,dico_ref):
    """but : Calculer le RMSD entre chaque conformation et la proteine de reference
    input : dico de la proteine et dico de la proteine de la structure dorigine
    output : dico contenant pour chaque res de chaque conformation le RMSD
    """
    dico_RMSD={}
    for key in dico:
        list_delta=[]
        for chain in dico[key]["chains"]:
            for res in dico[key][chain]["reslist"]:
                list_delta.append(Distance(float(dico[key][chain][res]['CA']['x']),float(dico[key][chain][res]['CA']['y']),float(dico[key][chain][res]['CA']['z']),float(dico_ref['0'][chain][res]['CA']['x']),float(dico_ref['0'][chain][res]['CA']['y']),float(dico_ref['0'][chain][res]['CA']['z']))) # compare dist entre les Ca
        dico_RMSD[key]=RMSD(list_delta)
    return dico_RMSD


#Enfouissement de chaque residu dans une proteine (pour l'analyse locale)
def Enfouissement(dico,dico_ref):#c'est le dico qu'on passe ici
    """but : etude de l'enfouissement pour tous les conformations d'une proteine
    input : dico de proteine et dico de proteine de la structure d'origine
    output : une tuple contenant un dictionnaire de centre de masse d'une residus (pour tous conformations) ainsi qu'un dictionnaire de valeur de l'enfouissement moyenne
    """
    glob={}
    glob_ref=CMglob(dico_ref["0"])
    for conformation in dico:
        dico_CM=CMglob(dico[conformation])
        for res in dico_CM:
            if res not in glob:
                glob[res]=[]
            glob[res].append(dico_CM[res]) #centre de masse de chaque residu de chaque conformation
    glob["residulist"]=dico_CM["residulist"] #les residus d'une proteine dans l'ordre

    for res in glob.keys():

        if res != "residulist":
            for i in range (len(glob[res])):
                glob[res][i]=(Distance(glob[res][i][0],glob[res][i][1],glob[res][i][2],glob_ref["prot"][0],glob_ref["prot"][1],glob_ref["prot"][2]))

        #2. Calcule de la distance moyenne
    moy={}
    for res in sorted(glob.keys()): #pk sorted? sachant que tes cles sont str, qui fait de ton sorted n'a pas de sens si tu veux des chiffres en ordres

        if res != "residulist" and res!="prot":
            moy[res]= sum(glob[res])/len(glob[res])

    return (glob,moy)


#calcul de RMSD local
def RMSDlocal(dico1,ref): #Calcule pour chaque residu le RMSD et renvoie le RMSD moyen pour chaque residu
    """but : calculer le RMSD moyen de chaque residu
    input : un dico de proteine et un dico de structure d'origine
    output : un dictionnaire contenant une liste de residus et une liste de RMSD correspondant
    """
    #les dicos comme arguments
    flex_res={} #dico ayant nom de residu comme cle et son RMSD de chaque conformation comme valeur
    flex_res["residulist"]=[]
    dico_moy={}
    
    for conformation in range(len(dico1)):
        for chain in dico1["%s"%conformation]["chains"]:
            for res in dico1["%s"%conformation][chain]["reslist"]:
                delta=[]
                for atom in dico1["%s"%conformation][chain][res]["atomlist"]:
                    x1=float(dico1["%s"%conformation][chain][res][atom]["x"])
                    x2=float(ref["0"][chain][res][atom]["x"])
                    y1=float(dico1["%s"%conformation][chain][res][atom]["y"])
                    y2=float(ref["0"][chain][res][atom]["y"])
                    z1=float(dico1["%s"%conformation][chain][res][atom]["z"])
                    z2=float(ref["0"][chain][res][atom]["z"])
                    delta.append(Distance(x1,y1,z1,x2,y2,z2))#l'ensemble de distances des atomes d'un residu
                 
                if res not in flex_res.keys():#si res nest pas deja une cle
                    flex_res[res]=[]
                flex_res[res].append(RMSD(delta)) # pour un residu donne, on ajoute son rmsd dans chaque conformation
                
                if conformation==0:
                    flex_res["residulist"].append(dico1["%s"%conformation][chain][res]["resname"])

    for res in flex_res.keys():
        if res !="residulist":
            dico_moy[res]=sum(flex_res[res])/len(flex_res[res]) #Pour chaque residu, on calcule le RMSD moyen


    return (flex_res,dico_moy)


#calcul de giration
def giration(dico):#c'est le dico[key] qu'on passe ici
    """but : calculer le rayon de giration d'une proteine
    input : un dico de proteine d'une conformation donnee
    output : la valeur du rayon de giration
    """
    dico_Giration={}
    for key in dico:
        dico_CM=CMglob(dico[key])
        list_dist=[]
        for res in dico_CM:
            if res != "residulist":
                list_dist.append(Distance(dico_CM["prot"][0],dico_CM["prot"][1],dico_CM["prot"][2],dico_CM[res][0],dico_CM[res][1],dico_CM[res][2]))
        dico_Giration[key]=max(list_dist)
    return dico_Giration
#################################################################
#Permet de representer des graphs
#Et mettre dans un pdf

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

##########################################################################
def writefile_glob(dico,dRMSD,dGiration):
    """but : ecrire un fichier de sortie contenant pour chaque conformation et pour la structure dorigine les rayon de giration et le RMSD correspondant
    input : le dico de la proteine, le dico de RMSD et le dico du rayon de giration
    output: un fichier texte qui possede pour chaque conformation, le RMSD et le rayon de giration
    """
    out=open("%s/PythonProgResults/GlobalAnalysis_%s"%(path,os.path.basename(fichier)),"w")
    out.write("\nConformation \t RMSD results \t\t Giration\n")
    for i in range(len(dico)):
        out.write("\n%s: \t\t %.12f \t %.12f"%(i,dRMSD["%s"%i],dGiration["%s"%i]))
    out.close()


def writefile_local(residulist,dRMSD_moy,dEnf,dclasse,dclasseEnf):
    """but : ecrire un fichier contenant pour chaque residu le RMSD moyen ainsi que la distance moyenne de chacun des residus par rapport au centre de masse
    input: dico de RMSD moyen et dico de lenfouissement
    output: un fichier texte
    """
    out=open("%s/PythonProgResults/LocalAnalysis_%s"%(path,os.path.basename(fichier)),"w")
    out.write("Residue number \t\t Residue \t Mean RMSD \t\t\t Residue depth \n")
    for i in range(1,len(dRMSD_moy)+1):
        out.write("%s \t\t\t %s \t\t %.12f"%(i,residulist[i-1],dRMSD_moy["%s"%i]))
        if dclasse["%s"%i]==1:
            out.write(" *")
        if dclasseEnf["%s"%i]==1:
            out.write(" \t\t %.12f *\n"%dEnf["%s"%i])
        else:
            out.write(" \t\t %.12f \n"%dEnf["%s"%i])
        

    out.close()

#######################################################################################
def Global(fichier):
    """but: Analyse des changements conformationnels globaux de la proteine
    input:lensemble des fichiers de donnees presents dans le repertoire
    """
    print "Parsing:",fichier
    dico=ParsingPDB(fichier) #1/ Parse le fichier pdb
    list_temps=Temps(fichier)

    #2.calcul RMSD de chaque conformation par rapport a la structure d'origine
    #3. et aussi calcul du rayon de giration de chaque conformation
    dico_RMSD=RMSDglobal(dico,dico_ref)
    dico_Giration=giration(dico)
    
    writefile_glob(dico,dico_RMSD,dico_Giration)

    list_conformation=sorted(int(i) for i in dico.keys()) #numero de conformation trie

    #Variation du RMSD en fonction du temps
    title='Evolution du RMSD en fonction du temps'
    l2=[]
    not l2
    graph(dico_RMSD,list_temps,l2,title,"RMSD","temps (ps)",'line')

    #Variation du rayon de giration en fonction du temps
    title='Evolution du rayon de giration en fonction du temps'
    graph(dico_Giration,list_temps,l2,title,"rayon de giration","temps (ps)",'line')

    #Variation du RMSD et Giration en fonction de la conformation
    title='Variation RMSD/Giration en fonction de la conformation'
    graph(dico_RMSD,list_conformation,dico_Giration,title,"[RMSD (rouge),Giration (bleu)]","conformation","line")


def Local(fichier):
    """but :Analyse des changements conformationnels locaux de la proteine
    input: lensemble des fichiers de donnees presents dans le repertoire
    """
    print "Parsing:",fichier
    dico=ParsingPDB(fichier)
    list_temps=Temps(fichier)
    
    (dico_RMSD,dico_RMSD_moy)=RMSDlocal(dico,dico_ref) #2. Calcul du RMSD moyen pour chaque residu
    (dicoCM,dicoEnf)=Enfouissement(dico,dico_ref) #3. enfouissement pour chaque res selon les conformations et lenfouissement moyen
    
    #On classe les residus selon leur valeurs de RMSD moyens : on veut 2 classes (1 pour les residus avec RMSD inferieur au seuil
    # et une autre pour les residus avec RMSD superieur au seuil
    
    dclasse=createClass(dico_RMSD_moy,max(dico_RMSD_moy.values())+min(dico_RMSD_moy.values()),2)
    dclasseEnf=createClass(dicoEnf,max(dicoEnf.values())+min(dicoEnf.values()),3) #j'avais essaye 2 classes mais il me semble qu'il y a trop de residus significatifs, du coup 3, apres il faut verifier avec pymol

    writefile_local(dico_RMSD["residulist"],dico_RMSD_moy,dicoEnf,dclasse,dclasseEnf)
   
  

    ##########################Analyses graphiques###########################################
    
    #Enfouissement des residus : Enfouissement moyen en fonction du numero de residus
    num=[] #Liste qui contient le numero de residu
    l2=[]
    not l2

    for cle in range(len(dicoEnf)):
        num.append(cle)
    del num[0]
    graph(dicoEnf,num,l2,"Enfouissement moyen en fonction du numero de res","Enfouissement ","residu","line")

    #Identification des residus presents dans les regions flexibles : RMSD moyen en fonction du numero de residu
    num=[]
    l2=[]
    
    for cle in range(len(dico_RMSD_moy)):
        num.append(cle)
    graph(dico_RMSD_moy.values(),num,l2,"RMSD moyen en fonction du numero de residus","RMSDmoyen(A)","residu","line")
    
   
    #Comparaison enfouissement des residus avec le RMSD
    graph(dico_RMSD_moy.values(),num,dicoEnf.values(),"Comparaison de Enfouissement et RMSD moyen en fonction du residu","[RMSDmoyen(rouge),Enfouissement(bleu)]","residu","line")
    #Les residus dont la valeur de lenfouissement diminue ont aussi une augmentation de la valeur moyenne de leur RMSD
    
    
    #Regarder levolution de quelques residus au cours du temps (cf residu pris dans la publi)
    #recuperer pour un residu toutes les valeurs du RMSD et les representer en fonction du temps
    l2=[]
   
    for cle in range(1,len(dico_RMSD)+1):
        if cle==76:
            graph(dico_RMSD["%s"%cle],list_temps,l2,"Evolution du RMSD du residu glu76 en fonction du temps","RMSD (A)","temps(ps)","line")
        if cle==39:
            graph(dico_RMSD["%s"%cle],list_temps,l2,"Evolution du RMSD du residu asp 39 en fonction du temps","RMSD (A)","temps(ps)","line")
        if cle==15:
            graph(dico_RMSD["%s"%cle],list_temps,l2,"Evolution du RMSD du residu Lys 15 en fonction du temps","RMSD (A)","temps(ps)","line")
        if cle==17:
            graph(dico_RMSD["%s"%cle],list_temps,l2,"Evolution du RMSD du residu arg 17 en fonction du temps","RMSD (A)","temps(ps)","line")
        if cle==68:
            graph(dico_RMSD["%s"%cle],list_temps,l2,"Evolution du RMSD du residu arg 68 en fonction du temps","RMSD (A)","temps(ps)","line")
        if cle==45:
            graph(dico_RMSD["%s"%cle],list_temps,l2,"Evolution du RMSD du residu arg 45 en fonction du temps","RMSD (A)","temps(ps)","line")
            



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








