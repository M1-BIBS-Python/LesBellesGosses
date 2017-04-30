#!/usr/bin/env python
#-*- coding : utf8 -*-

"""
Author : LesBellesGosses
Description : Projet Barstar
Analyse locale : identification des regions flexibles
"""
from math import sqrt
import WriteFile as writefile
import StructureTools as ST


#Enfouissement de chaque residu dans une proteine (pour l'analyse locale)
def Enfouissement(dico,dico_ref):#c'est le dico qu'on passe ici
    """but : etude de l'enfouissement pour tous les conformations d'une proteine
    input : dico de proteine et dico de proteine de la structure d'origine
    output : une tuple contenant un dictionnaire de centre de masse d'une residus (pour tous conformations) ainsi qu'un dictionnaire de valeur de l'enfouissement moyenne
    """
    glob={}
    glob_ref=ST.CMglob(dico_ref["0"])
    for conformation in dico:
        dico_CM=ST.CMglob(dico[conformation])
        for res in dico_CM:
            if res not in glob:
                glob[res]=[]
            glob[res].append(dico_CM[res]) #centre de masse de chaque residu de chaque conformation
    glob["residulist"]=dico_CM["residulist"] #les residus d'une proteine dans l'ordre


	#1. Calcul de lenfouissement pour chaque residu pour chaque conformation
    for res in glob.keys():
        if res != "residulist":
            for i in range (len(glob[res])):
                glob[res][i]=(ST.Distance(glob[res][i][0],glob[res][i][1],glob[res][i][2],glob_ref["prot"][0],glob_ref["prot"][1],glob_ref["prot"][2]))

    #2. Calcule de lenfouissement moyen
    moy={}
    for res in glob.keys(): 

        if res != "residulist" and res!="prot":
            moy[res]= sum(glob[res])/len(glob[res])

    return (glob,moy)
    
    
#calcul de RMSD local
def RMSDlocal(dico1,ref): #Calcule pour chaque residu le RMSD et renvoie le RMSD moyen pour chaque residu
    """but : calculer le RMSD moyen de chaque residu
    input : un dico de proteine et un dico de structure d'origine
    output : un dico contenant pour chaque residu les RMSD de chaque conformation et un dico de RMSD moyen (pour chaque residu, un RMSD moyen)
    """
    #les dicos comme arguments
    flex_res={} #dico ayant nom de residu comme cle et son RMSD de chaque conformation comme valeur
    flex_res["residulist"]=[]
    dico_moy={}
    
    for conformation in range(len(dico1)):
        for chain in dico1["%s"%conformation]["chains"]:
            for res in dico1["%s"%conformation][chain]["reslist"]:
                delta=[]
                for atom in dico1["%s"%conformation][chain][res]["atomlist"]: #on recupere les coordonnees
                    x1=float(dico1["%s"%conformation][chain][res][atom]["x"])
                    x2=float(ref["0"][chain][res][atom]["x"])
                    y1=float(dico1["%s"%conformation][chain][res][atom]["y"])
                    y2=float(ref["0"][chain][res][atom]["y"])
                    z1=float(dico1["%s"%conformation][chain][res][atom]["z"])
                    z2=float(ref["0"][chain][res][atom]["z"])
                    
                    delta.append(ST.Distance(x1,y1,z1,x2,y2,z2))#l'ensemble de distances des atomes d'un residu
                 
                if res not in flex_res.keys():#si res nest pas deja une cle, on cree une nouvelle cle
                    flex_res[res]=[]
                flex_res[res].append(ST.RMSD(delta)) # pour un residu donne, on ajoute son rmsd dans chaque conformation
                
                if conformation==0:
                    flex_res["residulist"].append(dico1["%s"%conformation][chain][res]["resname"])

    for res in flex_res.keys():
        if res !="residulist":
            dico_moy[res]=sum(flex_res[res])/len(flex_res[res]) #Pour chaque residu, on calcule le RMSD moyen


    return (flex_res,dico_moy)

def Local(fichier,dico_ref,path):
    """but :Analyse des changements conformationnels locaux de la proteine
    input: les fichiers que lon veut analyser, le fichier de la structure de reference et le chemin
    """
    print "Parsing:",fichier
    dico=ST.ParsingPDB(fichier)
    list_temps=ST.Temps(fichier)
    
    (dico_RMSD,dico_RMSD_moy)=RMSDlocal(dico,dico_ref) #2. Calcul du RMSD moyen pour chaque residu
    (dicoCM,dicoEnf)=Enfouissement(dico,dico_ref) #3. enfouissement pour chaque res selon les conformations et lenfouissement moyen
    
    #On classe les residus selon leur valeurs de RMSD moyens : on veut 2 classes (1 pour les residus avec RMSD inferieur au seuil
    # et une autre pour les residus avec RMSD superieur au seuil
    
    dclasse=ST.createClass(dico_RMSD_moy,max(dico_RMSD_moy.values())+min(dico_RMSD_moy.values()),2)
    dclasseEnf=ST.createClass(dicoEnf,max(dicoEnf.values())+min(dicoEnf.values()),3) #j'avais essaye 2 classes mais il me semble qu'il y a trop de residus significatifs, du coup 3, apres il faut verifier avec pymol

    writefile.writefile_local(dico_RMSD["residulist"],dico_RMSD_moy,dicoEnf,dclasse,dclasseEnf,fichier,path)
   
  

    ##########################Analyses graphiques###########################################
    
    #Enfouissement des residus : Enfouissement moyen en fonction du numero de residus
    num=[] #Liste qui contient le numero de residu
    l2=[]
    not l2

    for cle in range(len(dicoEnf)):
        num.append(cle)
    del num[0]
    ST.graph(dicoEnf,num,l2,"Enfouissement moyen en fonction du numero de res","Enfouissement (A) ","residu")

    #Identification des residus presents dans les regions flexibles : RMSD moyen en fonction du numero de residu
    num=[]
    l2=[]
    
    for cle in range(len(dico_RMSD_moy)):
        num.append(cle)
    
    ST.graph(dico_RMSD_moy.values(),num,l2,"RMSD moyen en fonction du numero de residus","RMSDmoyen(A)","residu")
    
    #Comparaison enfouissement des residus avec le RMSD
    ST.graph(dico_RMSD_moy.values(),num,dicoEnf.values(),"Comparaison de Enfouissement et RMSD moyen en fonction du residu","[RMSDmoyen(rouge),Enfouissement(bleu)]","residu")
    
    #Regarder levolution de quelques residus au cours du temps (cf residu pris dans la publi)
    #recuperer pour un residu toutes les valeurs du RMSD et les representer en fonction du temps
    l2=[]
   
    for cle in range(1,len(dico_RMSD)+1):
        if cle==76:
            ST.graph(dico_RMSD["%s"%cle],list_temps,l2,"Evolution du RMSD du residu glu76 en fonction du temps","RMSD (A)","temps(ps)")
        if cle==39:
            ST.graph(dico_RMSD["%s"%cle],list_temps,l2,"Evolution du RMSD du residu asp 39 en fonction du temps","RMSD (A)","temps(ps)")
        if cle==15:
            ST.graph(dico_RMSD["%s"%cle],list_temps,l2,"Evolution du RMSD du residu Lys 15 en fonction du temps","RMSD (A)","temps(ps)")
        if cle==17:
            ST.graph(dico_RMSD["%s"%cle],list_temps,l2,"Evolution du RMSD du residu arg 17 en fonction du temps","RMSD (A)","temps(ps)")
        if cle==68:
            ST.graph(dico_RMSD["%s"%cle],list_temps,l2,"Evolution du RMSD du residu arg 68 en fonction du temps","RMSD (A)","temps(ps)")
        if cle==45:
            ST.graph(dico_RMSD["%s"%cle],list_temps,l2,"Evolution du RMSD du residu arg 45 en fonction du temps","RMSD (A)","temps(ps)")
            

