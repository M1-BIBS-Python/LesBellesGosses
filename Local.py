#!/usr/bin/env python
#-*- coding : utf8 -*-

"""
Author : LesBellesGosses
Description : Projet Barstar
Analyse locale : identification des regions flexibles
"""

import WriteFile as writefile
import StructureTools as ST


#Enfouissement de chaque residu dans une proteine
def Enfouissement(dico,dico_ref):
    """but : etude de l'enfouissement pour toutes les conformations d'une proteine
    input : dictionnaire de proteine et dictionnaire de proteine de la structure d'origine
    output : un dictionnaire de valeur de l'enfouissement moyenne
    """
    dGlob={}
    dGlob_ref=ST.CMglob(dico_ref["0"])
    for conformation in dico:
        dico_CM=ST.CMglob(dico[conformation])
        for res in dico_CM:
            if res not in dGlob:
                dGlob[res]=[]
            dGlob[res].append(dico_CM[res]) #centre de masse de chaque residu de chaque conformation
    dGlob["residulist"]=dico_CM["residulist"] #les residus d'une proteine dans l'ordre


	#1. Calcul de l'enfouissement pour chaque residu pour chaque conformation
    for res in dGlob.keys():
        if res != "residulist":
            for i in range (len(dGlob[res])):
                dGlob[res][i]=(ST.Distance(dGlob[res][i][0],dGlob[res][i][1],dGlob[res][i][2],dGlob_ref["prot"][0],dGlob_ref["prot"][1],dGlob_ref["prot"][2]))

    #2. Calcul de l'enfouissement moyen
    dMoy={}
    for res in dGlob.keys(): 

        if res != "residulist" and res!="prot":
            dMoy[res]= sum(dGlob[res])/len(dGlob[res])

    return dMoy
    
    
#Calcul de RMSD local
def RMSDlocal(dico1,ref): #Calcule pour chaque residu le RMSD et renvoie le RMSD moyen pour chaque residu
    """but : calculer le RMSD moyen de chaque residu
    input : un dictionnaire de proteine et un dictionnaire de structure d'origine
    output : un dictionnaire contenant pour chaque residu les RMSD de chaque conformation et un dictionnaire de RMSD moyen (pour chaque residu, un RMSD moyen)
    """
    #les dictionnaires comme arguments
    dFlex_res={} #dictionnaire ayant nom de residu comme cle et son RMSD de chaque conformation comme valeur
    dFlex_res["residulist"]=[]
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
                 
                if res not in dFlex_res.keys():#si res n'est pas deja une cle, on cree une nouvelle cle
                    dFlex_res[res]=[]
                dFlex_res[res].append(ST.RMSD(delta)) # pour un residu donne, on ajoute son rmsd de chaque conformation
                
                if conformation==0:
                    dFlex_res["residulist"].append(dico1["%s"%conformation][chain][res]["resname"])

    for res in dFlex_res.keys():
        if res !="residulist":
            dico_moy[res]=sum(dFlex_res[res])/len(dFlex_res[res]) #Pour chaque residu, on calcule le RMSD moyen

    return (dFlex_res,dico_moy)

def Local(fichier,dico_ref,path):
    """but :Analyse des changements conformationnels locaux de la proteine
    input: les fichiers que l'on veut analyser, le fichier de la structure de reference et le chemin
    """
    print "Parsing:",fichier
    dico=ST.ParsingPDB(fichier)
    list_temps=ST.Temps(fichier)
    
    (dico_RMSD,dico_RMSD_moy)=RMSDlocal(dico,dico_ref) #1. Calcul du RMSD moyen pour chaque residu
    dicoEnf_moy=Enfouissement(dico,dico_ref) #2. Calcul d'enfouissement moyen
    
    #On classe en 2 classes les residus selon les valeurs de RMSD moyen, en 2 classes pour les valeurs d'enfouissement
    dclasse=ST.createClass(dico_RMSD_moy,max(dico_RMSD_moy.values())+min(dico_RMSD_moy.values()),2)
    dclasseEnf=ST.createClass(dicoEnf_moy,max(dicoEnf_moy.values())+min(dicoEnf_moy.values()),2) 

    writefile.writefile_local(dico_RMSD["residulist"],dico_RMSD_moy,dicoEnf_moy,dclasse,dclasseEnf,fichier,path)
   
  

    ##########################Analyses graphiques###########################################
    
    #Enfouissement des residus : Enfouissement moyen en fonction du numero de residus
    num=[] #Liste qui contient le numero de residu
    l2=[]

    for cle in range(len(dicoEnf_moy)):
        num.append(cle)
    del num[0]
    ST.graph(dicoEnf_moy,num,l2,"Analyse locale : Enfouissement moyen en fonction du numero de residus","Enfouissement (A) ","residu")

    #Identification des residus presents dans les regions flexibles : RMSD moyen en fonction du numero de residu
    num=[]
    l2=[]
    
    for cle in range(len(dico_RMSD_moy)):
        num.append(cle)
    del num[0]
    ST.graph(dico_RMSD_moy,num,l2,"Analyse locale : RMSD moyen en fonction du numero de residus","RMSDmoyen(A)","residu")
    
    #Comparaison enfouissement des residus avec le RMSD
    ST.graph(dico_RMSD_moy,num,dicoEnf_moy,"Comparaison de Enfouissement et RMSD moyen en fonction du residu","[RMSD moyen(rouge),Enfouissement(bleu)]","residu")
    
    #Regarder l'evolution de quelques residus au cours du temps (cf residu pris dans la publi)
    #recuperer pour un residu toutes les valeurs du RMSD et les representer en fonction du temps
    l2=[]
   
    for cle in range(1,len(dico_RMSD)+1):
        
        if cle==72:
            ST.graph(dico_RMSD["%s"%cle],list_temps,l2,"Analyse locale: Evolution du RMSD du residu 72 en fonction du temps","RMSD (A)","temps(ps)")
        if cle==42:
            ST.graph(dico_RMSD["%s"%cle],list_temps,l2,"Analyse locale: Evolution du RMSD du residu 42 en fonction du temps","RMSD (A)","temps(ps)")
        
        if cle==39:
            ST.graph(dico_RMSD["%s"%cle],list_temps,l2,"Analyse locale: Evolution du RMSD du residu 39 en fonction du temps","RMSD (A)","temps(ps)")
        if cle==76:
            ST.graph(dico_RMSD["%s"%cle],list_temps,l2,"Analyse locale: Evolution du RMSD du residu 76(GLU) en fonction du temps","RMSD (A)","temps(ps)")
            
        if cle==38:
            ST.graph(dico_RMSD["%s"%cle],list_temps,l2,"Analyse locale: Evolution du RMSD du residu 38(LYS) en fonction du temps","RMSD (A)","temps(ps)")
       
       
