#!/usr/bin/env python
#-*- coding : utf8 -*-

"""
Author : LesBellesGosses
Description : Projet Barstar
"""

from math import sqrt
import WriteFile as writefile
import StructureTools as ST


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
                list_delta.append(ST.Distance(float(dico[key][chain][res]['CA']['x']),float(dico[key][chain][res]['CA']['y']),float(dico[key][chain][res]['CA']['z']),float(dico_ref['0'][chain][res]['CA']['x']),float(dico_ref['0'][chain][res]['CA']['y']),float(dico_ref['0'][chain][res]['CA']['z']))) # compare dist entre les Ca
        dico_RMSD[key]=ST.RMSD(list_delta)
    return dico_RMSD

 
#calcul de giration
def giration(dico):#c'est le dico[key] qu'on passe ici
    """but : calculer le rayon de giration d'une proteine
    input : un dico de proteine d'une conformation donnee
    output : la valeur du rayon de giration
    """
    dico_Giration={}
    print dico.keys()
    for key in dico:
        dico_CM=ST.CMglob(dico[key])
        print dico_CM
        list_dist=[]
        for res in dico_CM:
            if res != "residulist":
                list_dist.append(ST.Distance(dico_CM["prot"][0],dico_CM["prot"][1],dico_CM["prot"][2],dico_CM[res][0],dico_CM[res][1],dico_CM[res][2]))
        dico_Giration[key]=max(list_dist)
    return dico_Giration

    
def Global(fichier,dico_ref,path):
    """but: Analyse des changements conformationnels globaux de la proteine
    input:lensemble des fichiers de donnees presents dans le repertoire
    """
    print "Parsing:",fichier
    dico=ST.ParsingPDB(fichier) #1/ Parse le fichier pdb
    list_temps=ST.Temps(fichier)

    #2.calcul RMSD de chaque conformation par rapport a la structure d'origine
    #3. et aussi calcul du rayon de giration de chaque conformation
    dico_RMSD=RMSDglobal(dico,dico_ref)
    dico_Giration=giration(dico)
    
    writefile.writefile_glob(dico,dico_RMSD,dico_Giration,fichier,path)

    list_conformation=sorted(int(i) for i in dico.keys()) #numero de conformation trie

    #Variation du RMSD en fonction du temps
    title='Evolution du RMSD en fonction du temps'
    l2=[]
    not l2
    ST.graph(dico_RMSD,list_temps,l2,title,"RMSD","temps (ps)",'line')

    #Variation du rayon de giration en fonction du temps
    title='Evolution du rayon de giration en fonction du temps'
    ST.graph(dico_Giration,list_temps,l2,title,"rayon de giration","temps (ps)",'line')

    #Variation du RMSD et Giration en fonction de la conformation
    title='Variation RMSD/Giration en fonction de la conformation'
    ST.graph(dico_RMSD,list_conformation,dico_Giration,title,"[RMSD (rouge),Giration (bleu)]","conformation","line")

