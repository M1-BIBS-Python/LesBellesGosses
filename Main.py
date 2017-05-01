#!/usr/bin/env python
#-*- coding : utf8 -*-

"""
Author : LesBellesGosses
Description : Projet Barstar
"""

import argparse,os,glob,shutil,sys
import Global as Global
import Local as Local
import StructureTools as ST

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
dico_ref=ST.ParsingPDB("start_prot_only.pdb")#dictionnaire de reference de la structure d'origine

if (analyse == "global"): #si vous voulez seulement une analyse globale
    for fichier in fichiers:
        if fichier!="start_prot_only.pdb":
            Global.Global(fichier,dico_ref,path)

elif (analyse == "local"): #si vous voulez seulement une analyse locale
    for fichier in fichiers:
        if fichier!="start_prot_only.pdb":
            Local.Local(fichier,dico_ref,path)
            
else : #si vous voulez a la fois une analyse globale et une analyse locale
    for fichier in fichiers:
        if fichier!="start_prot_only.pdb":
            Global.Global(fichier,dico_ref)
            Local.Local(fichier)









