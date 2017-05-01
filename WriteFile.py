#!/usr/bin/env python
#-*- coding : utf8 -*-

"""
Author : LesBellesGosses
Description : Projet Barstar
Ecriture dans un fichier de sortie
"""

import os

def writefile_glob(dico,dRMSD,dGiration,fichier,path):
    """but : ecrire un fichier de sortie contenant pour chaque conformation et pour la structure d'origine le rayon de giration et le RMSD correspondant
    input : le dictionaire de la proteine, le dictionnaire de RMSD et le dictionnaire du rayon de giration,le nom du fichier pdb et le chemin du repertoire
    output: un fichier texte qui possede pour chaque conformation, le RMSD et le rayon de giration
    """
    out=open("%s/PythonProgResults/GlobalAnalysis_%s"%(path,os.path.basename(fichier)),"w")
    out.write("\nConformation \t RMSD results \t\t Giration\n")
    for i in range(len(dico)):
        out.write("\n%s: \t\t %.12f \t %.12f"%(i,dRMSD["%s"%i],dGiration["%s"%i]))
    out.close()


def writefile_local(residulist,dRMSD_moy,dEnf,dclasse,dclasseEnf,fichier,path):
    """but : ecrire un fichier contenant pour chaque residu le RMSD moyen, ainsi que la distance moyenne de chacun des residus par rapport au centre de masse (enfouissement)
    input: le nom des residus, dictionnaire de RMSD moyen et dictionnaire de l'enfouissement, le nom du fichier et le chemin du repertoire
    dclasse et dclasseEnf sont ajoutes pour identifier par * les residus probablement les plus enfouis ou les plus flexibles
    output: un fichier texte
    """
    out=open("%s/PythonProgResults/LocalAnalysis_%s"%(path,os.path.basename(fichier)),"w")
    out.write("Residue number \t\t Residue \t Mean RMSD \t\t\t Residue depth \n")
    for i in range(1,len(dRMSD_moy)+1):
        out.write("%s \t\t\t %s \t\t %.12f"%(i,residulist[i-1],dRMSD_moy["%s"%i]))
        if dclasse["%s"%i]==1: #identifier les residus les plus flexibles (ils sont de classe 1)
            out.write(" *")
        if dclasseEnf["%s"%i]==max(dclasseEnf.values()): #identifier les residus les plus enfouis (= ceux qui ont un numero de classe maximal)
            out.write(" \t\t %.12f *\n"%dEnf["%s"%i])
        else:
            out.write(" \t\t %.12f \n"%dEnf["%s"%i])
        

    out.close()

