#!/usr/bin/env python
#-*- coding : utf8 -*-

import os

def writefile_glob(dico,dRMSD,dGiration,fichier,path):
    """but : ecrire un fichier de sortie contenant pour chaque conformation et pour la structure dorigine les rayon de giration et le RMSD correspondant
    input : le dico de la proteine, le dico de RMSD et le dico du rayon de giration
    output: un fichier texte qui possede pour chaque conformation, le RMSD et le rayon de giration
    """
    out=open("%s/PythonProgResults/GlobalAnalysis_%s"%(path,os.path.basename(fichier)),"w")
    out.write("\nConformation \t RMSD results \t\t Giration\n")
    for i in range(len(dico)):
        out.write("\n%s: \t\t %.12f \t %.12f"%(i,dRMSD["%s"%i],dGiration["%s"%i]))
    out.close()


def writefile_local(residulist,dRMSD_moy,dEnf,dclasse,dclasseEnf,fichier,path):
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

