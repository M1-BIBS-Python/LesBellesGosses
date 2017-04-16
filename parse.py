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
            flex_res["%s"%res]=RMSD(delta)
    return flex_res #dico ayant ResSeq comme cle, et RMSD local comme valeur
                    
    
def giration(dico):#c'est le dico[key] qu'on passe ici
    dico_CM=CMglob(dico)
    list_dist=[]
    for res in dico_CM.keys():
        if res != "residulist":
            list_dist.append(Distance(dico_CM["prot"][0],dico_CM["prot"][1],dico_CM["prot"][2],dico_CM[res][0],dico_CM[res][1],dico_CM[res][2]))
    return max(list_dist)  

    
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
            
            list_RMSD=[]
            list_conformation=[]
            list_Giration=[]
            
  
            for key in dico:
                list_delta=[]
                dico_Giration[key]=giration(dico[key])
                list_Giration.append(giration(dico[key]))
                list_conformation.append(float(key))               #On recupere le numero de la conformation
                for chain in dico[key]["chains"]:
                     
                     for res in dico[key][chain]["reslist"]:
                         list_delta.append(Distance(float(dico[key][chain][res]['CA']['x']),float(dico[key][chain][res]['CA']['y']),float(dico[key][chain][res]['CA']['z']),float(dico['0'][chain][res]['CA']['x']),float(dico['0'][chain][res]['CA']['y']),float(dico['0'][chain][res]['CA']['z']))) # compare dist entre les Ca         
            
                            
                dico_RMSD[key]=RMSD(list_delta)
                list_RMSD.append(RMSD(list_delta))  
                
            writefile_glob()     
            
            
            #Je propose pour interpreter les resultats (pour loral):
            #une analyse visuelle des resultats : graph representant RMSD/Giration en fonction de la conformation
            
            y=np.array(list_RMSD) #RMSD 
            y2=np.array(list_Giration)
            x=np.array(list_conformation) #La conformation=numero du modele
            plt.scatter(x,y,c='red')
            plt.scatter(x,y2,c='blue')
            axes = plt.gca()
            axes.set_xlim(-30, 2100)
            axes.set_ylim(-1,25)
            plt.title('Variations RMSD (ou rayons de Giration) en fonction de la conformation')
            plt.legend(['RMSD','Giration'])
            plt.show()
            
            #On peut aussi regarder la variation de RMSD en fonction du temps : 
            
            x=np.array(list_temps)
            plt.plot(x,y)
            plt.title('Evolution du RMSD en fonction du temps')
            plt.show()
            
            #Idem pour la variation de Giration en fonction du temps: 
            x=np.array(list_temps)
            plt.plot(x,y2)
            plt.title('Evolution du rayon en fonction du temps')
            plt.show()
            #Il me reste les legendes a mettre sur les boucles
            
    elif (analyse == "local"): #si vous voulez seulement une analyse locale
        for fichier in fichiers:
            print "Parsing:",fichier
            dico=ParsingPDB(fichier)
            
            #################################
            #Donnees pour les graphs
            list_dist=[]
            list_res=[] #Contient le numero de chaque residus
            list_temps=Temps(fichier) 
            ############################################
            
            #dist de chaque residu de chaque conformation par rapport au CM de la prot
            #et calcul de RMSD local
            dico_dist={}
            dico_RMSD={} #Va contenir les RMSD de chaque residus par rapport a la position dans la structure
            for key in dico: #pour chaque conformation
                dico_Enfouissement=CMglob(dico[key])
                dico_RMSD[key]=RMSDlocal(dico["0"],dico[key]) #comparaison de deux conformations (dans la fct, comparaison de chaque residu des deux conformations)
                    
                dico_dist[key]={} #dico pour une conformation, dist de chaque residu de chaque conformation seront stockee dedans
                
                for res in dico_Enfouissement:
                    
                    if res != "prot" and res != "residulist":
                        dico_dist[key][res]=Distance(dico_Enfouissement["prot"][0],dico_Enfouissement["prot"][1],dico_Enfouissement["prot"][2],dico_Enfouissement[res][0],dico_Enfouissement[res][1],dico_Enfouissement[res][2])
                        list_dist.append(Distance(dico_Enfouissement["prot"][0],dico_Enfouissement["prot"][1],dico_Enfouissement["prot"][2],dico_Enfouissement[res][0],dico_Enfouissement[res][1],dico_Enfouissement[res][2]))
                    if res == "residulist":
                        dico_dist[key]["residulist"]=dico_Enfouissement["residulist"]
            
                        #On attribue des numeros a chacun des AA:
                        #~ if cle == 'CYS':
                            #~ list_res.append(1)
                        #~ if cle == 'GLU':
                            #~ list_res.append(2)
                        #~ if cle == 'VAL':
                            #~ list_res.append(3)
                        #~ if cle == 'PRO':
                            #~ list_res.append(4)
                        #~ if cle == 'ALA':
                            #~ list_res.append(5)
                        #~ if cle == 'THR':
                            #~ list_res.append(6)
                        #~ if cle == 'PHE':
                            #~ list_res.append(7)
                        #~ if cle == 'ILE':
                            #~ list_res.append(8)
                        #~ if cle == 'LYS':
                            #~ list_res.append(9)
                        #~ if cle == 'GLY':
                            #~ list_res.append(10)
                        #~ if cle == 'ASP':
                            #~ list_res.append(11)
                        #~ if cle == 'LEU':
                            #~ list_res.append(12)
                        #~ if cle == 'HIS':
                            #~ list_res.append(13)
                        #~ if cle == 'SER':
                            #~ list_res.append(14)
                        #~ if cle == 'TRP':
                            #~ list_res.append(15)
                        #~ if cle == 'GLN':
                            #~ list_res.append(16)
                        #~ if cle == 'TYR':
                            #~ list_res.append(17)
                        #~ if cle == 'ASN':
                            #~ list_res.append(18)
                        #~ if cle == 'ARG':
                            #~ list_res.append(19)
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
            
    
    


     
