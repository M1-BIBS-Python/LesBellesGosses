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
                dico_models[dico_molecule][chain] = {}
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
	print temps
	infile.close()
	return(temps)
	
#Calcule de la distance entre 2 residus selon 2 methodes de calcule:
#1. Distance la plus courte entre les 2 residus
#2. Distance entre les centres de masse

def Distance_res_min(dico1,dico2):   #Mode de calcul 1 : distance la plus courte
#Dico 1 et Dico2 sont des dico de residus
    min=200 #On fixe une distance min
    
    for atom1 in dico1['atomlist']:#Calcul de la distance entre les 2 residus
        coor1=[dico1[atom1]["x"],dico1[atom1]["y"],dico1[atom1]["z"]]
        for atom2 in dico2["atomlist"] :
            coor2=[dico2[atom2]["x"],dico2[atom2]["y"],dico2[atom2]["z"]]
            distance=Distance(float(coor1[0]),float(coor1[1]),float(coor1[2]),float(coor2[0]),float(coor2[1]),float(coor2[2]))
            
            #print(distance)
            if distance>0:   #A mettre ? car si=0, les residus sont a la meme position...
                if distance<min: #On en deduit la distance min
                    min=distance
    print(min)
    return(min)


def CM(listx,listy,listz):
    x=sum(listx)/float(len(listx))
    y=sum(listy)/float(len(listy))
    z=sum(listz)/float(len(listy))
    coords=[x,y,z]
    return coords
   
    
def Distance(x1,y1,z1,x2,y2,z2): #Calcule la distance entre deux points
    return(sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2))


def RMSD(list_delta): #calcul de RMSD
    distcarre=[]
    for delta in list_delta:
        distcarre.append(delta**2)
    return (sqrt((sum(distcarre))/float(len(list_delta))))


def CMglob(dico):#c'est le dico[key] qu'on passe ici
    globx=[] #list permet de stocker le x de tous les atomes d'une prot
    globy=[]
    globz=[]
    glob={}
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
            glob["%s-%s"%(res,dico[chain][res]["resname"])]=CM(listx,listy,listz)
    glob["prot"]=CM(globx,globy,globz)
    return glob
    
    
def giration(dico):#c'est le dico[key] qu'on passe ici
    dico_CM=CMglob(dico)
    list_dist=[]
    for key in dico_CM.keys():
        list_dist.append(Distance(dico_CM["prot"][0],dico_CM["prot"][1],dico_CM["prot"][2],dico_CM[key][0],dico_CM[key][1],dico_CM[key][2]))
    return max(list_dist)  

    
def writefile_glob():
    if os.path.isdir(path): #si le chemin est vers un dossier
        out=open("%s/PythonProgResults/output_analyse_global_%s"%(path,fichier),"w")
    else: #si le chemin est vers un fichier
        out=open("%s/PythonProgResults/output_analyse_global_%s"%(os.path.dirname(path),os.path.basename(fichier)),"w")
    
    #est-ce utile ?  On veut juste la conformation, la giration et le RMSD ?
    #ecrire la structure d'origine (il s'agit des lignes pour le modele 0)
    #~ with open(fichier, "r") as filin:
        #~ line=filin.readline()
        #~ while line != "ENDMDL\n":
            #~ line=filin.next()
            #~ out.write(line)
    #~ filin.close()
    
    #ecrire dans l'ordre les resultats de la comparaison des conformations avec la structure d'origine
    out.write("\nConformation \t RMSD results \t\t Giration\n")
    for i in range(len(dico)):
        out.write("\n%s: \t\t %.12f \t %.12f"%(i,dico_RMSD["%s"%i],dico_Giration["%s"%i]))
        
    out.close()
    

def writefile_local():
    if os.path.isdir(path): #si le chemin est vers un dossier
        out=open("%s/PythonProgResults/output_analyse_local_%s"%(path,fichier),"w")
    else: #si le chemin est vers un fichier
        out=open("%s/PythonProgResults/output_analyse_local_%s"%(os.path.dirname(path),os.path.basename(fichier)),"w")
        
    for i in range(len(dico_dist)):
		out.write("MODEL %s\n"%i)
		out.write("Residu \t\t Distance\n")
		for key in dico_dist["%s"%i].keys():
			out.write("%s \t\t %.12f\n"%(key,dico_dist["%s"%i][key]))
        
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

    ########################################################
    #               Etapes preliminaires
    ########################################################

    #on supprime le dossier contenant les resultats de l'execution precedente (s'il existe) afin d'eviter les chevauchements
    #le dossier "PythonProgResults" sert a stocker les fichiers sortis
    
    if os.path.isdir(path): # si le chemin est vers un dossier
        shutil.rmtree("%s/PythonProgResults"%path) 
        os.mkdir("%s/PythonProgResults"%path)
    else: #si le chemin est vers un fichier
        try:
            shutil.rmtree("%s/PythonProgResults"%os.path.dirname(path)) 
            os.mkdir("%s/PythonProgResults"%os.path.dirname(path))
        except: #s'il s'agit de la premiere execution
            os.mkdir("%s/PythonProgResults"%os.path.dirname(path)) 
    
    try:
        os.chdir(path)
        fichiers=glob.iglob("*.pdb") # si le path est vers un dossier, on lira tous les fichiers pdb
    except OSError:
        fichiers=glob.iglob("%s"%path) # si le path est vers un fichier, on lira ce ficher

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
                        #pk comparer les distances que entre les Ca ?
							
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
            plt.subplot(2,1,1) #Partage la fenetre pour les emplacements des graphs
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
            
            #dist de chaque residu de chaque conformation par rapport au CM de la prot
            dico_dist={}
            for key in dico: #pour chaque conformation
                dico_Enfouissement=CMglob(dico[key])
                dico_dist[key]={} #dico pour une conformation, dist de chaque residu de chaque conformation seront stockee dedans
                for cle in dico_Enfouissement:
					if cle != "prot":
						dico_dist[key][cle]=Distance(dico_Enfouissement["prot"][0],dico_Enfouissement["prot"][1],dico_Enfouissement["prot"][2],dico_Enfouissement[cle][0],dico_Enfouissement[cle][1],dico_Enfouissement[cle][2])

            writefile_local()
                
            #RMSD moyen pour chaque residu de chaque conformations (a ajouter)
                
    #~ else : #si vous voulez a la fois une analyse globale et une analyse locale
        #~ for fichier in fichiers:
            #~ #do sth
            
    
    


     
