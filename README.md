# LesBellesGosses
Projet Barstar de BEI WANYING et BERTRAND MATHILDE

## Introduction : 
On cherche à étudier la stabilité d'une structure protéique et les changements conformationnels qu'elle subit au cours du temps.
2 analyses ont été réalisées : 
	- Une analyse globale : identifier des mouvements de boucles de grandes amplitudes, des mouvements allostériques, du dépliement ou encore de la compaction
	- Une analyse locale : identifier des petits mouvements de boucles ou encore des mouvements de chaîne latérales.

Le dossier d'étude doit contenir le fichier de la structure protéique de référence nommé "start_prot_only.pdb" ainsi que le fichier de la protéine avec ses conformations

## Pré-requis : 
Version python utilisée : python 2.7
#### Packages: 
-numpy
-matplotlib

## Exécution du code : 

usage: Main.py [-h] [-p P] [-a A]

optional arguments:
  -h, --help  show this help message and exit
  -p P        Path to directory where pdb files are stored
  -a A        Type of analysis: global or local or both

Exemple : python Main.py -p /Users/mathildebertrand/Desktop/M1-S2/python/gitt/LesBellesGosses -a global


