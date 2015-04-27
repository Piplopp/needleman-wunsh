@author : Jerome PIVERT p1308594
@resume : Le programme permet d'Aligner 2 séquences de type nucléotidique ou protéique
	  de façon globale -algorithme de Needleman- ou de façon locale -algorithme de Smith & Waterman-
	  les extension de Gap sont prises en compte mais pas les ex aequo

	  Une option help dans le menu du programme permet d'afficher ce README.txt

------------------
Comment lancer : -
------------------
./alignement_main.py [-h] [-p] [infile infile]
./alignement_main.py -h pour afficher l'aide

ex: ./alignement_main.py --prot seq_prot1.txt seq_prot2.txt		--> Séquences Protéiques par fichiers
ex: ./alignement_main.py						--> Séquences Nucléiques à la main

Les fichiers sont optionnels mais doivent être donnés par 2, une structure de contrôle prévient si le nombre
d'argument est erroné.


------------------
Modifier valeurs:-
------------------
Pour modifier les valeurs de match, mismatch, gap, extension de gap:

Modifier le fichier alignements.py
	lignes 12 à 15.



------------------
Les fichiers :   -
------------------
@file alignement.py	: fichier contenant les differentes fonctions permettant les alignements
			  contient également les valeurs des gap/matchs etc
@file blossum62.py 	: fichier contenant uniquement la fonction pour utiliser la matrice blossum62
@file alignement_main.py: fichier d'execution contenant les verifications premieres ainsi que le menu principal


------------------
Les modules :    -
------------------
@import sys, os: gestion des commandes système
@import string
@import argparse: gestion des options de la ligne de commande
@import numpy: gestion des matrices -array-
@import re: gestion des regexp

@import alignement
@import blossum62


------------------
Les fonctions :  -
------------------
alignement.py
	@fct init_matrice(lignes, colonnes, data_type)
		@resume: créé une matrice vide/nulle et la retourne
		@param lignes: <type int> nombre de ligne de la matrice
		@param colonnes: <type int> nombre de colonnes de la matrice
		@param data_type: <type str> type de données contenues dans la matrice
		@return: matrice de data_type vide/nulle de taille lignes x colonnes

	@fct check(a, b)
		@resume: regarde si a et b sont identiques ou non et retourne la valeur de match ou mismatch
		@param a,b: <type char> nucléotide en cours de traitement des séquences à aligner
		@return: valeur de match ou mismatch selon les conditions
	@fct needleman(seq1, seq2, boolean)
		@resume: effectue l'alignement de Needleman -global-
			 En premier lieu elle appelle la fonction "init_matrice" afin de créér les
			 matrices score et traceback, on les initialise ensuite suivant l'algorithme de
			 Needleman puis on remplit en parallèle les matrices "score" et "traceback"

			 Si on ne gère pas l'extension de gap on peut aisément se passer d'une des 2
			 matrices. Mais cette extension demande de constamment vérifier si notre case
			 précédente correspond à un GAP ou non, d'où la necessité de remplir les 2 matrices

			 On fait ensuite le traceback en lisant les données de la matrice depuis la fin
			 et en concaténant dans 2 variables les résultats de l'alignement.

			 On print finalement le reverse des 2 séquences dans un fichier
			 *Dans le cas de séquences protéiques (boolean=True) on utilise @fct blossum62()
			  plutôt que @fct check()
		@param seq1,2: <type str> sequences nucléotidiques ou protéiques
		@param boolean: <type bool> si oui ou non les sequences sont proteiques
		@out: fichier alignement_needleman.txt comprenant les 2 séquences alignées
	@fct waterman(seq1, seq2, boolean)
		@resume: effectue l'alignement de Smith & Waterman -local-
			 En premier lieu elle appelle la fonction "init_matrice" afin de créér les
			 matrices score et traceback, ici inutile de les initialiser puisqu'elles sont créées
			 vides. On remplit en parallèle les matrices "score" et "traceback" toujours à cause
			 de la gestion des extensions de GAP.

			 La principale différence ici viens du fait que l'on préfère ne pas aligner plutôt que mal
			 aligner (d'ou max(0...))

			 On fait ensuite le traceback en lisant les données de la matrice depuis la première
			 occurence de la plus grande valeur et enn concaténant dans 2 variables les résultats de
			 l'alignement.

			 On print finalement le reverse des 2 séquences dans un fichier
			 *Dans le cas de séquences protéiques (boolean=True) on utilise @fct blossum62()
			  plutôt que @fct check()
		@param seq1,2: <type str> sequences nucléotidiques ou protéiques
		@param boolean: <type bool> si oui ou non les sequences sont proteiques
		@out: fichier alignement_waterman.txt comprenant les 2 séquences alignées

blossum62.py
	@fct blossum62(aa1, aa2)
		@resume: contient un dictionnaire de correspondance d'AA selon Blossum62
		@param aa1,2 : <type char> aa en cours de traitement dans l'algorithme d'alignement
		@return: valeur de pénalité entre les 2 aa données en paramètres

alignement_main.py
	@fct main()
		@resume: permet de lancer le programme, elle contient les commandes de base ainsi qu'un menu
			 C'est cette fonction qui va permettre à l'utilisateur de choisir de rentrer ses séquences
			 d'utiliser des fichiers, d'indiquer les protéines ainsi qu'offrir le choix des 2
			 algorithmes.
					
			 Une boucle permet l'affichage du menu tant que l'utilisateur ne rentre pas de bonne options
			 -1,2,q,help-
			 
			 La structure try/except permet de gérer les options string/int du menu en effet les 
			 conditions if/elif/else utilisent tous la variable associée au choix utilisateur, 
			 celui-ci pouvant entrer ce que bon lui semble.
	@fct verif_seq(seq1,seq2)
		@resume: verifie si les sequences sont justes et ne contiennent pas de caractères
			 non reconnus par le programme.
			 
			 Si l'option -p est utilisée alors on vérifie que la séquence ne contient que 
			 les caractères utilisés par Blossum62
			 Sinon on vérifie la présence de ACGT
		@param seq1,2: <type str> les 2 séquences
		@return: <type bool> return True ou False si les séquences sont justes ou non
	@fct open_file(file1, file2)
		@resume: ouvre et entre en mémoire les fichiers passés en arguments.
			 Appelle ensuite la fonction @fct verif_seq() afin de vérifier les séquences.
			 Si les séquences sont erronées le programme se ferme.
		@param file1, file2: <type str> noms des fichiers donnés en arguments, chacun doit contenir
				     une séquence uniligne.
		@return: les 2 séquences des fichiers sous forme de string
	@fct user_inpupt_seq()
		@resume: fonction utilisée pour que l'utilisateur entre ses propres séquences.
			 Celles-ci sont vérifiées par @fct verif_seq() et une boucle permet de recommencer
			 si les séquences sont erronnées.
		@return: renvoie les 2 séquences entrées pour l'utilisateur



q pour quitter (README.TXT open with less)
