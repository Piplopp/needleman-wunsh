#!/usr/bin/python

import sys, os
import string
import alignements
import argparse
import re


parser=argparse.ArgumentParser(usage='./%(prog)s [-h] [-p] [infile infile]\n', description='Alignement de sequences selon "Needleman" ou "Smith & Waterman"')

parser.add_argument('-p', '--prot',
	help = 'specifie que vous traitez des sequences proteiques',
	action="store_true"
)
parser.add_argument('infile',
	nargs='*',
	help='specifier 0 ou 2 fichiers sequences',
	default=[False, False],
	type=argparse.FileType('r')
)

args=parser.parse_args()


if len(args.infile) != 2: # Si 1 fichier
	parser.error('2 fichiers attendus')
if len(sys.argv) > 1:
	if sys.argv[1]=='-p' or sys.argv[1]=='--p' or sys.argv[1]=='--prot': # pour conserver les fichiers en position 1 et 2
		sys.argv.remove(sys.argv[1])# dans le sys.argv (utile par la suite)





# Verification sequence
def verif_seq(seq1, seq2):
	if args.prot==True and re.match('^[ARNDCQEGHILKMFPSTWYVBZX]+$', seq1.upper()) and re.match('^[ARNDCQEGHILKMFPSTWYVBZX]+$', seq2.upper()):
		return True  #seq proteique OK
	elif args.prot==True:
		return False #seq proteique non OK
	elif args.prot == False and re.match('^[atgc]+$', seq1.lower()) and re.match('^[atgc]+$', seq2.lower()): 
		return True  #seq nucleique OK
	else: 
		return False #seq nucleique non OK


# Ouverture des fichiers
def open_file(file1,file2):
	seq1, seq2 = '', ''

	f1 = open(file1, 'r')
	seq1 = f1.readline()
	seq1 = string.strip(seq1).lower()
	f1.close()

	f2 = open(file2, 'r')
	seq2 = f2.readline()
	seq2 = string.strip(seq2).lower()
	f2.close()
	
	if verif_seq(seq1, seq2):
		return seq1, seq2
	else:
		sys.exit('erreur dans la sequence')


# Entree manuelle des sequences
def user_input_seq():
	seq1, seq2 = '', ''
	while verif_seq(seq1, seq2) == False:
		seq1=string.strip(raw_input('Sequence 1 : ')).lower()
		seq2=string.strip(raw_input('Sequence 2 : ')).lower()

		print 'Si les sequences contiennent de mauvais caracteres elles vous seront redemandees'
		
		raw_input()
		#os.system('clear')

	return seq1,seq2





def main():
	os.system('clear')
	seq1, seq2 = '', ''

	if args.infile[0] != False and args.infile[1] != False:
		seq1,seq2 = open_file(sys.argv[1], sys.argv[2])
	else:
		seq1,seq2 = user_input_seq()

	i = -9


	# Menu
	while i == -9:
		os.system('clear')
		if args.infile[0] != False and args.infile[1] != False:
			print '<'+sys.argv[1]+' read><'+sys.argv[2]+' read>'
		else:
			print 'seq1: '+seq1
			print 'seq2: '+seq2

		print 'proteine: '+str(args.prot)
		print "\n\t1 - Alignement local - Smith & Waterman"
		print "\t2 - Alignement global - Needleman"
		print "\tq - Exit\n"
		a = raw_input('Choix de l\'alignement : ')
		a = a.lower()
		
		try:
			if a == 'quit' or a == 'q' or a == 'exit':
				sys.exit()
			elif a == 'help':
				os.system('less README.txt')
			elif int(a) == 1:
				alignements.waterman(seq1, seq2, args.prot)			
				raw_input('Smith & Waterman - Done - out:alignement_waterman.txt - Appuyez sur une <Enter>')
				i = -9
			elif int(a) == 2:
				alignements.needleman(seq1, seq2, args.prot)
				raw_input('Needleman - Done - out:alignement_needleman.txt - Appuyez sur une <Enter>')
				i = -9
			else:
				raw_input('erreur mauvais param')
				i = -9
		except ValueError: #Gestion des entree [A-Za-z]*[0-9]*
				raw_input('erreur mauvais param')

main()
