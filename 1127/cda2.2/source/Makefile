#
# Makefile f"ur die endg"ultige Version der CDA
# 
# Stefan Dapprich, 25.07.1995
# 
# Corrections for GAUSSIAN03 by Andreas Krapp, 06.03.2004
# Corrections for simultaneous GAUSSIAN03 and GAUSSIAN09 compatibility by Moritz von Hopffgarten, 13.11.2011
#

CC      = cc
OPTIONS = -O -Wno-unused-result

cda2_2: main.o readgaus.o matrix.o output.o
	$(CC) $(OPTIONS) main.o readgaus.o matrix.o output.o -o cda2_2

main.o: main.c
	$(CC) $(OPTIONS) -c main.c

readgaus.o: readgaus.c
	$(CC) $(OPTIONS) -c readgaus.c

matrix.o: matrix.c
	$(CC) $(OPTIONS) -c matrix.c

output.o: output.c
	$(CC) $(OPTIONS) -c output.c

