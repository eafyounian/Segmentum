all:
	mkdir -p bin
	ln -sf ../Segmentum.py bin/Segmentum
	chmod u+x ./Segmentum.py
	gcc -O3 -std=c99 -o compiled/spileup compiled/spileup.c