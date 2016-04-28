all:
	mkdir -p bin
	ln -sf ../Segmentum.py bin/Segmentum
	gcc -O3 -std=c99 -o compiled/spileup compiled/spileup.c