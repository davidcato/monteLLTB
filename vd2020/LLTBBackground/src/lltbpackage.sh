#!/bin/sh


	mkdir -p tmppackage
	rm -rf tmppackage/*
	cp ../lib/*a tmppackage/.
	cd tmppackage 
    rm -f liblltbpack.a
	for y in `ls *.a`; do
	    ar x ${y}
	done
	ar -csrv liblltbpack.a *.o
	mv liblltbpack.a ../../lib/.
	cd ..
	rm -rf tmppackage
