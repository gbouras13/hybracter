CC ?=				gcc
CFLAGS ?=			-g -O3
SAMTOOLS_PATH ?=	samtools-0.1.19

INSTALL ?=	install -c
MKDIR ?=	mkdir -p

DESTDIR ?=	../stage
PREFIX ?=	/usr/local
BINDIR :=	$(PREFIX)/bin

UNAME := $(shell uname)

ifeq ($(UNAME), Linux)
GNU_FLAG ?= -D_GNU_SRC
else ifeq ($(UNAME), Solaris)
else ifeq ($(UNAME), Darwin)
endif

DEFAULT: all

.PHONY: all install clean cleantest test example cleanexample

samlib:
	make -C $(SAMTOOLS_PATH) lib

$(SAMTOOLS_PATH)/bam.o: samlib
	true

samtools_helper.o: samlib
	$(CC) -c $(CFLAGS) $(GNU_FLAG) samtools_helper.c -I$(SAMTOOLS_PATH) 

ALElike.o: ALElike.c ALElike.h
	$(CC) -c $(CFLAGS) $(GNU_FLAG) ALElike.c -I$(SAMTOOLS_PATH)

ALEhelpers.o: ALEhelpers.c ALEhelpers.h
	$(CC) -c $(CFLAGS) $(GNU_FLAG) ALEhelpers.c -I$(SAMTOOLS_PATH)

ALE: ALE.c ALE.h ALEhelpers.o ALElike.o samtools_helper.o $(SAMTOOLS_PATH)/bam.o ssw.o
	$(CC) $(CFLAGS) $(GNU_FLAG) ALE.c ALElike.o ALEhelpers.o ssw.o -o ALE -lz -lm -I$(SAMTOOLS_PATH) -L$(SAMTOOLS_PATH) -lbam samtools_helper.o -lpthread

synthReadGen: synthReadGen.c
	$(CC) $(CFLAGS) synthReadGen.c -o synthReadGen -lz -lm -I$(SAMTOOLS_PATH) -L$(SAMTOOLS_PATH)

readFileSplitter: readFileSplitter.c
	$(CC) $(CFLAGS) readFileSplitter.c -o readFileSplitter

GCcompFinder: GCcompFinder.c
	$(CC) $(CFLAGS) GCcompFinder.c -o GCcompFinder

all: GCcompFinder readFileSplitter synthReadGen ALE

install: all
	$(MKDIR) $(DESTDIR)$(BINDIR)
	$(INSTALL) GCcompFinder readFileSplitter synthReadGen ALE $(DESTDIR)$(BINDIR)

example: ALE synthReadGen
	./image_maker.py

clean: cleantest cleanexample
	rm -f *.o GCcompFinder readFileSplitter synthReadGen2 ALE 
	make -C $(SAMTOOLS_PATH) clean

cleantest: 
	rm -f ../test/*-testPhiX.ale* ../test/*-tinytestPhiX.ale* ../test/*-tinytestPhiX-broken.ale* ../test/*-testRealignPhiX.ale* ../test/*-testSNPPhiX.ale*
		
cleanexample:
	rm -rf ../example/Ecoli_first350k_* ../example/*.fna* ../example/*.fastq

test: ALE cleantest test1 test2 test3 test4 test5 test-single-relations test-pair-relations test-pair-chimer-baseline test-pair-chimer

%-testPhiX.ale: % ALE
	./ALE $< ../test/PhiX174.fsa $@ > $@.log 2>&1 && ls -l $@

%-tinytestPhiX.ale: % ALE ../test/PhiX-tiny.sam-testPhiX.ale
	./ALE --pm ../test/PhiX-tiny.sam-testPhiX.ale.param $< ../test/PhiX174.fsa $@ > $@.log 2>&1 && ls -l $@

%-tinytestPhiX-broken.ale: % ALE ../test/PhiX-tiny.sam-testPhiX.ale
	./ALE --pm ../test/PhiX-tiny.sam-testPhiX.ale.param $< ../test/PhiX174-broken.fsa $@ > $@.log 2>&1 && ls -l $@

%-testRealignPhiX.ale: % ALE
	./ALE --realign $< ../test/PhiX174.fsa $@ > $@.log 2>&1 && ls -l $@

%-testSNPPhiX.ale: % ALE
	./ALE --SNPreport $@.SNPreport.txt $< ../test/PhiX174-SNP.fsa $@ > $@.log 2>&1 && ls -l $@

test1: ../test/PhiX-sort.bam-testPhiX.ale
	diff -q $< ../test/PhiX-sort.ale > $<.diff && echo "$@ PASSED"

test2: ../test/PhiX.bam-testPhiX.ale
	diff -q $< ../test/PhiX.ale > $<.diff && echo "$@ PASSED"

test3: ../test/PhiX.sam-testPhiX.ale
	diff -q $< ../test/PhiX.ale > $<.diff && echo "$@ PASSED"

test4: ../test/PhiX-tiny.sam-testPhiX.ale
	diff -q $< ../test/PhiX-tiny.ale > $<.diff && echo "$@ PASSED"

test5: ../test/PhiX.sam-testRealignPhiX.ale
	diff -q $< ../test/PhiX.ale > $<.diff && echo "$@ PASSED"

test6: ../test/PhiX.sam-testSNPPhiX.ale
	diff -q $< ../test/PhiX-SNP.ale > $<.diff \
	&& diff -q $<.SNPreport.txt ../test/PhiX-SNP.ale.SNPreport.txt \
	&& echo "$@ PASSED"

test-single-relations: ../test/PhiX-tiny.sam-testPhiX.ale  ../test/PhiX-tiny-partial-single.sam-testPhiX.ale  ../test/PhiX-tiny-partial2-single.sam-testPhiX.ale  ../test/PhiX-tiny-partial3-single.sam-testPhiX.ale  ../test/PhiX-tiny-partial4-single.sam-testPhiX.ale  ../test/PhiX-tiny-partial5-single.sam-testPhiX.ale  ../test/PhiX-tiny-+unmapped-single.sam-testPhiX.ale  ../test/PhiX-tiny-unmapped-single.sam-testPhiX.ale
	awk '/ALE_score/ {print}' $^ > /tmp/test-relations.alescore && sort -k3 -n -r /tmp/test-relations.alescore | diff - /tmp/test-relations.alescore && rm -f /tmp/test-relations.alescore && echo "$@ PASSED"
		
test-pair-relations: ../test/PhiX-tiny-stretch0-pair.sam-tinytestPhiX.ale ../test/PhiX-tiny.sam-testPhiX.ale ../test/PhiX-tiny-stretch1-pair.sam-tinytestPhiX.ale ../test/PhiX-tiny-stretch2-pair.sam-tinytestPhiX.ale ../test/PhiX-tiny-stretch3-pair.sam-tinytestPhiX.ale ../test/PhiX-tiny-stretch4-pair.sam-tinytestPhiX.ale ../test/PhiX-tiny-stretch5-pair.sam-tinytestPhiX.ale ../test/PhiX-tiny-stretch6-pair.sam-tinytestPhiX.ale ../test/PhiX-tiny-stretch7-pair.sam-tinytestPhiX.ale ../test/PhiX-tiny-stretch8-pair.sam-tinytestPhiX.ale ../test/PhiX-tiny-removed2-pair.sam-tinytestPhiX.ale
	awk '/ALE_score/ {print}' $^ > /tmp/test-relations.alescore && sort -k3 -n -r /tmp/test-relations.alescore | diff - /tmp/test-relations.alescore && rm -f /tmp/test-relations.alescore && echo "$@ PASSED"

test-pair-relations2: ../test/PhiX-tiny-stretch0-pair.sam-tinytestPhiX.ale ../test/PhiX-tiny-removed-pair.sam-tinytestPhiX.ale ../test/PhiX-tiny-removed2-pair.sam-tinytestPhiX.ale
	awk '/ALE_score/ {print}' $^ > /tmp/test-relations.alescore && sort -k3 -n -r /tmp/test-relations.alescore | diff - /tmp/test-relations.alescore && rm -f /tmp/test-relations.alescore && echo "$@ PASSED"

test-pair-chimer-baseline: ../test/PhiXbroken-tiny-stretch0-pair.sam-tinytestPhiX-broken.ale ../test/PhiXbroken-tiny.sam-tinytestPhiX-broken.ale ../test/PhiXbroken-tiny-stretch1-pair.sam-tinytestPhiX-broken.ale ../test/PhiXbroken-tiny-stretch2-pair.sam-tinytestPhiX-broken.ale ../test/PhiXbroken-tiny-stretch3-pair.sam-tinytestPhiX-broken.ale ../test/PhiXbroken-tiny-stretch4-pair.sam-tinytestPhiX-broken.ale ../test/PhiXbroken-tiny-stretch5-pair.sam-tinytestPhiX-broken.ale ../test/PhiXbroken-tiny-stretch6-pair.sam-tinytestPhiX-broken.ale ../test/PhiXbroken-tiny-stretch7-pair.sam-tinytestPhiX-broken.ale ../test/PhiXbroken-tiny-stretch8-pair.sam-tinytestPhiX-broken.ale
	awk '/ALE_score/ {print}' $^ > /tmp/test-relations.alescore && sort -k3 -n -r /tmp/test-relations.alescore | diff - /tmp/test-relations.alescore && rm -f /tmp/test-relations.alescore && echo "$@ PASSED"

test-pair-chimer: ../test/PhiXbroken-tiny-stretch0-pair.sam-tinytestPhiX-broken.ale ../test/PhiXbroken-tiny.sam-tinytestPhiX-broken.ale ../test/PhiXbroken-tiny-stretch1-pair.sam-tinytestPhiX-broken.ale ../test/PhiXbroken-tiny-stretch2-pair.sam-tinytestPhiX-broken.ale ../test/PhiXbroken-tiny-stretch3-pair.sam-tinytestPhiX-broken.ale ../test/PhiXbroken-tiny-stretch4-pair.sam-tinytestPhiX-broken.ale ../test/PhiXbroken-tiny-stretch5-pair.sam-tinytestPhiX-broken.ale ../test/PhiXbroken-tiny-stretch6-pair.sam-tinytestPhiX-broken.ale ../test/PhiXbroken-tiny-chimer-pair.sam-tinytestPhiX-broken.ale # ../test/PhiXbroken-tiny-removed-pair.sam-tinytestPhiX-broken.ale ../test/PhiXbroken-tiny-removed2-pair.sam-tinytestPhiX-broken.ale
	awk '/ALE_score/ {print}' $^ > /tmp/test-relations.alescore && sort -k3 -n -r /tmp/test-relations.alescore | diff - /tmp/test-relations.alescore && rm -f /tmp/test-relations.alescore && echo "$@ PASSED"
