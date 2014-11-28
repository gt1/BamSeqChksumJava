#! /usr/bin/make
SOURCES=\
	BamCRC32ProdChecksumRecordHandler.java \
	BamHeaderParser.java \
	BamParser.java \
	BamRecordHandler.java \
	BamSeqChksum.java \
	BamSHA512ChecksumRecordHandler.java \
	BamSHA512PrimeSumsChecksumRecordHandler.java \
	BgzfInput.java \
	UnsignedNumber.java
OBJECTS=BamCRC32ProdChecksumRecordHandler*.class \
	BamHeaderParser*.class \
	BamParser*.class \
	BamRecordHandler*.class \
	BamSeqChksum*.class \
	BamSHA512ChecksumRecordHandler*.class \
	BamSHA512PrimeSumsChecksumRecordHandler*.class \
	BgzfInput*.class \
	UnsignedNumber*.class

all: BamSeqChksum.jar

BamSeqChksum.jar: ${SOURCES}
	javac ${SOURCES}
	echo "Main-Class: BamSeqChksum" > manifest.txt
	jar cmf  manifest.txt $@ ${OBJECTS}
	rm -f manifest.txt

clean:
	rm -f ${OBJECTS}
distclean: clean
	rm -f BamSeqChksum.jar
