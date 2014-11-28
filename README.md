BamSeqChksumJava
================

This project contains a prototype implementation of BAM checksumming in Java.

# Compilation

The make file can be used to create a JAR file from the Java source files. This requires a JDK.

# Running

The program can be run using

```
java -jar BamSeqChksum.jar < in.bam
```

It reads the input BAM file from the standard input channel and produces the checksum information on standard output.

The implementation supports three checksums types:

* crc32prod: checksums for single alignments are produced by crc32 and combined over multiple records by multiplication in
  a prime number field. The prime used is 2^31-1.
* sha512: checksums for single alignments are produced by SHA-512. They are combined over multiple records by adding up and truncating to 512 bit numbers.
* sha512primesums512: checksums for single alignments are produced by SHA-512. They are combined over multiple records by adding up in a prime number field.
  The prime used is 2^512+75.

  