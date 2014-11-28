#! /bin/bash
PATH=/software/hpag/biobambam/latest/bin:$PATH

function finish
{
	rm -f log_biobambam
	rm -f log_java
	rm -f time_biobambam
	rm -f time_java
}

trap finish EXIT SIGINT SIGTERM

for i in $* ; do
	for hash in crc32prod sha512 sha512primesums512 ; do
		LANG=C /usr/bin/time bamseqchksum hash=${hash} < $i 2>&1 >log_biobambam | awk '/elapsed/ {V=$3 ; sub(/elapsed/,"",V) ; print V }' >time_biobambam &
		BIOBAMBAMPID=$!
		LANG=C /usr/bin/time java -Xmx64m -jar BamSeqChksum.jar  hash=${hash} < $i 2>&1 >log_java | awk '/elapsed/ {V=$3 ; sub(/elapsed/,"",V) ; print V }' >time_java &
		JAVAPID=$!
		
		wait ${BIOBAMBAMPID}
		RETBIOBAMBAM=$?
		
		wait ${JAVAPID}
		RETJAVA=$?
		
		if test ${RETBIOBAMBAM} -ne 0 ; then
			echo "biobambam bamseqchksum failed"
		fi
		if test ${RETJAVA} -ne 0 ; then
			echo "Java BamSeqChksum failed"
		fi
		
		BIOBAMBAMTIME=`cat time_biobambam`
		JAVATIME=`cat time_java`
		printf "$i\t$hash\t${BIOBAMBAMTIME}\t${JAVATIME}\n"
		diff -c log_biobambam log_java
	done
done
