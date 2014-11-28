/**
    BamSeqChksum
    Copyright (C) 2009-2014 German Tischler
    Copyright (C) 2011-2014 Genome Research Limited

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
**/
public class BamCRC32ProdChecksumRecordHandler implements BamRecordHandler
{
	BamHeaderParser BHP;
	BamParser.DecodedSequence seq;
	BamParser.DecodedSequence qual;
	BamParser.DecodedSequence name;
	int [] tagstart;
	int [] taglen;
	java.util.zip.CRC32 crc32;
	
	long crc_all_name_flags_seq;
	long crc_all_flags_seq;
	long crc_all_flags_seq_qual;
	long crc_all_flags_seq_tags;
	long allcnt;
	long crc_pass_name_flags_seq;
	long crc_pass_flags_seq;
	long crc_pass_flags_seq_qual;
	long crc_pass_flags_seq_tags;
	long passcnt;

	long [] rg_crc_all_name_flags_seq;
	long [] rg_crc_all_flags_seq;
	long [] rg_crc_all_flags_seq_qual;
	long [] rg_crc_all_flags_seq_tags;
	long [] rg_allcnt;
	long [] rg_crc_pass_name_flags_seq;
	long [] rg_crc_pass_flags_seq;
	long [] rg_crc_pass_flags_seq_qual;
	long [] rg_crc_pass_flags_seq_tags;
	long [] rg_passcnt;
	
	final int maskflags = BamParser.FPAIRED | BamParser.FREAD1 | BamParser.FREAD2;
	
	public static long [] getInitArray(int size, long init)
	{
		long [] A = new long[size];
		for ( int i = 0; i < A.length; ++i )
			A[i] = init;
		return A;
	}
	
	public void setupReadGroups() throws Exception
	{
		int numrg = BHP.getNumReadGroups();
		rg_crc_all_name_flags_seq = getInitArray(numrg+1,1);
		rg_crc_all_flags_seq = getInitArray(numrg+1,1);
		rg_crc_all_flags_seq_qual = getInitArray(numrg+1,1);
		rg_crc_all_flags_seq_tags = getInitArray(numrg+1,1);
		rg_crc_pass_name_flags_seq =  getInitArray(numrg+1,1);
		rg_crc_pass_flags_seq = getInitArray(numrg+1,1);
		rg_crc_pass_flags_seq_qual = getInitArray(numrg+1,1);
		rg_crc_pass_flags_seq_tags = getInitArray(numrg+1,1);
		rg_allcnt = getInitArray(numrg+1,0);
		rg_passcnt = getInitArray(numrg+1,0);
	}

	public BamCRC32ProdChecksumRecordHandler(BamHeaderParser BHP)
	{
		this.BHP = BHP;
		seq = new BamParser.DecodedSequence();
		qual = new BamParser.DecodedSequence();
		name = new BamParser.DecodedSequence();
		tagstart = new int[5];
		taglen = new int[5];
		crc32 = new java.util.zip.CRC32();
		crc_all_name_flags_seq = 1;
		crc_all_flags_seq = 1;
		crc_all_flags_seq_qual = 1;
		crc_all_flags_seq_tags = 1;
		allcnt = 0;
		crc_pass_name_flags_seq = 1;
		crc_pass_flags_seq = 1;
		crc_pass_flags_seq_qual = 1;
		crc_pass_flags_seq_tags = 1;
		passcnt = 0;
	}

	public void handleRecord(byte [] B, int offset, int length) throws Exception
	{
		int flags = BamParser.getFlags(B,offset,length);
		
		if ( ! ( BamParser.isSecondary(flags) || BamParser.isSupplementary(flags) ) )
		{	
			int flagslow = (flags & maskflags) & 0xFF;
			boolean isqcpass = ! BamParser.isQCFail(flags);
			
			allcnt += 1;
			if ( isqcpass )
				passcnt += 1;
			
			if ( BamParser.isReverse(flags) )
			{
				BamParser.decodeSequenceRC(B,offset,length,seq);
				BamParser.decodeQualityRC(B,offset,length,qual);
			}
			else
			{
				BamParser.decodeSequence(B,offset,length,seq);
				BamParser.decodeQuality(B,offset,length,qual);
			}
			
			BamParser.getReadName(B,offset,length,name);

			int rgid = -1;
			int suboff = BamParser.getAuxOffset(B,offset,length);
			for ( int i = 0; i < tagstart.length; ++i )
				tagstart[i] = -1;
			while ( suboff < length )
			{
				int ltaglen = BamParser.getAuxLength(B,offset+suboff);
				
				if ( 
					B[offset+suboff]   == 'R' && 
					B[offset+suboff+1] == 'G' && 
					B[offset+suboff+2] == 'Z' 
				)
				{
					int rgstart = offset+suboff+3;
					int rgend = rgstart;
					while ( B[rgend] != 0 )
						rgend++;
					rgid = BHP.getReadGroupId(B,rgstart,rgend-rgstart);
				}

				switch ( (char)B[offset+suboff] )
				{
					case 'B':
					{
						switch ( (char)B[offset+suboff+1] )
						{
							case 'C':
							{
								tagstart[0] = suboff;
								taglen[0] = ltaglen;
								break;
							}
							default:
								break;
						}			
					
						break;
					}
					case 'F':
					{
						switch ( (char)B[offset+suboff+1] )
						{
							case 'I':
							{
								tagstart[1] = suboff;
								taglen[1] = ltaglen;
								break;
							}
							default:
								break;
						}

						break;
					}
					case 'Q':
					{
						switch ( (char)B[offset+suboff+1] )
						{
							case 'T':
							{
								tagstart[2] = suboff;
								taglen[2] = ltaglen;
								break;
							}
							default:
								break;
						}			
					
						break;
					}
					case 'R':
					{
						switch ( (char)B[offset+suboff+1] )
						{
							case 'T':
							{
								tagstart[3] = suboff;
								taglen[3] = ltaglen;
								break;
							}
							default:
								break;
						}			
					
						break;
					}
					case 'T':
					{
						switch ( (char)B[offset+suboff+1] )
						{
							case 'C':
							{
								tagstart[4] = suboff;
								taglen[4] = ltaglen;
								break;
							}
							default:
								break;
						}			
						break;
					}
					default:
					{
						break;
					}
				}
				suboff += ltaglen;
			}

			// name + flags + sequence
			crc32.reset();
			crc32.update(name.sequence,0,name.length);
			crc32.update(flagslow);
			crc32.update(seq.sequence,0,seq.length);		
			long crcval = crc32.getValue();		

			crcval = crcval & 0x7FFFFFFFL;
			if ( crcval == 0 || crcval == 0x7FFFFFFFL )
				crcval = 1;
			crc_all_name_flags_seq *= crcval;
			crc_all_name_flags_seq %= 0x7FFFFFFFL;
			rg_crc_all_name_flags_seq[rgid+1] *= crcval;
			rg_crc_all_name_flags_seq[rgid+1] %= 0x7FFFFFFFL;
			if ( isqcpass )
			{
				crc_pass_name_flags_seq *= crcval;
				crc_pass_name_flags_seq %= 0x7FFFFFFFL;
				rg_crc_pass_name_flags_seq[rgid+1] *= crcval;
				rg_crc_pass_name_flags_seq[rgid+1] %= 0x7FFFFFFFL;
			}
			
			// flags + sequence
			crc32.reset();
			crc32.update(flagslow);
			crc32.update(seq.sequence,0,seq.length);		
			crcval = crc32.getValue();		
			crcval = crcval & 0x7FFFFFFFL;
			if ( crcval == 0 || crcval == 0x7FFFFFFFL )
				crcval = 1;
			crc_all_flags_seq *= crcval;
			crc_all_flags_seq %= 0x7FFFFFFFL;
			rg_crc_all_flags_seq[rgid+1] *= crcval;
			rg_crc_all_flags_seq[rgid+1] %= 0x7FFFFFFFL;
			if ( isqcpass )
			{
				crc_pass_flags_seq *= crcval;
				crc_pass_flags_seq %= 0x7FFFFFFFL;
				rg_crc_pass_flags_seq[rgid+1] *= crcval;
				rg_crc_pass_flags_seq[rgid+1] %= 0x7FFFFFFFL;
			}

			// flags + sequence + quality
			crc32.reset();
			crc32.update(flagslow);
			crc32.update(seq.sequence,0,seq.length);
			crc32.update(qual.sequence,0,qual.length);
			crcval = crc32.getValue();	
			crcval = crcval & 0x7FFFFFFFL;
			if ( crcval == 0 || crcval == 0x7FFFFFFFL )
				crcval = 1;
			crc_all_flags_seq_qual *= crcval;
			crc_all_flags_seq_qual %= 0x7FFFFFFFL;
			rg_crc_all_flags_seq_qual[rgid+1] *= crcval;
			rg_crc_all_flags_seq_qual[rgid+1] %= 0x7FFFFFFFL;
			if ( isqcpass )
			{
				crc_pass_flags_seq_qual *= crcval;
				crc_pass_flags_seq_qual %= 0x7FFFFFFFL;
				rg_crc_pass_flags_seq_qual[rgid+1] *= crcval;
				rg_crc_pass_flags_seq_qual[rgid+1] %= 0x7FFFFFFFL;
			}

			// flags + sequence + tags
			crc32.reset();
			crc32.update(flagslow);
			crc32.update(seq.sequence,0,seq.length);
			for ( int i = 0; i < tagstart.length; ++i )
				if ( tagstart[i] != -1 )
					crc32.update(B, offset+tagstart[i], taglen[i]);
			crcval = crc32.getValue();	
			crcval = crcval & 0x7FFFFFFFL;
			if ( crcval == 0 || crcval == 0x7FFFFFFFL )
				crcval = 1;
			crc_all_flags_seq_tags *= crcval;
			crc_all_flags_seq_tags %= 0x7FFFFFFFL;
			rg_crc_all_flags_seq_tags[rgid+1] *= crcval;
			rg_crc_all_flags_seq_tags[rgid+1] %= 0x7FFFFFFFL;
			if ( isqcpass )
			{
				crc_pass_flags_seq_tags *= crcval;
				crc_pass_flags_seq_tags %= 0x7FFFFFFFL;
				rg_crc_pass_flags_seq_tags[rgid+1] *= crcval;
				rg_crc_pass_flags_seq_tags[rgid+1] %= 0x7FFFFFFFL;
			}
			
			rg_allcnt[rgid+1] ++;
			if ( isqcpass )
			{
				rg_passcnt[rgid+1] ++;
			}
		}
	}

	public String toString()
	{
		StringBuffer SB = new StringBuffer();

		SB.append("###\tset\tcount\t\tb_seq\tname_b_seq\tb_seq_qual\tb_seq_tags(BC,FI,QT,RT,TC)\n");
		SB.append("all\tall\t" +
			allcnt + "\t\t" +
			Long.toHexString(crc_all_flags_seq) + "\t" +
			Long.toHexString(crc_all_name_flags_seq) + "\t" +
			Long.toHexString(crc_all_flags_seq_qual) + "\t" +
			Long.toHexString(crc_all_flags_seq_tags) + "\n"
		);
		SB.append("all\tpass\t" +
			passcnt + "\t\t" +
			Long.toHexString(crc_pass_flags_seq) + "\t" +
			Long.toHexString(crc_pass_name_flags_seq) + "\t" +
			Long.toHexString(crc_pass_flags_seq_qual) + "\t" +
			Long.toHexString(crc_pass_flags_seq_tags) + "\n"		
		);

		if ( rg_crc_all_name_flags_seq.length > 1 )		
			for ( int i = 0; i < rg_crc_all_name_flags_seq.length; ++i )
			{
				SB.append(
					BHP.getReadGroupId(i-1)+"\tall\t" +
					rg_allcnt[i] + "\t\t" +
					Long.toHexString(rg_crc_all_flags_seq[i]) + "\t" +
					Long.toHexString(rg_crc_all_name_flags_seq[i]) + "\t" +
					Long.toHexString(rg_crc_all_flags_seq_qual[i]) + "\t" +
					Long.toHexString(rg_crc_all_flags_seq_tags[i]) + "\n"			
				);		
				SB.append(
					BHP.getReadGroupId(i-1)+"\tpass\t" +
					rg_passcnt[i] + "\t\t" +
					Long.toHexString(rg_crc_pass_flags_seq[i]) + "\t" +
					Long.toHexString(rg_crc_pass_name_flags_seq[i]) + "\t" +
					Long.toHexString(rg_crc_pass_flags_seq_qual[i]) + "\t" +
					Long.toHexString(rg_crc_pass_flags_seq_tags[i]) + "\n"			
				);		
			}

		return SB.toString();		
	}
}
