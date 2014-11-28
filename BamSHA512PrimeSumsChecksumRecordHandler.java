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
public class BamSHA512PrimeSumsChecksumRecordHandler implements BamRecordHandler
{
	BamHeaderParser BHP;
	BamParser.DecodedSequence seq;
	BamParser.DecodedSequence qual;
	BamParser.DecodedSequence name;
	int [] tagstart;
	int [] taglen;
	java.security.MessageDigest sha512;
	byte [] bdigest;

	private static final java.math.BigInteger prime = getPrime();
	
	UnsignedNumber crc_all_name_flags_seq;
	UnsignedNumber crc_all_flags_seq;
	UnsignedNumber crc_all_flags_seq_qual;
	UnsignedNumber crc_all_flags_seq_tags;
	long allcnt;
	UnsignedNumber crc_pass_name_flags_seq;
	UnsignedNumber crc_pass_flags_seq;
	UnsignedNumber crc_pass_flags_seq_qual;
	UnsignedNumber crc_pass_flags_seq_tags;
	long passcnt;

	UnsignedNumber [] rg_crc_all_name_flags_seq;
	UnsignedNumber [] rg_crc_all_flags_seq;
	UnsignedNumber [] rg_crc_all_flags_seq_qual;
	UnsignedNumber [] rg_crc_all_flags_seq_tags;
	long [] rg_allcnt;
	UnsignedNumber [] rg_crc_pass_name_flags_seq;
	UnsignedNumber [] rg_crc_pass_flags_seq;
	UnsignedNumber [] rg_crc_pass_flags_seq_qual;
	UnsignedNumber [] rg_crc_pass_flags_seq_tags;
	long [] rg_passcnt;
	
	final int maskflags = BamParser.FPAIRED | BamParser.FREAD1 | BamParser.FREAD2;
	
	private void applyModPrime()
	{
		applyModPrime(crc_all_name_flags_seq);
		applyModPrime(crc_all_flags_seq);
		applyModPrime(crc_all_flags_seq_qual);
		applyModPrime(crc_all_flags_seq_tags);
		applyModPrime(crc_pass_name_flags_seq);
		applyModPrime(crc_pass_flags_seq);
		applyModPrime(crc_pass_flags_seq_qual);
		applyModPrime(crc_pass_flags_seq_tags);
		
		for ( int i = 0; i < rg_crc_all_name_flags_seq.length; ++i )
		{
			applyModPrime(rg_crc_all_name_flags_seq[i]);
			applyModPrime(rg_crc_all_flags_seq[i]);
			applyModPrime(rg_crc_all_flags_seq_qual[i]);
			applyModPrime(rg_crc_all_flags_seq_tags[i]);
			applyModPrime(rg_crc_pass_name_flags_seq[i]);
			applyModPrime(rg_crc_pass_flags_seq[i]);
			applyModPrime(rg_crc_pass_flags_seq_qual[i]);
			applyModPrime(rg_crc_pass_flags_seq_tags[i]);		
		}
	}
	
	private static void applyModPrime(UnsignedNumber U)
	{
		String s = U.toString();
		int i = 0;
		while ( i < s.length() && s.charAt(i) == '0' )
			++i;
		s = s.substring(i);
		
		java.math.BigInteger num = (s.length() == 0) ? new java.math.BigInteger("0") : new java.math.BigInteger(s,16);
		num = num.mod(prime);
		
		s = num.toString(16);
		int add = (512+32)/4 - s.length();
		StringBuffer addbuf = new StringBuffer();
		for ( i = 0; i < add; ++i )
			addbuf.append('0');
		s = addbuf.toString() + s;
				
		for ( i = 0; i < (512+32)/32; ++i )
		{
			int charsperword = 32/4;
			String t = s.substring(
				s.length() - (i+1)*charsperword,
				s.length() - (i+0)*charsperword
			);
						
			U.setDigit(i,Long.parseLong(t,16));
		}
	}
	
	public static UnsignedNumber [] getUnsignedNumberArray(int size, int length, long init)
	{
		UnsignedNumber [] A = new UnsignedNumber[size];
		for ( int i = 0; i < A.length; ++i )
			A[i] = new UnsignedNumber(length,init);
		return A;
	}

	private static java.math.BigInteger getPrime()
	{
		java.math.BigInteger BI = new java.math.BigInteger("1");
		BI = BI.shiftLeft(512);
		BI = BI.add(new java.math.BigInteger("75"));
		return BI;
	}
		
	private static int getNumberLength()
	{
		return 512/32 + 1;
	}
	
	public void setupReadGroups() throws Exception
	{
		int numrg = BHP.getNumReadGroups();
		rg_crc_all_name_flags_seq = getUnsignedNumberArray(numrg+1,getNumberLength(),0);
		rg_crc_all_flags_seq = getUnsignedNumberArray(numrg+1,getNumberLength(),0);
		rg_crc_all_flags_seq_qual = getUnsignedNumberArray(numrg+1,getNumberLength(),0);
		rg_crc_all_flags_seq_tags = getUnsignedNumberArray(numrg+1,getNumberLength(),0);
		rg_crc_pass_name_flags_seq =  getUnsignedNumberArray(numrg+1,getNumberLength(),0);
		rg_crc_pass_flags_seq = getUnsignedNumberArray(numrg+1,getNumberLength(),0);
		rg_crc_pass_flags_seq_qual = getUnsignedNumberArray(numrg+1,getNumberLength(),0);
		rg_crc_pass_flags_seq_tags = getUnsignedNumberArray(numrg+1,getNumberLength(),0);
		rg_allcnt = new long[numrg+1];
		rg_passcnt = new long[numrg+1];
	}

	public BamSHA512PrimeSumsChecksumRecordHandler(BamHeaderParser BHP) throws Exception
	{
		this.BHP = BHP;
		seq = new BamParser.DecodedSequence();
		qual = new BamParser.DecodedSequence();
		name = new BamParser.DecodedSequence();
		tagstart = new int[5];
		taglen = new int[5];
		sha512 = java.security.MessageDigest.getInstance("SHA-512");
		bdigest = new byte[512/8];
		crc_all_name_flags_seq = new UnsignedNumber(getNumberLength(),0);
		crc_all_flags_seq = new UnsignedNumber(getNumberLength(),0);
		crc_all_flags_seq_qual = new UnsignedNumber(getNumberLength(),0);
		crc_all_flags_seq_tags = new UnsignedNumber(getNumberLength(),0);
		allcnt = 0;
		crc_pass_name_flags_seq = new UnsignedNumber(getNumberLength(),0);
		crc_pass_flags_seq = new UnsignedNumber(getNumberLength(),0);
		crc_pass_flags_seq_qual = new UnsignedNumber(getNumberLength(),0);
		crc_pass_flags_seq_tags = new UnsignedNumber(getNumberLength(),0);
		passcnt = 0;
	}

	public void handleRecord(byte [] B, int offset, int length) throws Exception
	{
		int flags = BamParser.getFlags(B,offset,length);
		
		if ( ! ( BamParser.isSecondary(flags) || BamParser.isSupplementary(flags) ) )
		{	
			int flagslow = (flags & maskflags) & 0xFF;
			boolean isqcpass = ! BamParser.isQCFail(flags);
			
			if ( ((++allcnt) & ((1l << 24)-1)) == 0 )
				applyModPrime();
			
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
			sha512.reset();
			sha512.update(name.sequence,0,name.length);
			sha512.update((byte)flagslow);
			sha512.update(seq.sequence,0,seq.length);		
			sha512.digest(bdigest,0,64);
			crc_all_name_flags_seq.addDigest(bdigest,0,64);
			rg_crc_all_name_flags_seq[rgid+1].addDigest(bdigest,0,64);
			if ( isqcpass )
			{
				crc_pass_name_flags_seq.addDigest(bdigest,0,64);
				rg_crc_pass_name_flags_seq[rgid+1].addDigest(bdigest,0,64);
			}

			// flags + sequence			
			sha512.reset();
			sha512.update((byte)flagslow);
			sha512.update(seq.sequence,0,seq.length);		
			sha512.digest(bdigest,0,64);
			crc_all_flags_seq.addDigest(bdigest,0,64);
			rg_crc_all_flags_seq[rgid+1].addDigest(bdigest,0,64);
			if ( isqcpass )
			{
				crc_pass_flags_seq.addDigest(bdigest,0,64);
				rg_crc_pass_flags_seq[rgid+1].addDigest(bdigest,0,64);
			}
						
			// flags + sequence + quality
			sha512.reset();
			sha512.update((byte)flagslow);
			sha512.update(seq.sequence,0,seq.length);
			sha512.update(qual.sequence,0,qual.length);
			sha512.digest(bdigest,0,64);			
			crc_all_flags_seq_qual.addDigest(bdigest,0,64);
			rg_crc_all_flags_seq_qual[rgid+1].addDigest(bdigest,0,64);
			if ( isqcpass )
			{
				crc_pass_flags_seq_qual.addDigest(bdigest,0,64);
				rg_crc_pass_flags_seq_qual[rgid+1].addDigest(bdigest,0,64);
			}

			// flags + sequence + tags
			sha512.reset();
			sha512.update((byte)flagslow);
			sha512.update(seq.sequence,0,seq.length);
			for ( int i = 0; i < tagstart.length; ++i )
				if ( tagstart[i] != -1 )
					sha512.update(B, offset+tagstart[i], taglen[i]);
			sha512.digest(bdigest,0,64);

			crc_all_flags_seq_tags.addDigest(bdigest,0,64);
			rg_crc_all_flags_seq_tags[rgid+1].addDigest(bdigest,0,64);
			if ( isqcpass )
			{
				crc_pass_flags_seq_tags.addDigest(bdigest,0,64);
				rg_crc_pass_flags_seq_tags[rgid+1].addDigest(bdigest,0,64);
			}
			
			rg_allcnt[rgid+1]++;
			if ( isqcpass )
				rg_passcnt[rgid+1]++;
		}
	}
	
	private String truncateDigest(UnsignedNumber U)
	{
		String s = U.toString();
		return s.substring(s.length() - ((512/32)*8) );
	}

	public String toString()
	{
		StringBuffer SB = new StringBuffer();
		
		applyModPrime();
		
		SB.append("###\tset\tcount\t\tb_seq\tname_b_seq\tb_seq_qual\tb_seq_tags(BC,FI,QT,RT,TC)\n");
		SB.append("all\tall\t" +
			allcnt + "\t\t" +
			truncateDigest(crc_all_flags_seq) + "\t" +
			truncateDigest(crc_all_name_flags_seq) + "\t" +
			truncateDigest(crc_all_flags_seq_qual) + "\t" +
			truncateDigest(crc_all_flags_seq_tags) + "\n"
		);
		SB.append("all\tpass\t" +
			passcnt + "\t\t" +
			truncateDigest(crc_pass_flags_seq) + "\t" +
			truncateDigest(crc_pass_name_flags_seq) + "\t" +
			truncateDigest(crc_pass_flags_seq_qual) + "\t" +
			truncateDigest(crc_pass_flags_seq_tags) + "\n"		
		);

		if ( rg_crc_all_name_flags_seq.length > 1 )	
			for ( int i = 0; i < rg_crc_all_name_flags_seq.length; ++i )
			{
				SB.append(
					BHP.getReadGroupId(i-1)+"\tall\t" +
					rg_allcnt[i] + "\t\t" +
					truncateDigest(rg_crc_all_flags_seq[i]) + "\t" +
					truncateDigest(rg_crc_all_name_flags_seq[i]) + "\t" +
					truncateDigest(rg_crc_all_flags_seq_qual[i]) + "\t" +
					truncateDigest(rg_crc_all_flags_seq_tags[i]) + "\n"			
				);		
				SB.append(
					BHP.getReadGroupId(i-1)+"\tpass\t" +
					rg_passcnt[i] + "\t\t" +
					truncateDigest(rg_crc_pass_flags_seq[i]) + "\t" +
					truncateDigest(rg_crc_pass_name_flags_seq[i]) + "\t" +
					truncateDigest(rg_crc_pass_flags_seq_qual[i]) + "\t" +
					truncateDigest(rg_crc_pass_flags_seq_tags[i]) + "\n"			
				);		
			}
				
		return SB.toString();		
	}
}
