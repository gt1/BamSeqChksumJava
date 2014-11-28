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
public class BamParser
{
	enum bamparserstate_enum
	{
		read_len,
		read_record
	};

	public static final int FPAIRED = (1 << 0);
	public static final int FPROPER_PAIR = (1 << 1);
	public static final int FUNMAP = (1 << 2);
	public static final int FMUNMAP = (1 << 3);
	public static final int FREVERSE = (1 << 4);
	public static final int FMREVERSE = (1 << 5);
	public static final int FREAD1 = (1 << 6);
	public static final int FREAD2 = (1 << 7);
	public static final int FSECONDARY = (1 << 8);
	public static final int FQCFAIL = (1 << 9);
	public static final int FDUP = (1 << 10);
	public static final int FSUPPLEMENTARY = (1 << 11);
	public static final int FCIGAR32 = (1 << 15);
	
	public static boolean isPaired(int flags) { return (flags & FPAIRED) != 0; }
	public static boolean isProperPair(int flags) { return (flags & FPROPER_PAIR) != 0; }
	public static boolean isUnmapped(int flags) { return (flags & FUNMAP) != 0; }
	public static boolean isMapped(int flags) { return (flags & FUNMAP) == 0; }
	public static boolean isNextUnmapped(int flags) { return (flags & FMUNMAP) != 0; }
	public static boolean isNextMapped(int flags) { return (flags & FMUNMAP) == 0; }
	public static boolean isReverse(int flags) { return (flags & FREVERSE) != 0; }
	public static boolean isNextReverse(int flags) { return (flags & FMREVERSE) != 0; }
	public static boolean isRead1(int flags) { return (flags & FREAD1) != 0; }
	public static boolean isRead2(int flags) { return (flags & FREAD2) != 0; }
	public static boolean isSecondary(int flags) { return (flags & FSECONDARY) != 0; }
	public static boolean isQCFail(int flags) { return (flags & FQCFAIL) != 0; }
	public static boolean isDup(int flags) { return (flags & FDUP) != 0; }
	public static boolean isSupplementary(int flags) { return (flags & FSUPPLEMENTARY) != 0; }
	public static boolean isCigar32(int flags) { return (flags & FCIGAR32) != 0; }

	static final byte decarray[] =
	{
		'=','A','C','M',
		'G','R','S','V',
		'T','W','Y','H',
		'K','D','B','N'
	};
	
	static final byte rcdecarray[] =
	{
		'=','T','G','K',
		'C','Y','S','B',
		'A','W','R','D',
		'M','H','V','N'
	};
	
	private boolean haveHeader;
	private	BamHeaderParser BHP;
	private BamHeaderParser.AddBlockResult headerabr;
	
	private int recordlen;
	private int recordlenread;
	
	private int recordread;
	private byte recordtext[];
	
	bamparserstate_enum state;

	public BamParser()
	{
		haveHeader = false;
		BHP = new BamHeaderParser();
		headerabr = new BamHeaderParser.AddBlockResult(false,0);
		state = bamparserstate_enum.read_len;
		recordlen = 0;
		recordlenread = 0;
	}
	
	public BamHeaderParser getHeader()
	{
		return BHP;
	}
	
	public static int decodeLE32(byte [] B, int offset)
	{
		return
			((B[offset+0] & 0xFF) << 0) |
			((B[offset+1] & 0xFF) << 8) |
			((B[offset+2] & 0xFF) << 16) |
			((B[offset+3] & 0xFF) << 24);
	}
	
	public static int getRefId(byte [] B, int offset, int len) throws Exception
	{
		if ( 0 + 4 > len )
			throw new Exception("BamParser.getRefId(): out of range");
		return decodeLE32(B,offset+0);
	}

	public static int getPos(byte [] B, int offset, int len) throws Exception
	{
		if ( 4 + 4 > len )
			throw new Exception("BamParser.getPos(): out of range");
		return decodeLE32(B,offset+4);
	}
	
	public static int getLReadName(byte [] B, int offset, int len) throws Exception
	{
		if ( 8 + 1 > len )
			throw new Exception("BamParser.getLReadName(): out of range");
		return (B[offset+8] & 0xFF);
	}

	public static int getMapQ(byte [] B, int offset, int len) throws Exception
	{
		if ( 9 + 1 > len )
			throw new Exception("BamParser.getMapQ(): out of range");
		return (B[offset+9] & 0xFF);
	}

	public static int getBin(byte [] B, int offset, int len) throws Exception
	{
		if ( 10 + 2 > len )
			throw new Exception("BamParser.getBin(): out of range");
		return (B[offset+10] & 0xFF) | ((B[offset+11] & 0xFF)<<8);
	}

	public static int getNCigar(byte [] B, int offset, int len) throws Exception
	{
		if ( 12 + 2 > len )
			throw new Exception("BamParser.getNCigar(): out of range");
		return (B[offset+12] & 0xFF) | ((B[offset+13] & 0xFF)<<8);
	}

	public static int getFlags(byte [] B, int offset, int len) throws Exception
	{
		if ( 14 + 2 > len )
			throw new Exception("BamParser.getFlags(): out of range");
		return (B[offset+14] & 0xFF) | ((B[offset+15] & 0xFF)<<8);
	}

	public static int getLSeq(byte [] B, int offset, int len) throws Exception
	{
		if ( 16 + 4 > len )
			throw new Exception("BamParser.getLSeq(): out of range");
		return decodeLE32(B,offset+16);
	}

	public static int getNextRefId(byte [] B, int offset, int len) throws Exception
	{
		if ( 20 + 4 > len )
			throw new Exception("BamParser.getNextRefId(): out of range");
		return decodeLE32(B,offset+20);
	}

	public static int getNextPos(byte [] B, int offset, int len) throws Exception
	{
		if ( 24 + 4 > len )
			throw new Exception("BamParser.getNextPos(): out of range");
		return decodeLE32(B,offset+24);
	}

	public static int getTLen(byte [] B, int offset, int len) throws Exception
	{
		if ( 28 + 4 > len )
			throw new Exception("BamParser.getTLen(): out of range");
		return decodeLE32(B,offset+28);
	}
	
	public static int getReadNameOffset()
	{
		return 32;
	}

	public static class DecodedSequence
	{
		byte [] sequence;
		int length;
	}

	public static String getReadName(byte [] B, int offset, int len) throws Exception
	{
		int lreadname = getLReadName(B,offset,len);

		if ( getReadNameOffset() + lreadname > len )
			throw new Exception("BamParser.getReadName(): out of range");
			
		return new String(B,offset+getReadNameOffset(),lreadname-1);
	}

	// get read name including terminating null byte
	public static void getReadName(byte [] B, int offset, int len, DecodedSequence seq) throws Exception
	{
		int lreadname = getLReadName(B,offset,len);

		if ( getReadNameOffset() + lreadname > len )
			throw new Exception("BamParser.getReadName(): out of range");
			
		if ( seq.sequence == null || seq.sequence.length < lreadname )
			seq.sequence = new byte[lreadname];
		
		int sub = getReadNameOffset() + offset;
		for ( int i = 0; i < lreadname; ++i, ++sub )
			seq.sequence[i] = B[sub];
			
		seq.length = lreadname;
	}
	
	public static int getCigarOffset(byte [] B, int offset, int len) throws Exception
	{
		return getReadNameOffset() + getLReadName(B,offset,len);
	}

	public static int getSeqOffset(byte [] B, int offset, int len) throws Exception
	{
		return getCigarOffset(B,offset,len) + 4 * getNCigar(B,offset,len);
	}

	public static int getQualOffset(byte [] B, int offset, int len) throws Exception
	{
		return getSeqOffset(B,offset,len) + ((getLSeq(B,offset,len) + 1)>>1);
	}

	public static int getAuxOffset(byte [] B, int offset, int len) throws Exception
	{
		return getQualOffset(B,offset,len) + getLSeq(B,offset,len);
	}
	
	
	private static byte decodeSeqSymbol(int b)
	{
		return decarray[b];
	}

	private static byte decodeSeqSymbolRC(int b)
	{
		return rcdecarray[b];
	}
	
	public static void decodeSequence(byte [] B, int offset, int len, DecodedSequence seq) throws Exception
	{
		int seqoff = getSeqOffset(B,offset,len);
		int lseq = getLSeq(B,offset,len);
		if ( seqoff + ((lseq+1)>>1) > len )
			throw new Exception("BamParser.decodeSequence(): out of range");
		if ( seq.sequence == null || seq.sequence.length < lseq )
			seq.sequence = new byte[lseq];
		seq.length = lseq;
		
		int o = 0;
		int z = offset + seqoff;
		int zz = z + (lseq>>1);
		for ( ; z != zz; ++z )
		{
			seq.sequence[o++] = decodeSeqSymbol((B[z] >> 4)&0xF);
			seq.sequence[o++] = decodeSeqSymbol((B[z] >> 0)&0xF);
		}
		if ( (lseq&1) != 0 )
		{
			seq.sequence[o++] = decodeSeqSymbol((B[z] >> 4)&0xF);
		}
	}

	public static void decodeSequenceRC(byte [] B, int offset, int len, DecodedSequence seq) throws Exception
	{
		int seqoff = getSeqOffset(B,offset,len);
		int lseq = getLSeq(B,offset,len);
		if ( seqoff + ((lseq+1)>>1) > len )
			throw new Exception("BamParser.decodeSequence(): out of range");
		if ( seq.sequence == null || seq.sequence.length < lseq )
			seq.sequence = new byte[lseq];
		seq.length = lseq;
		
		int o = lseq;
		int z = offset + seqoff;
		int zz = z + (lseq>>1);
		for ( ; z != zz; ++z )
		{
			seq.sequence[--o] = decodeSeqSymbolRC((B[z] >> 4)&0xF);
			seq.sequence[--o] = decodeSeqSymbolRC((B[z] >> 0)&0xF);
		}
		if ( (lseq&1) != 0 )
		{
			seq.sequence[--o] = decodeSeqSymbolRC((B[z] >> 4)&0xF);
		}
	}

	public static void decodeQuality(byte [] B, int offset, int len, DecodedSequence seq) throws Exception
	{
		int qualoff = getQualOffset(B,offset,len);
		int lseq = getLSeq(B,offset,len);
		if ( qualoff + lseq > len )
			throw new Exception("BamParser.decodeQuality(): out of range");
		if ( seq.sequence == null || seq.sequence.length < lseq )
			seq.sequence = new byte[lseq];
		seq.length = lseq;
		
		if ( lseq != 0 )
		{
			if ( (B[offset+qualoff] & 0xFF) == 255 )
				for ( int i = 0; i < lseq; ++i )
					seq.sequence[i] = B[offset+qualoff+i];
			else
				for ( int i = 0; i < lseq; ++i )
					seq.sequence[i] = (byte)(B[offset+qualoff+i] + 33);
		}
	}

	public static void decodeQualityRC(byte [] B, int offset, int len, DecodedSequence seq) throws Exception
	{
		int qualoff = getQualOffset(B,offset,len);
		int lseq = getLSeq(B,offset,len);
		if ( qualoff + lseq > len )
			throw new Exception("BamParser.decodeQuality(): out of range");
		if ( seq.sequence == null || seq.sequence.length < lseq )
			seq.sequence = new byte[lseq];
		seq.length = lseq;
		
		if ( lseq != 0 )
		{
			if ( (B[offset+qualoff] & 0xFF) == 255 )
			{
				int o = lseq;
				for ( int i = 0; i < lseq; ++i )
					seq.sequence[--o] = B[offset+qualoff+i];
			}
			else
			{
				int o = lseq;
				for ( int i = 0; i < lseq; ++i )
					seq.sequence[--o] = (byte)(B[offset+qualoff+i] + 33);
			}
		}
	}
	
	private static int getPrimLengthByType(byte b) throws Exception
	{
		switch ( b )
		{
			case 'A': case 'c': case 'C': 
				return 1;
			case 's': case 'S': 
				return 2;
			case 'i': case 'I': 
				return 4;
			case 'f': 
				return 4;
			default:
				throw new Exception("Unable to handle type "+b+" in BamParser.getPrimLengthByType()");
		}
	}

	/**
	 * get length of auxiliary field at D
	 *
	 * @param D encoded auxiliary field
	 * @return length of auxiliary field at D
	 **/
	public static int getAuxLength(byte [] D, int offset) throws Exception
	{
		switch ( D[offset+2] )
		{
			case 'A': case 'c': case 'C': 
				return 2+1+1;
			case 's': case 'S': 
				return 2+1+2;
			case 'i': case 'I': 
				return 2+1+4;
			case 'f': 
				return 2+1+4;
			case 'Z':
			case 'H':
			{
				int len = 2+1;
				offset += len;
				while ( D[offset] != 0 )
				{
					len++;
					offset++;
				}
				len++;
				return len;
			}
			case 'B':
			{
				byte eltype = D[offset+3];
				int numel = decodeLE32(D,offset+4);
				return 2/*tag*/+1/*B*/+1/*type*/+4/* array length */+numel *getPrimLengthByType(eltype);
			}
			default:
			{
				throw new Exception("Unable to handle type "+D[2]+" in BamParser.getAuxLength()");
			}
		}
	}

	public void addBlock(byte [] B, int offset, int len, BamRecordHandler handler) throws Exception
	{
		if ( ! haveHeader )
		{
			haveHeader = BHP.addBlock(B,len,headerabr);
			offset = headerabr.offset;
			
			if ( haveHeader )
			{
				BHP.extractReadGroups();
				handler.setupReadGroups();
			}
		}

		if ( haveHeader && offset != len )
		{
			while ( offset != len )
			{
				switch ( state )
				{
					case read_len:
					{
						if ( recordlenread == 0 )
						{
							int rl = -1;
							while (
								len-offset >= 4 &&
								len-offset >= 4+(rl = decodeLE32(B,offset))
							)
							{
								handler.handleRecord(B,offset+4,rl);
								offset += 4+rl;
							}
						}
						while ( recordlenread < 4 && offset < len )
							recordlen |= (B[offset++]&0xFF) << (8*(recordlenread++));
						if ( recordlenread == 4 )
						{
							state = bamparserstate_enum.read_record;
							recordread = 0;
							if ( recordtext == null || recordtext.length < recordlen )
								recordtext = new byte[recordlen];
						}
						break;
					}
					case read_record:
					{
						while ( recordread < recordlen && offset < len )
							recordtext[recordread++] = B[offset++];
							
						if ( recordread == recordlen )
						{
							handler.handleRecord(recordtext,0,recordlen);
							state = bamparserstate_enum.read_len;
							recordlen = 0;
							recordlenread = 0;		
						}
						break;	
					}
				}
			}
		}		
	}
}
