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
public class BamHeaderParser
{
	public static class ReferenceSequenceInfo
	{
		String name;
		int length;
		
		ReferenceSequenceInfo(String name, int length) 
		{
			this.name = name;
			this.length = length;
		}
	}

	enum state_enum {
		readmagic,
		readtextlength,
		readtext,
		readnumref,
		readrefnamelen,
		readrefnametext,
		readreflen,
		done
	};
	
	state_enum state;
	byte magic[];
	int magicread;

	int textlen;
	int textlenread;

	byte text[];
	int textread;
	
	int numref;
	int numrefread;
	int numrefparsed;
	
	int refnamelengthread;
	int refnamelength;
	int refnameread;
	byte refname[];
	int reflengthread;
	int reflen;
	
	public ReferenceSequenceInfo refseqs[];
	public ReadGroup readgroups[];
	int rgmask;
	int rghashcnt[];
	int rghashids[];
	
	public BamHeaderParser()
	{
		state = state_enum.readmagic;
		magicread = 0;
		magic = new byte[4];
		textlen = 0;
		textlenread = 0;
		textread = 0;
		numref = 0;
		numrefread = 0;
		numrefparsed = 0;
		rgmask = 0;
	}
	
	public int getNumReadGroups()
	{
		return readgroups.length;
	}
	
	public String getReadGroupId(int id)
	{
		if ( id < 0 || id >= readgroups.length )
			return new String();
		else
			return readgroups[id].getId();
	}
	
	public static class AddBlockResult
	{
		boolean finished;
		int offset;
		
		AddBlockResult(boolean finished, int offset)
		{
			this.finished = finished;
			this.offset = offset;
		}
	}
	
	public static class ReadGroupKeyValuePair
	{
		byte ida;
		byte idb;
		int vallow;
		int valhigh;
		
		ReadGroupKeyValuePair(byte ida, byte idb, int vallow, int valhigh)
		{
			this.ida = ida;
			this.idb = idb;
			this.vallow = vallow;
			this.valhigh = valhigh;
		}
	}

	public static final long FNV_prime = 0x100000001b3l;
	public static final long FNV_basis = 0xcbf29ce484222325l;
	
	public class ReadGroup
	{
		int idlow;
		int idhigh;
		ReadGroupKeyValuePair [] items;
				
		
		public int hashKey()
		{
			long h = FNV_basis;
			for ( int i = idlow; i < idhigh; ++i )
			{
				h = h ^ ((long)text[i]);
				h *= FNV_prime;
			}
			
			return (int)(h & 0x7FFFFFFFl);
		}
		
		
		ReadGroup(java.util.Vector<ReadGroupKeyValuePair> vitems)
		{
			boolean haveid = false;
			idlow = -1;
			idhigh = -1;
			
			for ( ReadGroupKeyValuePair item : vitems )
				if ( item.ida == 'I' && item.idb == 'D' )
				{
					haveid = true;
					idlow = item.vallow;
					idhigh = item.valhigh;	
				}
				
			items = haveid ? new ReadGroupKeyValuePair[vitems.size()-1] : new ReadGroupKeyValuePair[vitems.size()];

			int i = 0;
			for ( ReadGroupKeyValuePair item : vitems )
				if ( item.ida != 'I' || item.idb != 'D' )
					items[i++] = item;
		}
		
		public String getId()
		{
			return new String(text,idlow,idhigh-idlow);
		}
		
		public String toString()
		{
			StringBuffer SB = new StringBuffer();
			
			SB.append('@');
			SB.append('R');
			SB.append('G');
			
			if ( idlow != -1 )
			{
				SB.append('\t');
				SB.append('I');
				SB.append('D');
				SB.append(':');
				SB.append(new String(text,idlow,idhigh-idlow));
			}
			
			for ( int i = 0; i < items.length; ++i )
			{
				SB.append('\t');
				SB.append((char)items[i].ida);
				SB.append((char)items[i].idb);
				SB.append(':');
				SB.append(new String(text,items[i].vallow,items[i].valhigh-items[i].vallow));
			}
			
			return SB.toString();
		}
	}

	public static int hashKey(byte [] B, int off, int len)
	{
		long h = FNV_basis;
		for ( int i = off; i < off+len; ++i )
		{
			h = h ^ ((long)B[i]);
			h *= FNV_prime;
		}
		
		return (int)(h & 0x7FFFFFFFl);
	}
	
	public int getReadGroupId(byte [] B, int off, int len)
	{
		int h = hashKey(B,off,len) & rgmask;
		int c = rghashcnt[h+1]-rghashcnt[h];
		
		for ( int i = 0; i < c; ++i )
		{
			ReadGroup rg = readgroups[rghashids[rghashcnt[h]+i]];
			
			if ( len == rg.idhigh-rg.idlow )
			{
				boolean ok = true;
				
				for ( int j = 0; j < len; ++j )
					if ( text[rg.idlow+j] != B[off+j] )
					{
						ok = false;
						break;
					}
					
				if ( ok )
					return rghashids[rghashcnt[h]+i];
			}
		}
		
		return -1;
	}
	
	public void extractReadGroups() throws Exception
	{
		int offset = 0;
		java.util.Vector<ReadGroup> lreadgroups = new java.util.Vector<ReadGroup>();
		
		while ( offset != textlen )
		{
			int sol = offset;
			while ( offset < textlen && text[offset] != '\n' )
				++offset;
			int eol = offset;
			
			if ( eol-sol >= 4 && text[sol] == '@' && text[sol+1] == 'R' && text[sol+2] == 'G' && text[sol+3] == '\t' )
			{
				// System.err.println(new String(text,sol,eol-sol));
				
				int low = sol+4;
				
				java.util.Vector<ReadGroupKeyValuePair> items = new java.util.Vector<ReadGroupKeyValuePair>();
				
				while ( low != eol )
				{
					int high = low;
					while ( high != eol && text[high] != '\t' )
						++high;
						
					if ( 
						high-low >= 3 && text[low+2] == ':' 
					)
					{
						byte ida = text[low];
						byte idb = text[low+1];
						int vallow = low+3;
						int valhigh = high;
						
						items.add(
							new ReadGroupKeyValuePair(ida,idb,vallow,valhigh)
						);
						/*
						System.err.println(
							"" + (char)ida + (char)idb + "\t" + new String(text,vallow,valhigh-vallow));
						*/
					}
					
					if ( high != eol )
						++high;
					low = high;
				}
				
				ReadGroup rg = new ReadGroup(items);
				if ( rg.idlow < 0 )
					throw new Exception("invalid read group line "+new String(text,sol,eol-sol));
				lreadgroups.add(rg);
					
			}

			// skip over newline
			if ( offset < textlen )
				++offset;
		}
		
		readgroups = new ReadGroup[lreadgroups.size()];
		
		int s = 1;
		int shift = 0;
		while ( s < readgroups.length )
		{
			s <<= 1;
			shift += 1;
		}
		
		s <<= 2;
		shift += 2;
		rgmask = s-1;
		rghashcnt = new int[s+1];
				
		int o = 0;
		for ( ReadGroup rg : lreadgroups )
		{
			readgroups[o++] = rg;
			int hash = rg.hashKey() & rgmask;
			rghashcnt[hash]++;
			// System.err.println(rg.getId() + " " + hash );
		}
		
		// compute prefix sums over rghashcnt
		int sum = 0;
		int [] rghashcntcopy = new int[rghashcnt.length];
		for ( int i = 0; i < rghashcnt.length; ++i )
		{
			int t = rghashcnt[i];
			rghashcnt[i] = sum;
			sum += t;
			rghashcntcopy[i] = rghashcnt[i];
		}
		
		rghashids = new int[readgroups.length];
		for ( int i = 0; i < readgroups.length; ++i )
		{
			int hash = readgroups[i].hashKey() & rgmask;
			rghashids[rghashcntcopy[hash]++] = i;
		}

		/*		
		for ( int i = 0; i+1 < rghashcnt.length; ++i )
		{
			int off = rghashcnt[i];
			int len = rghashcnt[i+1]-off;
			
			for ( int j = 0; j < len; ++j )
			{
				System.err.println("*" + rghashids[off+j]);
			}
		}
		*/
	}
	
	public boolean addBlock(byte [] B, int len, AddBlockResult result) throws Exception
	{
		int offset = 0;
	
		while ( offset < len && state != state_enum.done )
		{
			switch ( state )
			{
				case readmagic:
				{
					while ( magicread < 4 && offset < len )
						magic[magicread++] = B[offset++];
						
					if ( magicread == 4 )
					{
						if ( 
							magic[0] == 'B' &&
							magic[1] == 'A' &&
							magic[2] == 'M' &&
							magic[3] == '\1'
						)
						{
							state = state_enum.readtextlength;
						}
						else
						{
							throw new Exception("Got wrong BAM magic");						
						}
					}
					break;
				}
				case readtextlength:
				{
					while ( textlenread < 4 && offset < len )
						textlen |= (B[offset++]&0xFF) << ((textlenread++)*8);
						
					if ( textlenread == 4 )
					{
						state = state_enum.readtext;
						text = new byte[textlen];
					}
					break;
				}
				case readtext:
				{
					while ( textread < textlen && offset < len )
						text[textread++] = B[offset++];
						
					if ( textread == textlen )
						state = state_enum.readnumref;
					break;
				}
				case readnumref:
				{
					while ( numrefread < 4 && offset < len )
						numref |= (B[offset++]&0xFF) << ((numrefread++)*8);
				
					if ( numrefread == 4 )
					{
						refseqs = new ReferenceSequenceInfo[numref];
					
						if ( numref == 0 )
							state = state_enum.done;
						else
						{
							state = state_enum.readrefnamelen;
							refnamelengthread = 0;
							refnamelength = 0;
						}
					}
					break;
				}
				case readrefnamelen:
				{
					while ( refnamelengthread < 4 && offset < len )
						refnamelength |= (B[offset++]&0xFF) << (refnamelengthread++ *8);
					
					if ( refnamelengthread == 4 )
					{
						state = state_enum.readrefnametext;
						refnameread = 0;
						if ( refname == null || refnamelength > refname.length )
							refname = new byte[refnamelength];
					}
					break;
				}
				case readrefnametext:
				{
					while ( refnameread < refnamelength && offset < len )
						refname[refnameread++] = B[offset++];
						
					if ( refnameread == refnamelength )
					{
						state = state_enum.readreflen;
						reflengthread = 0;
						reflen = 0;
					}
					break;
				}
				case readreflen:
				{
					while ( reflengthread < 4 && offset < len )
						reflen |= (B[offset++]&0xFF) << (8*(reflengthread++));
						
					if ( reflengthread == 4 )
					{
						refseqs[numrefparsed] = new ReferenceSequenceInfo(new String(refname,0,refnamelength-1),reflen);

						if ( ++numrefparsed == numref )
						{
							state = state_enum.done;
						}
						else
						{
							state = state_enum.readrefnamelen;
							refnamelengthread = 0;
							refnamelength = 0;						
						}
					}
				}
			}
		}
	
		result.finished = (state == state_enum.done);
		result.offset = offset;

		return result.finished;
	}
}
