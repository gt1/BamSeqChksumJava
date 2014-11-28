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
import java.io.DataInputStream;

public class BgzfInput
{
	private static int headersize = 18;
	private static int bgzfblockmax = 64*1024;
	private java.io.DataInputStream dis;
	
	private byte[] bgzfheader;
	private byte[] bgzfin;
	private byte[] bgzfout;

	private java.util.zip.Inflater inflater;

	BgzfInput(java.io.InputStream is)
	{
		dis = new java.io.DataInputStream(is);
		bgzfheader = new byte[headersize];
		bgzfin = new byte[bgzfblockmax];
		bgzfout = new byte[bgzfblockmax];
		inflater = new java.util.zip.Inflater(true);
	}
	
	public static class GetBlockResult
	{
		byte [] block;
		int blocksize;
	}
	
	boolean getNextBlock(GetBlockResult GBR) throws java.io.IOException, java.util.zip.DataFormatException
	{
		int read = dis.read(bgzfheader);
		
		switch ( read )
		{
			case 18: // header size
				break;
			case -1:
				return false;
			default:
				throw new java.io.IOException("Unable to read BGZF header");
		}
		
		if ( bgzfheader[0] != 31 || (bgzfheader[1]&0xFF) != 139 || bgzfheader[2] != 8 || bgzfheader[3] != 4 ||
			bgzfheader[10] != 6 || bgzfheader[11] != 0 ||
			bgzfheader[12] != 66 || bgzfheader[13] != 67 || bgzfheader[14] != 2 || bgzfheader[15] != 0
		)
		{
			throw new java.io.IOException("BgzfInput.getNextBlock(): invalid header data");
		}
	
		int cblocksize = 
			((bgzfheader[16]&0xFF)
			|
			((bgzfheader[17]&0xFF)<<8))+1;
		
		if ( cblocksize < 18 + 8 )
			throw new java.io.IOException("BgzfInput.getNextBlock(): invalid header data");
			
		// size of compressed data
		int payloadsize = cblocksize - (18 + 8);
		
		int inread = dis.read(bgzfin,0,payloadsize+8);
		
		inflater.reset();
		inflater.setInput(bgzfin,0,payloadsize);
		int inflated = inflater.inflate(bgzfout);
		
		if ( inflater.getTotalIn() != payloadsize )
			throw new java.io.IOException("BgzfInput.getNextBlock(): invalid compressed data");
		
		GBR.block = bgzfout;
		GBR.blocksize = inflated;
		
		return true;
	}
}
