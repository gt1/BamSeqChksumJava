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
public class BamSeqChksum
{
	private static boolean startsWith(String a, String prefix)
	{
		return
			a.length() >= prefix.length() && a.substring(0,prefix.length()).equals(prefix);
	}

	public static void main(String [] args) throws Exception
	{
		try
		{
			BgzfInput bgzfin = new BgzfInput(System.in);
			BgzfInput.GetBlockResult GBR = new BgzfInput.GetBlockResult();
			boolean done = false;
			BamParser BP = new BamParser(true,20);
			
			String hash = "crc32prod";
			for ( int i = 0; i < args.length; ++i )
				if ( startsWith(args[i],"hash=") )
				{
					hash = args[i].substring(new String("hash=").length());
				}

			BamRecordHandler BRH = null;
			if ( hash.equals("crc32prod") )
				BRH = new BamCRC32ProdChecksumRecordHandler(BP.getHeader());
			else if ( hash.equals("sha512") )
				BRH = new BamSHA512ChecksumRecordHandler(BP.getHeader());
			else if ( hash.equals("sha512primesums512") )
				BRH = new BamSHA512PrimeSumsChecksumRecordHandler(BP.getHeader());
			else
				throw new Exception("Unknown hash type "+hash);
			
			while ( ! (done = (!bgzfin.getNextBlock(GBR))) )
				BP.addBlock(GBR.block,0,GBR.blocksize,BRH);

			System.out.print(BRH);
		}
		catch(Exception ex)
		{
			System.err.println(ex);
			System.exit(1);
		}
	}
}
