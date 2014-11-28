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
public class UnsignedNumber
{
	private int length;
	private long [] A;
	
	private static String format32(long i)
	{
		String s = Long.toHexString(i & 0xFFFFFFFFL);
		StringBuffer sb = new StringBuffer();
		for ( int j = 0; j < 8-s.length(); ++j )
			sb.append('0');
		return sb.toString() + s;
	}
	
	public void setDigit(int i, long v)
	{
		A[i] = v;
	}
	
	public String toString()
	{
		StringBuffer sb = new StringBuffer();
		for ( int i = 0; i < A.length; ++i )
			sb.append(format32(A[A.length-i-1]));
		return sb.toString();
	}
	
	UnsignedNumber(int length)
	{
		this.length = length;
		this.A = new long[length];
	}

	UnsignedNumber(int length, long v)
	{
		this.length = length;
		this.A = new long[length];
		this.A[0] = v & 0xFFFFFFFFL;
	}
	
	void add(UnsignedNumber O) throws Exception
	{
		long sum = 0;

		int i = 0;		
		for ( ; i < Math.min(length,O.length); ++i )
		{
			sum = (A[i] + O.A[i] + (sum >> 32));
			A[i] = sum & 0xFFFFFFFFL;
		}

		while ( (sum != 0) && (i < length) )
		{
			sum = (sum >> 32) + A[i];
			A[i++] = sum & 0xFFFFFFFFL;
		}
	}
	
	void addDigest(byte [] B, int offset, int llength) throws Exception
	{
		if ( (llength & 3) != 0 )
			throw new Exception("UnsignedNumber.addDigest: length of digest is not a multiple of 4 bytes");
		if ( (llength >> 2) > length )
			throw new Exception("UnsignedNumber.addDigest: digest is too long for number format");
			
		int o = offset + llength;
		int i = 0;
		long sum = 0;
		
		// scan digest right to left (BE), number left to right (LE)
		while ( o != offset )
		{
			// construct number in BE order
			long v = 0;
			v |= ((long)(B[--o] & 0xFF)) << 0;
			v |= ((long)(B[--o] & 0xFF)) << 8;
			v |= ((long)(B[--o] & 0xFF)) << 16;
			v |= ((long)(B[--o] & 0xFF)) << 24;
			// compute sum
			sum = A[i] + v + (sum>>32);
			// set digit
			A[i++] = sum & 0xFFFFFFFFL;
		}
		
		// handle possible carry
		while ( (sum != 0) && (i < length) )
		{
			sum = (sum >> 32) + A[i];
			A[i++] = sum & 0xFFFFFFFFL;
		}
	}
}
