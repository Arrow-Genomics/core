/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 * 
 * Author: Tanveer Ahmad
 */

#include <seqan/bam_io.h>


using namespace seqan;
using namespace std;

int main(int argc, char* argv[])
{   

    if (argc > 1) {
        cout << "Parsing SAM/BAM File...-> " << argv[1] << endl; 
    } else {
        cout << "No file name entered. Exiting...";
        return -1;
    }

    CharString bamFileName = argv[1]; 
       
    // Open input file, BamFileIn can read SAM and BAM files.
    BamFileIn bamFileIn;
    if (!open(bamFileIn, toCString(bamFileName)))
    {
        std::cerr << "ERROR: Could not open " << bamFileName << std::endl;
        return 1;
    }
    // Open output file, BamFileOut accepts also an ostream and a format tag.
    BamFileOut bamFileOut(context(bamFileIn), std::cout, Sam());

    try
    {
        // Copy header.
        BamHeader header;
        readHeader(header, bamFileIn);
        writeHeader(bamFileOut, header);

        //alignments and records.
        BamAlignmentRecord record;
        BamAlignment alignment;
	
	//read and write records
	readRecord(record, bamFileIn);
	writeRecord(bamFileOut, record);

    }
    catch (Exception const & e)
    {
        std::cout << "ERROR: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}



