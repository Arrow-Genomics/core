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

//Basic
#include <iostream>
#include <memory>
#include <string>
#include <fstream>
#include <array>
#include <string>

//SeqAn
#include <seqan/bam_io.h>

//Apache Arrow
#include <arrow/api.h>

using namespace seqan;
using namespace arrow;
using namespace std;

int main(int argc, char* argv[])
{   
    // Read the input file name
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

        // alignments and records.
        BamAlignmentRecord record;
        BamAlignment alignment;
	
	// read and write records
	readRecord(record, bamFileIn);
	writeRecord(bamFileOut, record);

    }
    catch (Exception const & e)
    {
        std::cout << "ERROR: " << e.what() << std::endl;
        return 1;
    }

  //Arrow Apache table representation for alignment data.

  //Memory Pool
	arrow::MemoryPool* pool = arrow::default_memory_pool();

	arrow::StringBuilder qName_builder(pool);
	arrow::UInt32Builder flag_builder(pool);
	arrow::Int32Builder  rID_builder(pool);
	arrow::Int32Builder  beginPos_builder(pool);
	arrow::UInt32Builder mapQ_builder(pool);
	//Cigars
	arrow::UInt32Builder Cigar_count_builder(pool);
	arrow::Int8Builder Cigar_operation_builder(pool);

	arrow::Int32Builder  rNextId_builder(pool);
	arrow::Int32Builder  pNext_builder(pool);
	arrow::Int32Builder  tLen_builder(pool);
	arrow::StringBuilder seq_builder(pool);
	arrow::StringBuilder qual_builder(pool);
	arrow::StringBuilder tags_builder(pool);

  //Extract data from alignment structure and Append values to builders.
  for (uint32_t i = 0; i < length(alignment); i++) {

	auto qNames = toCString(alignment[i].qName);
	auto flags = alignment[i].flag;	
	auto rIDs = alignment[i].rID; 	
	auto beginPoss = alignment[i].beginPos;
	auto mapQs = alignment[i].mapQ;
	//Cigars
	auto rNextIds = alignment[i].rNextId;	
	auto pNexts = alignment[i].pNext;
	auto tLens = alignment[i].tLen;	
	CharString cstr = alignment[i].seq;
	auto seqs = toCString(cstr);	
	auto quals = toCString(alignment[i].qual);	
	auto tagss = toCString(alignment[i].tags);
        
	std::vector<uint32_t> Cigar_count_str = {0};
	std::vector<int8_t> Cigar_operation_str = {' '};

	for (uint32_t j = 0; j < length(alignment[i].cigar); ++j){
	Cigar_count_str.push_back( alignment[i].cigar[j].count );
	Cigar_operation_str.push_back( alignment[i].cigar[j].operation );
	}

	qName_builder.Append(qNames);
	flag_builder.Append(flags);
	beginPos_builder.Append(beginPoss);
	rID_builder.Append(rIDs);
	mapQ_builder.Append(mapQs);
	//Cigars
	Cigar_count_builder.AppendValues(Cigar_count_str);
	Cigar_operation_builder.AppendValues(Cigar_operation_str);

	rNextId_builder.Append(rNextIds);
	pNext_builder.Append(pNexts);
	tLen_builder.Append(tLens);
	seq_builder.Append(seqs);
	qual_builder.Append(quals);
	tags_builder.Append(tagss);
  }

  //Schema
  std::vector<std::shared_ptr<arrow::Field>> schema_vector = {
	arrow::field("qNames", arrow::utf8(), false),
	arrow::field("flags", arrow::uint32()),
	arrow::field("rIDs", arrow::int32()),
	arrow::field("beginPoss", arrow::int32()),
	arrow::field("mapQs", arrow::uint32()),
	arrow::field("Cigar_count_str", arrow::uint32()),
	arrow::field("Cigar_operation_str", arrow::int8()),
	arrow::field("rNextIds", arrow::int32()),
	arrow::field("pNexts", arrow::int32()),
	arrow::field("tLens", arrow::int32()),
	arrow::field("seqs", arrow::utf8(), false),
	arrow::field("quals", arrow::utf8(), false),
	arrow::field("tagss", arrow::utf8(), false)
  };
  auto schema = std::make_shared<arrow::Schema>(schema_vector);

  // Create an array and finish the builder
  shared_ptr<arrow::Array> qName_array;
  qName_builder.Finish(&qName_array);

   shared_ptr<arrow::Array> flag_array;
  flag_builder.Finish(&flag_array);

  shared_ptr<arrow::Array> rID_array;
  rID_builder.Finish(&rID_array);

  shared_ptr<arrow::Array> beginPos_array;
  beginPos_builder.Finish(&beginPos_array);

  shared_ptr<arrow::Array> mapQ_array;
  mapQ_builder.Finish(&mapQ_array);

  shared_ptr<arrow::Array> Cigar_count_array;
  Cigar_count_builder.Finish(&Cigar_count_array);

  shared_ptr<arrow::Array> Cigar_operation_array;
  Cigar_operation_builder.Finish(&Cigar_operation_array);

  shared_ptr<arrow::Array> rNextId_array;
  rNextId_builder.Finish(&rNextId_array);

  shared_ptr<arrow::Array> pNext_array;
  pNext_builder.Finish(&pNext_array);

  shared_ptr<arrow::Array> tLen_array;
  tLen_builder.Finish(&tLen_array);

  shared_ptr<arrow::Array> seq_array;
  seq_builder.Finish(&seq_array);

  shared_ptr<arrow::Array> qual_array;
  qual_builder.Finish(&qual_array);

  shared_ptr<arrow::Array> tags_array;
  tags_builder.Finish(&tags_array);

  // Create the table
  shared_ptr<arrow::Table> arrow_alignments_table = arrow::Table::Make(schema, { qName_array, flag_array, rID_array, beginPos_array, mapQ_array, Cigar_count_array, Cigar_operation_array, rNextId_array, pNext_array, tLen_array, seq_array, qual_array, tags_array});


//Reading back data test
/*
const shared_ptr<arrow::Column>& column = arrow_alignments_table->column(0);
const shared_ptr<arrow::Column>& column1 = arrow_alignments_table->column(1);
const shared_ptr<arrow::Column>& column10 = arrow_alignments_table->column(10);


// Get the StringArray representation
auto qNameS = static_pointer_cast<arrow::StringArray>(column->data()->chunk(0));
auto ids = std::static_pointer_cast<arrow::UInt32Array>(column1->data()->chunk(0));
auto seqS = std::static_pointer_cast<arrow::StringArray>(column10->data()->chunk(0));

//auto seqS = std::static_pointer_cast<arrow::UInt32Array>(column1->data()->chunk(0));

// Iterate over all strings
  for (uint32_t i = 0; i < ids->length(); i++) {
    // Get the string in a zero-copy manner
    int length;
    uint32_t id = 0;
    const char* str = (const char*) qNameS->GetValue(i, &length);
    id = ids->Value(i);
    //....
    //...
    //..
    //.
    const char* str1 = (const char*) seqS->GetValue(i, &length);


    std::cout << "qName(" << i << "): " << str << std::endl;
    std::cout << "flags(" << i << "): "  << id << std::endl;
    std::cout << "seq(" << i << "): "  << str1 << std::endl;

    std::cout << std::endl;
   }
*/

    return 0;
}


