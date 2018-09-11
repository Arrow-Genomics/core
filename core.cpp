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
#include <chrono>

//SeqAn
#include <seqan/bam_io.h>

//Apache Arrow
#include <arrow/api.h>

using namespace seqan;
using namespace arrow;
using namespace std;

typedef std::chrono::high_resolution_clock Clock;

////////////////////////////////////////////////////////////////////////
/* create_table_alignment() to create arrow table.
// 
// 
*/ 
////////////////////////////////////////////////////////////////////////
shared_ptr<arrow::Table> create_table_alignment(void)
{
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

	std::vector<uint32_t> Cigar_count_str;
	std::vector<int8_t> Cigar_operation_str;

        Cigar_count_str.push_back(length(alignment[i].cigar));
        Cigar_operation_str.push_back('x');

	for (uint32_t j = 0; j < length(alignment[i].cigar); j++){
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
  return move(arrow::Table::Make(schema, { qName_array, flag_array, rID_array, beginPos_array, mapQ_array, Cigar_count_array, Cigar_operation_array, rNextId_array, pNext_array, tLen_array, seq_array, qual_array, tags_array}));

}

////////////////////////////////////////////////////////////////////////
/* Main funtion() of core app.
// 
// 
*/ 
////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{   

long int timeVar,timeOverhead = 0;

    // Read the input file name
    if (argc > 1) {
        cout << "Parsing SAM/BAM File...-> " << argv[1] << endl; 
    } else {
        cout << "No file name entered. Exiting...";
        return -1;
    }

    CharString bamFileName = argv[1]; 

    auto t1 = Clock::now();  
   
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
        //writeHeader(bamFileOut, header);

        // alignments and records.
        BamAlignmentRecord record;
        BamAlignment alignment;
	
	// read and write records
	readRecord(record, bamFileIn);
	//writeRecord(bamFileOut, record);

        auto t2 = Clock::now(); timeVar=std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();

        writeHeader(bamFileOut, header);
	writeRecord(bamFileOut, record);

        auto t3 = Clock::now(); 
	
	timeOverhead = std::chrono::duration_cast<std::chrono::microseconds>(t3 - t1).count();
        std::cout << "Parsing, reading and Writing file time: " 
              << timeOverhead
              << " microseconds" << std::endl; 

    }
    catch (Exception const & e)
    {
        std::cout << "ERROR: " << e.what() << std::endl;
        return 1;
    }

    auto t4 = Clock::now(); 

//Arrow Apache table representation for alignment data.
    shared_ptr<arrow::Table> arrow_alignments_table = create_table_alignment();


//Reading back data test
	const shared_ptr<arrow::Column>& column = arrow_alignments_table->column(0);
	const shared_ptr<arrow::Column>& column1 = arrow_alignments_table->column(1);
	const shared_ptr<arrow::Column>& column2 = arrow_alignments_table->column(2);
	const shared_ptr<arrow::Column>& column3 = arrow_alignments_table->column(3);
	const shared_ptr<arrow::Column>& column4 = arrow_alignments_table->column(4);
	const shared_ptr<arrow::Column>& column5 = arrow_alignments_table->column(5);
	const shared_ptr<arrow::Column>& column6 = arrow_alignments_table->column(6);
	const shared_ptr<arrow::Column>& column7 = arrow_alignments_table->column(7);
	const shared_ptr<arrow::Column>& column8 = arrow_alignments_table->column(8);
	const shared_ptr<arrow::Column>& column9 = arrow_alignments_table->column(9);
	const shared_ptr<arrow::Column>& column10 = arrow_alignments_table->column(10);
	const shared_ptr<arrow::Column>& column11 = arrow_alignments_table->column(11);
	const shared_ptr<arrow::Column>& column12 = arrow_alignments_table->column(12);

// Get the StringArray representation
	auto qNameS = static_pointer_cast<arrow::StringArray>(column->data()->chunk(0));
	auto flagS = std::static_pointer_cast<arrow::UInt32Array>(column1->data()->chunk(0));
	auto rIDS = std::static_pointer_cast<arrow::Int32Array>(column2->data()->chunk(0));
	auto beginPosS = std::static_pointer_cast<arrow::Int32Array>(column3->data()->chunk(0));
	auto mapQS = std::static_pointer_cast<arrow::UInt32Array>(column4->data()->chunk(0));
	auto Cigar_countS = std::static_pointer_cast<arrow::UInt32Array>(column5->data()->chunk(0));
	auto Cigar_operationS = std::static_pointer_cast<arrow::Int8Array>(column6->data()->chunk(0));
	auto rNextIdS = std::static_pointer_cast<arrow::Int32Array>(column7->data()->chunk(0));
	auto pNextS = std::static_pointer_cast<arrow::Int32Array>(column8->data()->chunk(0));
	auto tLenS = std::static_pointer_cast<arrow::Int32Array>(column9->data()->chunk(0));
	auto seqS = std::static_pointer_cast<arrow::StringArray>(column10->data()->chunk(0));
	auto qualS = std::static_pointer_cast<arrow::StringArray>(column11->data()->chunk(0));
	auto tagsS = std::static_pointer_cast<arrow::StringArray>(column12->data()->chunk(0));

	static uint32_t offset = Cigar_countS->Value(0);
	static uint32_t start_pt = 1;

        BamAlignmentRecord record_arrow;
        clear(alignment);

//Iterate over all strings
  for (uint32_t i = 0; i < qNameS->length(); i++) {

// Get the alignment strings in a zero-copy manner
    int length;

    const char* buf = reinterpret_cast<const char*>(qNameS->GetValue(i, &length));
    record_arrow.qName = std::string(buf, length);

    record_arrow.flag = flagS->Value(i);
    record_arrow.rID = rIDS->Value(i);
    record_arrow.beginPos = beginPosS->Value(i);
    record_arrow.mapQ = mapQS->Value(i);

    // CIGAR
    CigarElement<> element;

    for (uint32_t j = start_pt; j < offset+1; j++){
        
        element.count = Cigar_countS->Value(j);
        element.operation =Cigar_operationS->Value(j);
        appendValue(record_arrow.cigar, element);
    }

    start_pt = offset + 2;
    offset = offset + Cigar_countS->Value(offset+1) +1;

    record_arrow.rNextId = rNextIdS->Value(i);
    record_arrow.pNext = pNextS->Value(i);
    record_arrow.tLen = tLenS->Value(i);

    buf = reinterpret_cast<const char*>(seqS->GetValue(i, &length));
    record_arrow.seq = std::string(buf, length);

    buf = reinterpret_cast<const char*>(qualS->GetValue(i, &length));
    record_arrow.qual = std::string(buf, length);

    buf = reinterpret_cast<const char*>(tagsS->GetValue(i, &length));
    record_arrow.tags= std::string(buf, length);

    appendValue(alignment, record_arrow);
    clear(record_arrow);
   }

//Writing back alignments from arrow memory
writeRecord(bamFileOut, record_arrow);

   auto t5 = Clock::now(); 
   std::cout << "Parsing, reading and Writing file through Arrow table time: " 
             << std::chrono::duration_cast<std::chrono::microseconds>(t5 - t4).count() + timeVar
             << " microseconds" << std::endl; 

   float Overhead =  ((float)timeOverhead / (std::chrono::duration_cast<std::chrono::microseconds>(t5 - t4).count() + timeVar))*100.0;
   std::cout << "Total overhead: " << Overhead << " %" << std::endl; 

 
    return 0;
}


