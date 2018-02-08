// Adapted from the SPT-3G collaboration's G3Timestream class (spt3g_software/core/src/G3Timestream.cxx @ https://github.com/CMB-S4/spt3g_software).
// As per the license distributed with the above software, the following MUST be included in any further usage of this code:
/* --- BEGIN LICENSE STATEMENT ---
Copyright (c) 2015-2017 South Pole Telescope Collaboration, except where
 otherwise noted.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:
1. Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
SUCH DAMAGE.
   --- END LICENSE STATEMENT --- */

#include <iostream>
#include <cstdint> 	// int32_t, uint8_t, etc
#include <fstream>	// std::fstream
#include <vector> 	// std::vector
#include "FLAC/stream_encoder.h"


// DO NOT USE THESE IF THIS CODE GETS INCLUDED IN ANYTHING ELSE. THIS IS ONLY
// FOR YOUR CONVENIENCE.
#define N_INPUT_SAMPLES 1000
#define FLAC_COMPRESSION_LEVEL 9

static FLAC__StreamEncoderWriteStatus flac_encoder_write_cb(
    const FLAC__StreamEncoder *encoder, const FLAC__byte buffer[], size_t bytes,
    unsigned samples, unsigned current_frame, void *client_data) {
  // Since the client_data comes in as a void*, recast to the proper accessor
	std::vector<uint8_t> *outbuf = static_cast<std::vector<uint8_t>*>(client_data);
  // Copy the data that comes out of the encoding stream to outbuf
	outbuf->insert(outbuf->end(), buffer, buffer + bytes);
  // Return okay
	return FLAC__STREAM_ENCODER_WRITE_STATUS_OK;
}

int main(int argc, char **argv) {
	std::vector<int32_t> inbuf(N_INPUT_SAMPLES);
  std::vector<uint8_t> outbuf;
  const int32_t *chanmap[1];

  // Open the input file with a stream
  std::fstream fs_in("noise1.txt", std::ios::in);
  if (!fs_in.good()) { // Check if the file opened okay
  	std::cout << "Error: File 'INSERT_FILENAME_HERE' cannot be opened." << std::endl;
    // Try to close the file anyway, even if it didn't open successfully
    fs_in.close();
    return 0;
  }

  // Populate the input buffer with the data from the file
  for (size_t i=0; i<inbuf.size(); ++i) fs_in >> inbuf[i];

  // Close the file since we don't need it anymore
  fs_in.close();

  // Apply the bitshift filter (NOTE: IF THE DATA IS NORMALIZED TO VALUES OUTSIDE
  // OF THE RANGE REPRESENTABLE BY 24 BITS THIS WILL DESTROY INFO)
  //for (size_t i=0; i<inbuf.size(); ++i) inbuf[i] = ((inbuf[i] & 0x00FFFFFF) << 8) >> 8;

	// Assign the input data "channel" to the data in inbuf
  chanmap[0] = inbuf.data();

	// Yay FLAC encoding stuff! Set up new FLAC encoder stream
  FLAC__StreamEncoder *encoder = FLAC__stream_encoder_new();
  // Set one channel of data to encode
  FLAC__stream_encoder_set_channels(encoder, 1);
  // Tell the encoder the individual sample size
	FLAC__stream_encoder_set_bits_per_sample(encoder, 24);
  // From the G3Timestream doc - compression levels range from 0-9 (0 == no FLAC)
  FLAC__stream_encoder_set_compression_level(encoder, FLAC_COMPRESSION_LEVEL);
  // Set up the encoder stream and set the callback function (AFTER it encodes the data,
  // but BEFORE it tries to do anything else with it, it'll go to the callback)
  FLAC__stream_encoder_init_stream(encoder, flac_encoder_write_cb, NULL, NULL, NULL, (void*)(&outbuf));
  // Encode the stream
  FLAC__stream_encoder_process (encoder, chanmap, inbuf.size());
  // Finalize the stream
  FLAC__stream_encoder_finish(encoder);
  // Delete the stream and any resources it allocated
  FLAC__stream_encoder_delete(encoder);

	// Open the output file with a stream
	std::fstream fs_out("noise.flac", std::ios::out);
  if (!fs_out.good()) { // Check if the file opened okay
  	std::cout << "Error: File 'INSERT_FILENAME_HERE' cannot be opened." << std::endl;
    // Try to close the file anyway, even if it didn't open successfully
    fs_out.close();
    return 0;
}

  // Write the output data stream to file (space-delimimted, change however you see fit)
  for (size_t i=0; i<outbuf.size(); ++i) fs_out << outbuf[i] << " ";

  // Close the file
  fs_out.close();

	return 0;
}

//Covering my ass
#undef N_INPUT_SAMPLES
#undef FLAC_COMPRESSION_LEVEL
