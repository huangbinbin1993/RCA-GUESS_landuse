///////////////////////////////////////////////////////////////////////////////////////
/// \file archive.cpp
/// \brief Classes to make (de)serializing to/from streams convenient
///
/// $Date: 2014-05-08 10:35:33 +0200 (Thu, 08 May 2014) $
///
///////////////////////////////////////////////////////////////////////////////////////

#include "config.h"
#include "archive.h"

ArchiveInStream::ArchiveInStream(std::istream& strm)
	: in(strm) {
}

bool ArchiveInStream::save() const {
	return false;
}

void ArchiveInStream::transfer(char* s, std::streamsize n) {
	in.read(s, n);
}

ArchiveOutStream::ArchiveOutStream(std::ostream& strm)
	: out(strm) {
}

bool ArchiveOutStream::save() const {
	return true;
}

void ArchiveOutStream::transfer(char* s, std::streamsize n) {
	out.write(s, n);
}
