"""
Bam.py

Pure python implementation of reading BAM(Binary sequence Alignment/Map format)

Author : Thomas Sunghoon Heo
"""

import gzip
import struct

BAM_MAGIC_STRING = "BAM\1"
BAM_EOF = "\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00"

def format_checker(fqname) :
	F = open(fqname,"rb")
	bytes = F.read(3)
	if bytes == "\x1f\x8b\x08" :
		F.seek(0)
		return gzip.GzipFile(fqname , fileobj=F)
	else:
		F.seek(0)
		return F

class BamHeader:
	def __init__(self):
		self.text = None
		self.references = {}
		
class BamRecord:
	def __init__(self):
		self.tid = None
		self.pos = None
		self.qname = None
		self.flag = None
		self.mapq = None
		self.seq = None
		self.tlen = None
		self.qual = None
		self.mpos = None
		self.mtid = None
		self.cigartuples = None
		self.aux = {}
		
class BamFile:
	def __init__(self, fname=None):
		self.fobj = None
		self.fname = None
		self.header = BamHeader()
	def open(self, fname=None):
		if self.fname == None:
			if fname != None:
				self.fname = fname
			else:
				raise Exception("fname parameter must be set")
		self.fobj = format_checker(self.fname)
		self.magic = self.fobj.read(4)
	def isBamFile(self):
		if self.fobj != None:
			if self.magic == BAM_MAGIC_STRING:
				return True
			else:
				return False
		else:
			Exception("Please open file before asking bam file")
	
	def readHeader(self):
		if self.fobj == None:
			raise Exception("File is not open for use")
			
		bdata = self.fobj.read(bytes_read)
		l_text = struct.unpack('<i', bdata)[0]
		bdata = self.fobj.read(l_text)
		self.header.text = bdata
		bdata = self.fobj.read(32/8)
		n_ref = struct.unpack('<i', bdata)[0]
		for i in range(n_ref):
			tid = i
			bdata = self.fobj.read(4)
			l_name = struct.unpack('<i', bdata)[0]
			name = self.fobj.read(l_name)
			bdata = self.fobj.read(4)
			l_ref = struct.unpack('<i', bdata)[0]
			self.header.references[tid] = name
		
	def readRecord(self):
		record = BamRecord()
		bdata = self.fobj.read(4)
		block_size = struct.unpack('<i', bdata)[0] # bytes
		consumed = 0
	
		bdata = self.fobj.read(4)
		consumed += 4
		ref_id = struct.unpack('<i', bdata)[0]
		record.tid = ref_id
		
		bdata = self.fobj.read(4)
		consumed += 4
		pos = struct.unpack('<i', bdata)[0]
		record.pos = pos
		
		bdata = self.fobj.read(1)
		consumed += 1
		l_read_name = struct.unpack('<B', bdata)[0] # uint8_t

		bdata = self.fobj.read(1)
		consumed += 1
		mapq = struct.unpack('<B', bdata)[0]
		record.mapq = mapq
		
		
		bdata = self.fobj.read(2)
		consumed += 2
		bai_index = struct.unpack('<H', bdata)[0]#uint16_t

		bdata = self.fobj.read(2)
		consumed += 2
		n_cigar_op = struct.unpack('<H', bdata)[0]

		bdata = self.fobj.read(2)
		consumed += 2
		flag = struct.unpack('<H', bdata)[0]
		record.flag = flag
		
		
		bdata = self.fobj.read(4)
		consumed += 4
		l_qseq = struct.unpack('<i', bdata)[0]

		bdata = self.fobj.read(4)
		consumed += 4
		next_ref_id = struct.unpack('<i', bdata)[0]
		record.mtid = next_ref_id
		
		bdata = self.fobj.read(4)
		consumed += 4
		next_ref_pos = struct.unpack('<i', bdata)[0]
		record.mpos = next_ref_pos
		
		bdata= self.fobj.read(4)
		consumed += 4
		tlen = struct.unpack('<i', bdata)[0]
		record.tlen = tlen
		
		bdata = self.fobj.read(l_read_name)
		consumed += l_read_name
		read_name = bdata
		record.qname = read_name
		
		
		cigartuples = []
		cigarletters="MIDNSHP=XB"
		for _ in range(n_cigar_op):
			bdata = self.fobj.read(4)
			consumed += 4
			cigar = struct.unpack('<I', bdata)[0]
			oplen = (cigar >> 4)
			op = cigarletters[(cigar & 0xf)]
		record.cigartuples = cigartuples
		
		seq = []
		seqletters = "=ACMGRSVTWYHKDBN"
		iters = (l_qseq + 1) / 2
		is_odd = True if l_qseq % 2 == 1 else False
		for _ in range(iters):
			bdata = self.fobj.read(1)
			consumed += 1
			encoded = struct.unpack('<B', bdata)[0]
			b1 = seqletters[(encoded & 0xf0) >> 4]
			b2 = seqletters[(encoded & 0xf)]
			seq.append(b1)
			seq.append(b2)

		if is_odd:
			seq.pop(-1)

		seq = "".join(seq)
		record.seq = seq
		

		qual = "".join([chr(ord(x) + 33) for x in self.fobj.read(l_qseq)])
		consumed += l_qseq
		record.qual = qual
		
		while consumed <= block_size:
			tag = self.fobj.read(2)
			consumed += 2
			val_type = self.fobj.read(1)
			if val_type == 'c':
				bdata = self.fobj.read(1)
				consumed += 1
				value = struct.unpack('<b', bdata)[0]
			elif val_type == 'C':
				bdata = self.fobj.read(1)
				consumed += 1
				value = struct.unpack('<B', bdata)[0]
			elif val_type == 'A':
				bdata = self.fobj.read(1)
				consumed+=1
				value = struct.unpack('<c', bdata)[0]
			elif val_type == 'f':
				bdata = self.fobj.read(4)
				consumed += 4
				value = struct.unpack('<f', bdata)[0]
			elif val_type == 'i':
				bdata = self.fobj.read(4)
				consumed += 4
				vlaue = struct.unpack('<i', bdata)[0]
			elif val_type == 'I':
				bdata = self.fobj.read(4)
				consumed += 4
				value = struct.unpack('<I', bdata)[0]
			elif val_type == 'S':
				bdata = self.fobj.read(2)
				consumed += 2
				value = struct.unpack('<h', bdata)[0]
			elif val_type == 's':
				bdata = self.fobj.read(2)
				consumed += 2
				value = struct.unpack('<H', bdata)[0]
			elif val_type == 'Z' or val_type == 'H':
				_tmp = []
				while 1:
					bdata = self.fobj.read(1)
					consumed += 1
					if bdata == '\x00':
						break
					_tmp.append(bdata)
				value = "".join(_tmp)	
			else:
				# B
				next_char = self.fobj.read(1)
				consumed += 1
				bdata = self.fobj.read(4)
				consumed += 4
				cnt = struct.unpack('<i', bdata)[0]
				_tmp = []
				for _ in range(cnt):
					if next_char == 'c':
						bdata = self.fobj.read(1)
						consumed += 1
						tmp = struct.unpack('<b',bdata)[0]
						_tmp.append(tmp)
					elif next_char == 'C':
						bdata = self.fobj.read(1)
						consumed += 1
						tmp = struct.unpack('<B', bdata)[0]
						_tmp.append(tmp)
					elif next_char == 's':
						bdata = self.fobj.read(2)
						consumed += 2
						tmp = struct.unpack('<h', bdata)[0]
						_tmp.append(tmp)
					elif next_char == 'S':
						bdata = self.fobj.read(2)
						consumed += 2
						tmp = struct.unpack('<H', bdata)[0]
						_tmp.append(tmp)
					elif next_char == 'i':
						bdata= self.fobj.read(4)
						consumed += 4
						tmp = struct.unpack('<i', bdata)[0]
						_tmp.append(tmp)
					elif next_char == 'I':
						bdata = self.fobj.read(4)
						consumed += 4
						tmp = struct.unpack('<I', bdata)[0]
						_tmp.append(tmp)
					elif next_char =='f':
						bdata = self.fobj.read(4)
						consumed += 4
						tmp = struct.unpack('<f', bdata)[0]
						_tmp.append(tmp)
				value = _tmp
			record.aux[tag] = value
		return record
