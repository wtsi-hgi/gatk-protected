## A standalone utility class for reading python bed/bim/fam files.
## Usage:
## >>> reader = new PlinkBinaryReader("/path/to/file_no_extension")
REQ_VERSION = (3,0)
import StingConversionUtils as utils
import collections
import os
import struct
import sys

cur_version = sys.version_info
if ( cur_version < REQ_VERSION ):
 raise BaseException("Must use Python 3.0 or greater")

class HeaderMode:
 INDIVIDUAL_MAJOR,SNP_MAJOR = range(2)

class FileIterator(collections.Iterable):
 '''A class for iterating directly over a file and parsing the relevant object'''
 def __init__(self,fileHandle,returnClass):
  self.returnClass = returnClass
  self.fileHandle = fileHandle

 def __iter__(self):
  return self

 def __next__(self):
  return self.next()

 def next(self):
  line = self.fileHandle.readline()
  if ( line == "" ):
   raise StopIteration
  return self.returnClass(line)

class ReusableCachedIterator(collections.Iterable):
 '''A class for reading in file information and iterating over the cached information,looping back to the top'''
 def __init__(self,fileHandle,returnClass):
  self.returnClass = returnClass
  self.objectCache = list(map(lambda z: self.returnClass(z),fileHandle.readlines()))
  self.cacheLength = len(self.objectCache)
  self.cacheIndex = 0 

 def __iter__(self):
  return self

 def __next__(self):
  return self.next()
 def next(self):
  if ( self.cacheIndex == self.cacheLength ):
   raise StopIteration
  object = self.objectCache[self.cacheIndex]
  self.cacheIndex += 1
  return object 

 def reset(self):
  self.cacheIndex = 0

 def size(self):
  return len(self.objectCache)

class BedGenotypeIterator(collections.Iterable):
 '''A class for iterating over bytes in a bed file and returning genotypes'''
 PLINK_MAGIC_BYTE1 = b"l" 
 PLINK_MAGIC_BYTE2 = b"\x1b" 
 
 def getMode(bedFileHandle):
  byte1 = bedFileHandle.read(1)
  byte2 = bedFileHandle.read(1) 
  assert byte1 == BedGenotypeIterator.PLINK_MAGIC_BYTE1, "Byte 1 does not match PLINK magic byte: %d %d" %(int(byte1,32),int(BedGenotypeIterator.PLINK_MAGIC_BYTE1,32))
  assert byte2 == BedGenotypeIterator.PLINK_MAGIC_BYTE2, "Byte 2 does not match PLINK magic byte: %d %d" %(int(byte2,32),int(BedGenotypeIterator.PLINK_MAGIC_BYTE2,32))
  mode = struct.unpack('b',bedFileHandle.read(1))[0]
  return mode

 def __init__(self,fileHandle):
  self.fileHandle = fileHandle
  self.currentByteDecoded = self.decode(self.fileHandle.read(1))
  self.currentGenotypeOffsetInByte = 0

 def __iter__(self):
  return self

 def __next__(self):
  return self.next()

 def next(self):
  if ( self.currentByteDecoded == None ):
   raise StopIteration
  if ( self.currentGenotypeOffsetInByte < 3 ):
   genotype = self.currentByteDecoded[self.currentGenotypeOffsetInByte]
   self.currentGenotypeOffsetInByte += 1
  else:
   genotype = self.currentByteDecoded[self.currentGenotypeOffsetInByte]
   nextByte = self.fileHandle.read(1)
   #print("Next byte: "+str(nextByte))
   if ( nextByte != b"" ):
    self.currentByteDecoded = self.decode(nextByte)
    self.currentGenotypeOffsetInByte = 0
   else:
    self.currentByteDecoded = None
    self.currentOffsetInByte = None
  return genotype

 def decode(self,genoByte):
  genoInt = struct.unpack('b',genoByte)[0]
  geno1 = 3 & genoInt
  geno2 = (12 & genoInt) >> 2
  geno3 = (48 & genoInt) >> 4
  geno4 = (192 & genoInt) >> 6
  return [geno1, geno2, geno3, geno4]

# The lowest level class, one which grants access to the raw genotypes, sample, and map data
class PlinkBinaryReader:

 def __init__(self,basePath):
  self.verifyFiles(basePath)
  self.famFile = open(basePath+".fam",'r')
  self.bimFile = open(basePath+".bim",'r')
  self.bedFile = open(basePath+".bed",'rb')
  self.mode = BedGenotypeIterator.getMode(self.bedFile) 
  self.sampleIterator = ReusableCachedIterator(self.famFile,utils.Plink.Sample)
  if ( self.mode == HeaderMode.SNP_MAJOR ):
   self.snpIterator = ReusableCachedIterator(self.bimFile,utils.Plink.Variant)
   self.numMajorGenotypes = self.snpIterator.size()
   self.currentMinor = self.sampleIterator.next()
  elif ( self.mode == HeaderMode.INDIVIDUAL_MAJOR ):
   self.snpIterator = FileIterator(self.bimFile,utils.Plink.Variant)
   self.numMajorGenotypes = self.sampleIterator.size()
   self.currentMinor = self.snpIterator.next()
  self.genotypeIterator = BedGenotypeIterator(self.bedFile)
  self.majorOffset = 0

 def __iter__(self):
  return self

 def verifyFiles(self,path):
  assert os.path.exists(path+".fam")
  assert os.path.exists(path+".bim")
  assert os.path.exists(path+".bed")

 def __next__(self):
  return self.next()

 def next(self):
  if ( self.majorOffset == None ):
   raise StopIteration
  try:
   if ( self.mode == HeaderMode.SNP_MAJOR ):
    genotype = utils.Plink.Genotype(self.snpIterator.next(),self.currentMinor,self.genotypeIterator.next())
   else:
    genotype = utils.Plink.Genotype(self.currentMinor,self.sampleIterator.next(),self.genotypeIterator.next())
  except StopIteration:
   raise BaseException("Iteration stopped prematurely!")
  self.majorOffset += 1
  if ( self.majorOffset == self.numMajorGenotypes ):
    # increment the minor
    try:
     if ( self.mode == HeaderMode.SNP_MAJOR ):
      self.currentMinor = self.sampleIterator.next()
      self.snpIterator.reset()
     else:
      self.currentMinor = self.snpIterator.next()
      self.sampleIterator.reset()
     self.majorOffset = 0
    except StopIteration:
     self.majorOffset = None
  return genotype

 def __getSamples__(self):
  # store the offset of the sample iterator
  iterOffset = self.sampleIterator.cacheIndex
  # set the index to 0
  self.sampleIterator.cacheIndex = 0
  # read the samples into a list
  samples = [s for s in self.sampleIterator]
  # reset the offset
  self.sampleIterator.cacheIndex = iterOffset
  # return samples
  return samples


def getReader(base):
 return PlinkBinaryReader(base)

def getSamples(reader):
 return reader.__getSamples__()
