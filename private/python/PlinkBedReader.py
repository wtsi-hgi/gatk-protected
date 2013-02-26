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
  self.nextLine = fileHandle.readline()
  self.nextLineCache = None

 def __iter__(self):
  return self

 def __next__(self):
  return self.next()

 def __peek__(self):
  return self.peek()

 def peek(self):
  if ( self.nextLine == "" ):
   return None
  if ( self.nextLineCache == None ):
   self.nextLineCache = self.returnClass(self.nextLine)
  return self.nextLineCache

 def isDone(self):
  return self.nextLine == ""

 def next(self):
  if ( self.nextLine == "" ):
   raise StopIteration
  if ( self.nextLineCache == None ):
   curObj = self.returnClass(self.nextLine)
  else:
   curObj = self.nextLineCache
  self.nextLine = self.fileHandle.readline()
  self.nextLineCache = None
  return curObj 

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

 def __peek__(self):
  return self.peek()

 def peek(self):
  if ( self.cacheIndex == self.cacheLength ):
   return None
  return self.objectCache[self.cacheIndex]

 def next(self):
  if ( self.cacheIndex == self.cacheLength ):
   raise StopIteration
  object = self.objectCache[self.cacheIndex]
  self.cacheIndex += 1
  return object 

 def isDone(self):
  return self.cacheIndex == self.cacheLength

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

 def __init__(self,fileHandle,nGenotypesPerMajor,opt=None):
  self.fileHandle = fileHandle
  self.nBytesPerMajor = int((3+nGenotypesPerMajor)/4)
  self.nGenotypesPerMajor = nGenotypesPerMajor
  if ( opt == None ):
   self.currentByteDecoded = self.decode(self.fileHandle.read(1))
   self.currentGenotypeOffsetInByte = 0
   self.nBytesPerMajor = int((3+nGenotypesPerMajor)/4)
   self.genotypesInFinalByte = nGenotypesPerMajor % 4
   if ( self.genotypesInFinalByte == 0 ):
    self.genotypesInFinalByte = 4
   self.currentByteOfMajor = 1
  if ( opt != None ):
   import itertools
   self.join = lambda z: list(itertools.chain.from_iterable(z)) 

 def __iter__(self):
  return self

 def __next__(self):
  return self.next()

 def __peek__(self):
  return self.peek()

 def peek(self):
  if ( self.currentByteDecoded == None ):
   return None
  return self.currentByteDecoded[self.currentGenotypeOffsetInByte]

 def nextMajor(self):
  # should not be called with opt = None in init, this will cause the decoded byte to go all wonky
  return self.join(map(self.decode2,self.fileHandle.read(self.nBytesPerMajor)))[0:self.nGenotypesPerMajor]

 def next(self):
  if ( self.currentByteDecoded == None ):
   raise StopIteration
  if ( self.currentByteOfMajor == self.nBytesPerMajor ):
   nGenotypesInByte = self.genotypesInFinalByte
  else:
   nGenotypesInByte = 4
  if ( self.currentGenotypeOffsetInByte < nGenotypesInByte - 1 ):
   genotype = self.currentByteDecoded[self.currentGenotypeOffsetInByte]
   self.currentGenotypeOffsetInByte += 1
  else:
   genotype = self.currentByteDecoded[self.currentGenotypeOffsetInByte]
   nextByte = self.fileHandle.read(1)
   if ( nextByte != b"" ):
    self.currentByteDecoded = self.decode(nextByte)
    self.currentGenotypeOffsetInByte = 0
   else:
    self.currentByteDecoded = None
    self.currentOffsetInByte = None
   if ( self.currentByteOfMajor == self.nBytesPerMajor ):
    self.currentByteOfMajor = 1
   else:
    self.currentByteOfMajor += 1
  return genotype

 def decode2(self,genoByte):
   geno1 = 3 & genoByte
   geno2 = (12 & genoByte ) >> 2
   geno3 = (48 & genoByte ) >> 4
   geno4 = (192 & genoByte ) >> 6
   return [geno1, geno2, geno3, geno4]

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
  if ( self.mode == HeaderMode.INDIVIDUAL_MAJOR ):
   self.majorIterator = ReusableCachedIterator(self.famFile,utils.Plink.Sample)
   self.minorIterator = ReusableCachedIterator(self.bimFile,utils.Plink.Variant)
  elif ( self.mode == HeaderMode.SNP_MAJOR ):
   self.majorIterator = FileIterator(self.bimFile,utils.Plink.Variant)
   self.minorIterator = ReusableCachedIterator(self.famFile,utils.Plink.Sample)
  self.numGenotypesPerMajor = self.minorIterator.size()
  self.genotypeIterator = BedGenotypeIterator(self.bedFile,self.numGenotypesPerMajor)
  self.majorOffset = 0

 def __iter__(self):
  return self

 def verifyFiles(self,path):
  assert os.path.exists(path+".fam")
  assert os.path.exists(path+".bim")
  assert os.path.exists(path+".bed")

 def __next__(self):
  return self.next()

 def __peek__(self):
  return self.peek()

 def peek(self):
  if ( self.majorOffset == None or self.majorIterator.peek() == None ):
   genotype = None
  elif ( self.mode == HeaderMode.SNP_MAJOR ):
   genotype = utils.Plink.Genotype(self.majorIterator.peek(),self.minorIterator.peek(),self.genotypeIterator.peek())
  else:
   genotype = utils.Plink.Genotype(self.minorIterator.peek(),self.majorIterator.peek(),self.genotypeIterator.peek())
  return genotype

 def next(self):
  if ( self.majorOffset == None or self.majorIterator.isDone() ):
   raise StopIteration
  try:
   if ( self.mode == HeaderMode.SNP_MAJOR ):
    genotype = utils.Plink.Genotype(self.majorIterator.peek(),self.minorIterator.next(),self.genotypeIterator.next())
   else:
    genotype = utils.Plink.Genotype(self.minorIterator.next(),self.majorIterator.peek(),self.genotypeIterator.next())
  except StopIteration:
   raise BaseException("Iteration stopped prematurely!")
  self.majorOffset += 1
  if ( self.majorOffset == self.numGenotypesPerMajor ):
    # done with the minor genotypes, so need to increment the major
    try:
     self.majorIterator.next()
     self.minorIterator.reset()
     self.majorOffset = 0
    except StopIteration:
     raise StopIteration 
  return genotype

 def snpMajor(self):
  return self.mode == HeaderMode.SNP_MAJOR

 def individualMajor(self):
  return self.mode == HeaderMode.INDIVIDUAL_MAJOR
 
 def __getSamples__(self):
  # store the offset of the sample iterator
  if ( self.mode == HeaderMode.SNP_MAJOR ):
   sampleIterator = self.minorIterator
  else:
   raise BaseError("Not implemented: getSamples for individual-major mode")
  iterOffset = sampleIterator.cacheIndex
  # set the index to 0
  sampleIterator.cacheIndex = 0
  # read the samples into a list
  samples = [s for s in sampleIterator]
  # reset the offset
  sampleIterator.cacheIndex = iterOffset
  # return samples
  return samples

## optimized reader for SNP-major binary files. Ultimately to replace reference version.
class SiteOptimizedPlinkBinaryReader(collections.Iterable):
 
 def __iter__(self):
  return self

 def __next__(self):
  return self.next

 def __init__(self,basePath):
  self.samples = list(map(lambda x: utils.Plink.Sample(x),open(basePath+".fam",'r').readlines()))
  self.sites = map(lambda x: utils.Plink.Variant(x),open(basePath+".bim",'r'))
  binaryReader = open(basePath+".bed",'rb')
  m = BedGenotypeIterator.getMode(binaryReader)
  if ( m != 1 ):
    raise BaseException("Optimized reader cannot be used in individual major mode.")
  self.genotypeReader = BedGenotypeIterator(binaryReader,len(self.samples),True) # turn on the optimization for nextMajor
  self.mode = HeaderMode.SNP_MAJOR
  self.numGenotypesPerMajor = len(self.samples)
  self.offset = 0

 def __getSamples__(self):
  return [x for x in self.samples] # clone it 

 def next(self):
  # return the standard ugly thing: utils.Plink.Genotypes for everything
  self.offset += 1
  return map(lambda x: utils.Plink.Genotype(self.sites.__next__(),x[0],x[1]), zip(self.samples,self.genotypeReader.nextMajor()))

 def nextDosage(self):
  # optimized: only return the variant and the dosage
  self.offset += 1
  return (self.sites.__next__(),list(map(lambda x: utils.Plink.Genotype.dosageFromEnc(x),self.genotypeReader.nextMajor())))

 def nextDosageFloat(self):
  # optimized: only return variant and dosage, with dosage as a float
  self.offset += 1
  return (self.sites.__next__(),list(map(utils.Plink.Genotype.dosageFromEncAsFloat,self.genotypeReader.nextMajor())))

 def snpMajor(self):
  return True

def getReader(base):
 return PlinkBinaryReader(base)

def getSamples(reader):
 return reader.__getSamples__()
