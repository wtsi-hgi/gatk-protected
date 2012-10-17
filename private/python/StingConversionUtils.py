def plinkGenotypeToVCF(plinkGenotype):
 vcfGenotype = VCF.Genotype()
 vcfGenotype.sample = plinkGenotype.sample.individual_id
 refAllele = VCF.Genotype.Allele(plinkGenotype.variant.ref,True)
 altAllele = VCF.Genotype.Allele(plinkGenotype.variant.alt,False)
 vcfGenotype.alleles,vcfGenotype.type = {
   Plink.Genotype.Type.HOM_REF : ([refAllele,refAllele],VCF.Genotype.Type.HOM_REF), 
   Plink.Genotype.Type.HET : ([refAllele,altAllele],VCF.Genotype.Type.HET),  
   Plink.Genotype.Type.HOM_VAR : ([altAllele,altAllele],VCF.Genotype.Type.HOM_VAR), 
   Plink.Genotype.Type.NO_CALL : ([VCF.Genotype.Allele.NO_CALL(),VCF.Genotype.Allele.NO_CALL()],VCF.Genotype.Type.NO_CALL)
 }[plinkGenotype.type]
 return vcfGenotype

def plinkGenotypeToVCFString(plinkGenotype):
 vcfGenotype = plinkGenotypeToVCF(plinkGenotype)
 return str(vcfGenotype)

class VCF:
 class Genotype:
  
  class Type:
   HOM_REF,HET,HOM_VAR,NO_CALL,POLYPLOID = range(5)

  class Allele:
   NO_CALL_BASES = "."

   def __init__(self,bases,isReference=False):
    self.bases = bases
    self.isReference = isReference

   def isNoCall(self):
    return self.bases == VCF.Genotype.Allele.NO_CALL_BASES

   def NO_CALL():
    return VCF.Genotype.Allele(VCF.Genotype.Allele.NO_CALL_BASES,False)

  def __init__(self):
   self.sampleID = None
   self.alleles = None
   self.likelihoods = None
   self.attributes = None
   self.type = None

  def getSample(self):
   return self.sampleID

  def getType(self):
   return self.type

  def isHomRef(self):
   return self.type == VCF.Genotype.Type.HOM_REF

  def isHomVar(self):
   return self.type == VCF.Genotype.Type.HOM_VAR

  def isHet(self):
   return self.type == VCF.Genotype.Type.HET

  def isNoCall(self):
   return self.type == VCF.Genotype.Type.NO_CALL

  def isDiploid(self):
   return len(self.alleles) == 2

  def getAlternateAlleleCount(self):
   if ( self.isDiploid() ):
    return { VCF.Genotype.Type.HOM_REF : 0,
             VCF.Genotype.Type.HOM_VAR : 2,
             VCF.Genotype.Type.HET : 1,
             VCF.Genotype.Type.NO_CALL : 0}[self.type]
   else:
    refAlelle = filter(lambda t: t.isReference,self.alleles)[0]
    return len(self.alleles) - self.alleles.count(refAllele)

  def getAlleleCount(self,allele):
   return self.alleles.count(allele)
  
  def __repr__(self):
   return "(%s,%s,%s,%s,%s)" % (repr(self.sampleID),repr(self.alleles),repr(self.attributes),repr(self.type))

  def __str__(self):
   # this assumes a bi-allelic context and only GT attributes
   if ( self.isDiploid() ):
    return { VCF.Genotype.Type.HOM_REF : "0/0",
             VCF.Genotype.Type.HOM_VAR : "1/1",
             VCF.Genotype.Type.HET : "0/1",
             VCF.Genotype.Type.NO_CALL : "./."}[self.type]
   else:
    refAllele = filter(lambda t: t.isReference,self.alleles)[0]
    refCt = self.alleles.count(refAllele)
    return "/".join(map(lambda x: "0" if x < refCt else "1",range(len(self.alleles))))

  def toString(self,alleles,format):
   if ( self.attributes == None ):
    return toString(self,alleles)
   else:
    self.attributes["GT"] = self.toString(alleles)
    return ":".join(map(lambda x: str(self.attributes[x]) if x in self.attributes else ".",format))

  def toString(self,alleles):
   if ( self.isDiploid() ):
    # find the alternate allele in the list
    altAllele = filter(lambda x: not x.isReference,self.alleles)
    altRepr = str(alleles.index(altAllele))
    return { VCF.Genotype.Type.HOM_REF : "0/"+altRepr,
             VCF.Genotype.Type.HOM_VAR : altRepr+"/"+altRepr,
             VCF.Genotype.Type.HET : "0/"+altRepr,
             VCF.Genotype.Type.NO_CALL : "./."}[self.type]
   else:
    raise NotImplementedError("Polyploid genotypes not implemented") 

class Plink:
 class Genotype:

  DOSAGE_ENC = [0,-1,1,2]
  DOSAGE_ENC_FLT = [0.,-1.,1.,2.]

  class Type:
   HOM_REF,NO_CALL,HET,HOM_VAR = range(4)

  def __init__(self,variant,sample,type):
   self.sample = sample
   self.variant = variant
   self.type = type
   self.dosage = Plink.Genotype.__getDosage__(self)

  def getDosage(self):
   return self.dosage

  def setType(self,newType):
   self.type = newType
   self.dosage = Plink.Genotype.__getDosage__(self)

  def isNoCall(genotype):
   return genotype.type == Plink.Genotype.Type.NO_CALL

  def dosageFromEnc(genotypeEnc):
   return Plink.Genotype.DOSAGE_ENC[genotypeEnc]

  def dosageFromEncAsFloat(genotypeEnc):
   return Plink.Genotype.DOSAGE_ENC_FLT[genotypeEnc]

  def __getDosage__(genotype):
   return { Plink.Genotype.Type.HOM_REF : 0,
            Plink.Genotype.Type.NO_CALL : -1,
            Plink.Genotype.Type.HET : 1,
            Plink.Genotype.Type.HOM_VAR : 2 }[genotype.type]

  def __str__(self):
   return "%s:%s:%s" % (str(self.sample),str(self.variant),str(self.type))

 ## todo -- there is no appropriate way to sort on PLINK variants unless the whole bim file is read in
 ## to establish contig order. Since we deal with human, assume standard 1-22,X,Y,MT
 class Variant:
  def __init__(self,bimLine):
   (chr,id,cm,pos,ref,alt) = bimLine.strip().split()
   self.chr = self.chr2int(chr)
   self.chrStr = chr
   self.pos = int(pos)
   self.id = id
   self.cm = cm
   self.ref = ref
   self.alt = alt
   self.frequency = None # none == uncalculated
   self.missing = None # uncalculated.

  def __repr__(self):
   return "Variant: (%d,%d)" %(self.chr,self.pos)

  def __lt__(self,other):
   if ( self.chr < other.chr ):
    return True
   elif ( self.pos < other.pos ):
    return True
   return False

  def __eq__(self,other):
   if ( other == None ):
    return False
   if ( self.chr == other.chr ):
    return self.pos == other.pos
   return False

  def __hash__(self):
   return self.chr.__hash__() ^ self.pos.__hash__()

  def distanceTo(self,other):
   if ( self.chr != other.chr ):
    return 500000000
   else:
    return abs(other.pos-self.pos)

  def chr2int(self,c):
   try:
    return int(c)
   except ValueError:
    if ( c == "X" ):
     return 23
    if ( c == "Y" ):
      return 24
    if ( c == "MT" ):
      return 25
   raise BaseError("Currently works only for contigs 1-22,X,Y,MT. Try removing the 'chr' prefix? Exceptioning Contig "+c)

 class Sample:
  def __init__(self,famLine):
   (famid,id,patid,matid,sex,pheno) = famLine.strip().split()[0:6]
   self.family_id = famid
   self.paternal_id = patid
   self.maternal_id = matid
   self.individual_id = id
   self.phenotype = pheno
   self.sex = sex

  def __repr__(self):
   return "Sample: (%s,%s)" % (self.family_id,self.individual_id)

  def __hash__(self):
   return self.individual_id.__hash__()

  def __eq__(self,other):
   return repr(self)==repr(other)

