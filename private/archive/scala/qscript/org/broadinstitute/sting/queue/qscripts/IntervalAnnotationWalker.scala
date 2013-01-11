/*
*  By downloading the PROGRAM you agree to the following terms of use:
*  
*  BROAD INSTITUTE - SOFTWARE LICENSE AGREEMENT - FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
*  
*  This Agreement is made between the Broad Institute, Inc. with a principal address at 7 Cambridge Center, Cambridge, MA 02142 (BROAD) and the LICENSEE and is effective at the date the downloading is completed (EFFECTIVE DATE).
*  
*  WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
*  WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
*  NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
*  
*  1. DEFINITIONS
*  1.1 PROGRAM shall mean copyright in the object code and source code known as GATK2 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute/GATK on the EFFECTIVE DATE.
*  
*  2. LICENSE
*  2.1   Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM. 
*  The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only.  For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
*  2.2  No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD.  LICENSEE shall ensure that all of its users agree to the terms of this Agreement.  LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
*  2.3  License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.  
*  
*  3. OWNERSHIP OF INTELLECTUAL PROPERTY 
*  LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies.  LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
*  Copyright 2012 Broad Institute, Inc.
*  Notice of attribution:  The GATK2 program was made available through the generosity of Medical and Population Genetics program at the Broad Institute, Inc.
*  LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
*  
*  4. INDEMNIFICATION
*  LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, (Indemnitees), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
*  
*  5. NO REPRESENTATIONS OR WARRANTIES
*  THE PROGRAM IS DELIVERED AS IS.  BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
*  IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
*  
*  6. ASSIGNMENT
*  This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
*  
*  7. MISCELLANEOUS
*  7.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
*  7.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
*  7.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
*  7.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested.  All notices under this Agreement shall be deemed effective upon receipt. 
*  7.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter. 
*  7.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
*  7.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
*/

package org.broadinstitute.sting.queue.qscripts.archive

import java.io.PrintStream
import org.broadinstitute.sting.gatk.contexts.{ReferenceContext, AlignmentContext}
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker
import collection.JavaConversions._
import collection.immutable.List
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature
import org.broadinstitute.sting.gatk.refdata.features.table.TableFeature
import org.broadinstitute.sting.gatk.refdata.features.refseq.RefSeqFeature
import org.broadinstitute.sting.gatk.walkers.{TreeReducible, RefWalker}
import org.broadinstitute.sting.commandline.{Output, Argument}
import org.broadinstitute.sting.utils.{GenomeLoc}
import collection.mutable.{ListBuffer, HashSet}
import java.lang.Math
import org.broadinstitute.variant.utils.BaseUtils

class IntervalAnnotationWalker extends RefWalker[AnnotationMetaData,List[IntervalInfoBuilder]] {
  @Argument(doc="Min proportion of bases overlapping between an interval of interest and an annotation interval for annotation to occur",shortName="mpb")
  var minPropOverlap : Double = 0.50 // default to 50%
  @Output var out : PrintStream = _

  override def reduceInit : List[IntervalInfoBuilder] = {
    Nil
  }

  override def map(tracker: RefMetaDataTracker, ref: ReferenceContext, context: AlignmentContext) : AnnotationMetaData = {
    val ivalsOfInterest : List[GATKFeature] = tracker.getGATKFeatureMetaData("interval_list",false).toList
    val refSeqData : List[Object] = tracker.getReferenceMetaData("refseq").toList
    val tableMetaData : List[Object] = tracker.getReferenceMetaData("table",false).toList
    var amd : AnnotationMetaData = new AnnotationMetaData(ref.getLocus,ref.getBase)
    if ( refSeqData != null ) {
      amd.refSeqFeatures = refSeqData.map( (o: Object) => o.asInstanceOf[RefSeqFeature])
    }
    if ( tableMetaData != null ) {
      amd.tableFeatures = tableMetaData.map( (o: Object) => o.asInstanceOf[TableFeature]).map( (u: TableFeature) =>
      new Tuple2(u.getLocation,u.getAllValues.toList.reduceLeft(_ + "," + _)))
    }
    amd.newIvals = ivalsOfInterest.map(_.getLocation)
    return amd
  }

  override def reduce(mapMetaData : AnnotationMetaData, rList : List[IntervalInfoBuilder]) : List[IntervalInfoBuilder] = {
    rList.filter( (u : IntervalInfoBuilder) => u.location.isBefore(mapMetaData.locus)).foreach(u => out.print("%s%n".format(u.done)))
    val newList : List[IntervalInfoBuilder] = rList.filter( ! _.finalized ) :::
      mapMetaData.newIvals.filter( (g: GenomeLoc) => g.getStart == mapMetaData.locus.getStart ).map( u => new IntervalInfoBuilder(u,minPropOverlap) )
    newList.foreach(_.update(mapMetaData))
    return newList
  }

  override def onTraversalDone(finalList : List[IntervalInfoBuilder]) = {
    finalList.foreach(u => out.print("%s%n".format(u.done)))
  }

}

class AnnotationMetaData(loc: GenomeLoc, ref: Byte) {
  val locus : GenomeLoc = loc
  var newIvals : List[GenomeLoc] = Nil
  var refSeqFeatures : List[RefSeqFeature] = Nil
  var tableFeatures : List[(GenomeLoc,String)] = Nil
  val refBase : Byte = ref
}

class IntervalInfoBuilder(loc : GenomeLoc, minProp : Double) {

  val OVERLAP_PROP : Double = minProp

  var metaData : HashSet[String] = new HashSet[String]
  var seenTranscripts : HashSet[String] = new HashSet[String]
  var baseContent : ListBuffer[Byte] = new ListBuffer[Byte]
  var location : GenomeLoc = loc
  var gcContent : Double = _
  var entropy : Double = _
  var finalized : Boolean = false


  def update(mmd : AnnotationMetaData) = {
    baseContent += mmd.refBase
    mmd.refSeqFeatures.filter( p => ! seenTranscripts.contains(p.getTranscriptUniqueGeneName)).view.foreach(u => addRefSeq(u))
    mmd.tableFeatures.filter( (t : (GenomeLoc,String)) =>
      ! metaData.contains(t._2) && (t._1.intersect(location).size.asInstanceOf[Double]/location.size() >= OVERLAP_PROP) ).foreach( t => metaData.add(t._2))
  }

  def addRefSeq( rs: RefSeqFeature ) : Unit = {
    seenTranscripts.add(rs.getTranscriptUniqueGeneName)
    val exons : List[(GenomeLoc,Int)] = rs.getExons.toList.zipWithIndex.filter( p => ((p._1.intersect(location).size+0.0)/location.size) >= minProp )
    if ( exons.size > 0 ) {
      metaData.add("%s[%s]".format(rs.getTranscriptUniqueGeneName,exons.map( p => "exon%d".format(p._2) ).reduceLeft(_ + "," + _)))
    } else {
      if ( location.isBefore(rs.getExons.get(0)) || location.isPast(rs.getExons.get(rs.getNumExons-1)) ) {
        metaData.add("%s[%s]".format(rs.getTranscriptUniqueGeneName,"UTR"))
      } else {
        metaData.add("%s[%s]".format(rs.getTranscriptUniqueGeneName,"Intron"))
      }
    }

  }

  def done : String = {
    finalized = true
    def isGC(b : Byte) : Int = if ( BaseUtils.gIndex == b || BaseUtils.cIndex == b ) { 1 } else { 0 }
    gcContent = baseContent.foldLeft[Int](0)( (a,b) => a + isGC(b)).asInstanceOf[Double]/location.size()
    entropy = calcEntropy(baseContent.map(b => ListBuffer(b))) + calcEntropy(baseContent.reverse.map(b => ListBuffer(b)))
    val meta : String = metaData.reduceLeft(_ + "\t" + _)
    return "%s\t%d\t%d\t%.2f\t%.2f\t%s".format(location.getContig,location.getStart,location.getStop,gcContent,entropy,meta)
  }

  def calcEntropy(byteList : ListBuffer[ListBuffer[Byte]]) : Double = {
    if(byteList.size == 1) return 0
    Math.log(1+byteList.tail.size-byteList.tail.dropWhile( u => u.equals(byteList(1))).size) +
      calcEntropy(byteList.tail.foldLeft(ListBuffer(byteList(0)))( (a,b) => {
        if ( b.equals(byteList(1)) ) {
          a.dropRight(1) :+ (a.last ++ b)
        } else {
          a :+ b
        }
      }))
  }

}