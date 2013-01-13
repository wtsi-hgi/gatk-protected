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

package org.broadinstitute.sting.utils.distributedutils;

import org.apache.log4j.Logger;
import org.apache.lucene.store.*;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.io.File;
import java.io.IOException;

/**
 * User: depristo
 * Date: 1/19/11
 * Time: 8:24 AM
 *
 * A reentrant lock for a shared file common file in the file system.  Relies on a a Lucene SimpleFSLock
 * to manage on disk file locking.
 */
public class SharedFileLock extends ClosableReentrantLock { // todo -- kinda gross inheritance.  The super lock is never used
    private static Logger logger = Logger.getLogger(SharedFileLock.class);

    private static final String VERIFY_HOST = System.getProperty("verify.host", "gsa1");
    private static final boolean VERIFY = false;
    private static final int VERIFY_PORT = 5050;

    // 5 minutes => 360 seconds of trying -> failure
    protected static final int DEFAULT_N_TRIES = 1000;
    protected static final long DEFAULT_MILLISECONDS_PER_TRY = 360;

    /** The file we are locking */
    private final File file;

    private final LockFactory lockFactory;
    private Lock fileLock = null;

    /**
     * A counter that indicates the number of 'locks' on this file.
     * If locks == 2, then two unlocks are required
     * before any resources are freed.
     */
    int fileLockReentrantCounter = 0;

    // type of locking
    private final int nRetries;
    private final long milliSecPerTry;

    /**
     * Create a SharedFileThreadSafeLock object locking the file
     * @param file
     */
    public SharedFileLock(File file, int nRetries, long milliSecPerTry, int ID) {
        super();
        this.file = file;
        this.nRetries = nRetries;
        this.milliSecPerTry = milliSecPerTry;

        File lockDir = new File(file.getParent() == null ? "./" : file.getParent());
        try {
            LockFactory factory = new SimpleFSLockFactory(lockDir);
            if ( VERIFY ) { // don't forget to start up the VerifyLockServer
                this.lockFactory = new VerifyingLockFactory((byte)ID, factory, VERIFY_HOST, VERIFY_PORT);
            } else {
                this.lockFactory = factory;
            }
        } catch (IOException e) {
            throw new UserException.CouldNotCreateOutputFile(lockDir, "Could not create coordination file locking directory " + lockDir, e);
        }
    }

    public SharedFileLock(File file, int ID) {
        this(file, DEFAULT_N_TRIES, DEFAULT_MILLISECONDS_PER_TRY, ID);
    }

    @Override
    public void close() {
        if ( ownsLock() ) throw new ReviewedStingException("closing SharedFileLock while still owned: ownership count " + fileLockReentrantCounter);
    }

    @Override
    public int getHoldCount() {
        return fileLockReentrantCounter;
    }

    @Override
    public boolean ownsLock() {
        return fileLockReentrantCounter > 0;
    }

    // ------------------------------------------------------------------------------------------
    //
    // workhorse routines -- acquiring file locks
    //
    // ------------------------------------------------------------------------------------------

    private boolean obtainFileLock() throws IOException {
        // annoying bug work around for verifylockserver
        if ( VERIFY )
            try {
                return fileLock.obtain(1);
            } catch ( LockObtainFailedException e ) {
                return false;
            }
        else
            return fileLock.obtain();
    }

    /**
     * Two stage [threading then file] locking mechanism.  Reenterant in that multiple lock calls will be
     * unwound appropriately.  Uses file channel lock *after* thread locking.
     */
    @Override
    public void lock() {
        if ( SharedFileThreadSafeLock.DEBUG ) logger.warn("            lock() " + Thread.currentThread().getName() + ", fileLockReentrantCounter = " + fileLockReentrantCounter);
        if ( fileLockReentrantCounter++ == 0 ) {
            // Precondition -- lock is always null while we don't have a lock
            if ( fileLock != null )
                throw new ReviewedStingException("BUG: lock() function called when a lock already is owned!");

            int i = 1;
            fileLock = lockFactory.makeLock(file.getName() + ".lock");
            try {
                boolean obtained = obtainFileLock(); // todo -- maybe use intrinsic lock features
                for ( ; ! obtained && i < nRetries; i++ ) {
                    try {
                        //logger.warn("tryLock failed on try " + i + ", waiting " + milliSecPerTry + " millseconds for retry");
                        Thread.sleep(milliSecPerTry);
                    } catch ( InterruptedException e ) {
                        throw new UserException("SharedFileThreadSafeLock interrupted during wait for file lock", e);
                    }
                    obtained = obtainFileLock(); // gross workaround for error in verify server
                }

                if ( i > 1 ) logger.warn("tryLock required " + i + " tries before completing, waited " + i * milliSecPerTry + " millseconds");

                if ( ! obtained ) {
                    fileLock = null;
                    // filelock == null -> we never managed to acquire the lock!
                    throw new UserException("SharedFileThreadSafeLock failed to obtain the lock after " + nRetries + " attempts");
                }

                if ( SharedFileThreadSafeLock.DEBUG ) logger.warn("            lock() " + Thread.currentThread().getName() + ", obtained = " + obtained + ", tries = " + i);
            } catch (IOException e) {
                fileLock = null;
                throw new ReviewedStingException("Coordination file could not be created because a lock could not be obtained.", e);
            }
        }
    }

    @Override
    public void unlock() {
        // update for reentrant unlocking
        if ( fileLock == null ) throw new ReviewedStingException("BUG: file lock is null -- file lock was not obtained");
        if ( fileLockReentrantCounter <= 0 ) throw new ReviewedStingException("BUG: file lock counter < 0");

        // this unlock counts as 1 unlock.  If this is our last unlock, actually do something
        if ( SharedFileThreadSafeLock.DEBUG ) logger.warn("            unlock() " + Thread.currentThread().getName() + ", count = " + fileLockReentrantCounter);
        if ( --fileLockReentrantCounter == 0 ) {
            try {
                if ( ! fileLock.isLocked() ) throw new ReviewedStingException("BUG: call to unlock() when we don't have a valid lock!");
                fileLock.release();
                if ( SharedFileThreadSafeLock.DEBUG ) logger.warn("            unlock() " + Thread.currentThread().getName() + ", actually releasing");
            } catch ( IOException e ) {
                throw new ReviewedStingException("Could not free file lock on file " + file, e);
            } finally {  // make sure we null out the filelock, regardless of our state
                fileLock = null;
            }
        } else {
            if ( SharedFileThreadSafeLock.DEBUG ) logger.warn("            unlock() " + Thread.currentThread().getName() + ", skipping, count = " + fileLockReentrantCounter);
        }
    }
}
