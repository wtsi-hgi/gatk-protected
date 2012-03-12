/*
 * Copyright (c) 2012, The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.tools;

import net.sf.picard.io.IoUtil;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.CommandLineProgram;
import org.broadinstitute.sting.utils.crypt.CryptUtils;
import org.broadinstitute.sting.utils.crypt.GATKKey;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.io.File;

/**
 * Standalone utility to generate a new GATK user key.
 *
 * By default, the newly-generated key is placed in the global GATK user key
 * directory defined in CryptUtils.GATK_USER_KEY_DIRECTORY. Existing keys
 * will not be overwritten unless -f is specified.
 *
 * The name of the new key will consist of a sanitized version of
 * the user's email address followed by a ".key" extension.
 *
 * This utility requires access to the GATK master private key,
 * and so can only be run by members of the group "gsagit", and
 * only while the Broad filesystem is accessible.
 *
 * Usage:
 *
 * ant private (if necessary)
 * java -cp dist/StingUtils.jar org.broadinstitute.sting.tools.GenerateGATKUserKey -e "user@email.address"
 */
public class GenerateGATKUserKey extends CommandLineProgram {

    private static Logger logger = Logger.getLogger(GenerateGATKUserKey.class);

    // REQUIRED ARGUMENTS:

    @Argument(fullName = "emailAddress", shortName = "e", doc = "The user email address for which to generate a GATK key", required = true)
    private String emailAddress = null;

    // OPTIONAL ARGUMENTS:

    @Argument(fullName = "outputDirectory", shortName = "o", doc = "Directory to which the GATK user key should be written (defaults to the global GATK user key directory)", required = false)
    private File outputDirectory = new File(CryptUtils.GATK_USER_KEY_DIRECTORY);

    @Argument(fullName = "forceOverwriteKeys", shortName = "f", doc = "If specified, allow existing GATK user keys to be overwritten", required = false)
    private boolean overwriteKeys = false;

    @Argument(fullName = "privateKey", shortName = "privateKey", doc = "File containing the private key with which to sign the new GATK user key (defaults to the GATK master private key)", required = false)
    private File privateKeyFile = null;

    @Argument(fullName = "publicKey", shortName = "publicKey", doc = "File containing the public key with which to validate the new GATK user key (defaults to the public key distributed with the GATK)", required = false)
    private File publicKeyFile = null;


    protected int execute() throws Exception {
        GATKKey gatkKey = new GATKKey(privateKeyFile == null ? CryptUtils.loadGATKMasterPrivateKey() : CryptUtils.readPrivateKey(privateKeyFile),
                                      publicKeyFile == null ? CryptUtils.loadGATKDistributedPublicKey() : CryptUtils.readPublicKey(publicKeyFile),
                                      emailAddress);

        // Sanitize the email address before using it as a file name:
        String keyFileName = IoUtil.makeFileNameSafe(String.format("%s.key", emailAddress));
        File keyFile = new File(outputDirectory, keyFileName);

        if ( keyFile.exists() && ! overwriteKeys ) {
            throw new UserException(String.format("File %s already exists -- refusing to overwrite unless -f is specified.",
                                                  keyFile.getAbsolutePath()));
        }

        logger.info(String.format("Writing new GATK user key to file %s", keyFile.getAbsolutePath()));
        gatkKey.writeKey(keyFile);
        logger.info(String.format("Successfully wrote %d bytes", keyFile.length()));

        return 0;
    }

    public static void main ( String[] args ) {
        try {
            GenerateGATKUserKey instance = new GenerateGATKUserKey();
            start(instance, args);
            System.exit(CommandLineProgram.result);
        } catch ( UserException e ) {
            exitSystemWithUserError(e);
        } catch ( Exception e ) {
            exitSystemWithError(e);
        }
    }
}
