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

import org.apache.log4j.Logger;
import org.broadinstitute.sting.commandline.Advanced;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.CommandLineProgram;
import org.broadinstitute.sting.utils.crypt.CryptUtils;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.io.File;
import java.security.*;

/**
 * Standalone utility to generate a new public/private key pair.
 *
 * By default, the new key pair is placed in the current directory
 * with the following naming scheme:
 *
 * GATK_private.key
 * GATK_public.key
 *
 * The utility will refuse to overwrite existing keys unless -f is
 * specified.
 *
 * It is not recommended to override the default encryption settings
 * (key length, algorithm, etc.) unless you know exactly what you're
 * doing.
 *
 * On very rare occasions we may want to use this utility to replace the
 * existing GATK master public/private key pair with a new one. Doing so
 * will revoke all existing GATK user keys from the current version onwards,
 * so this is a fairly drastic action! The procedure for doing this is
 * documented in the GSA internal wiki.
 *
 * Usage:
 *
 * To generate a new public/private key pair in the current directory with
 * the default encryption settings:
 *
 * ant private (if necessary)
 * java -cp dist/StingUtils.jar org.broadinstitute.sting.tools.GenerateKeyPair
 *
 * To generate a new public/private key pair in the directory /local/mykeys instead:
 *
 * ant private (if necessary)
 * java -cp dist/StingUtils.jar org.broadinstitute.sting.tools.GenerateKeyPair -o /local/mykeys
 */
public class GenerateKeyPair extends CommandLineProgram {

    private static Logger logger = Logger.getLogger(GenerateGATKUserKey.class);

    // All arguments are optional -- the defaults should work fine in most cases:

    @Argument(fullName = "outputDirectory", shortName = "o", doc = "Directory to which the key pair should be written (defaults to current working directory)", required = false)
    private File outputDirectory = new File(".");

    @Argument(fullName = "keyNamePrefix", shortName = "prefix", doc = "The public/private key names will start with this prefix", required = false)
    private String keyNamePrefix = DEFAULT_KEY_NAME_PREFIX;

    @Argument(fullName = "forceOverwriteKeys", shortName = "f", doc = "If specified, allow existing keys to be overwritten", required = false)
    private boolean overwriteKeys = false;

    @Advanced
    @Argument(fullName = "keyLength", shortName = "keyLength", doc = "Length of the key in bits (DO NOT SPECIFY UNLESS YOU KNOW EXACTLY WHAT YOU'RE DOING)", required = false)
    private int keyLength = CryptUtils.DEFAULT_KEY_LENGTH;

    @Advanced
    @Argument(fullName = "encryptionAlgorithm", shortName = "encryptionAlgorithm ", doc = "Encryption algorithm to use to generate the key pair (DO NOT SPECIFY UNLESS YOU KNOW EXACTLY WHAT YOU'RE DOING)", required = false)
    private String encryptionAlgorithm = CryptUtils.DEFAULT_ENCRYPTION_ALGORITHM;

    @Advanced
    @Argument(fullName = "randNumberAlgorithm", shortName = "randNumberAlgorithm", doc = "Random number generation algorithm to use to generate the key pair (DO NOT SPECIFY UNLESS YOU KNOW EXACTLY WHAT YOU'RE DOING)", required = false)
    private String randNumberAlgorithm = CryptUtils.DEFAULT_RANDOM_NUMBER_GENERATION_ALGORITHM;

    public static final String DEFAULT_KEY_NAME_PREFIX = "GATK";
    public static final String PRIVATE_KEY_SUFFIX = "_private.key";
    public static final String PUBLIC_KEY_SUFFIX = "_public.key";


    protected int execute() throws Exception {
        File privateKeyFile = new File(outputDirectory, keyNamePrefix + PRIVATE_KEY_SUFFIX);
        File publicKeyFile = new File(outputDirectory, keyNamePrefix + PUBLIC_KEY_SUFFIX);

        checkForPreexistingKeys(privateKeyFile, publicKeyFile);

        logger.info("Generating public/private key pair");
        logger.info(String.format("Key Length: %d bits", keyLength));
        logger.info(String.format("Encryption Algorithm: %s", encryptionAlgorithm));
        logger.info(String.format("Random Number Generation Algorithm: %s", randNumberAlgorithm));
        logger.info("...");

        KeyPair keyPair = CryptUtils.generateKeyPair(keyLength, encryptionAlgorithm, randNumberAlgorithm);
        CryptUtils.writeKeyPair(keyPair, privateKeyFile, publicKeyFile);

        logger.info("Successfully wrote key pair to files:");
        logger.info(String.format("Private Key: %s", privateKeyFile));
        logger.info(String.format("Public Key: %s", publicKeyFile));

        return 0;
    }

    private void checkForPreexistingKeys ( File ... keyFiles ) {
        for ( File keyFile : keyFiles ) {
            if ( keyFile.exists() && ! overwriteKeys ) {
                throw new UserException(String.format("Key file %s already exists -- refusing to overwrite unless -f is specified",
                                                      keyFile.getAbsolutePath()));
            }
        }
    }

    public static void main ( String[] args ) {
        try {
            GenerateKeyPair instance = new GenerateKeyPair();
            start(instance, args);
            System.exit(CommandLineProgram.result);
        } catch ( UserException e ) {
            exitSystemWithUserError(e);
        } catch ( Exception e ) {
            exitSystemWithError(e);
        }
    }
}
