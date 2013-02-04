#!/bin/tcsh

# documentation on using this code is https://iwww.broadinstitute.org/gsa/wiki/index.php/GATK_logs_and_Tableau

# download CLI tools
# http://aws.amazon.com/developertools/AWS-Identity-and-Access-Management/4143

# where does Java live?
setenv JAVA_HOME /Library/Java/Home

# the path to the IAM command line interface tools
setenv AWS_IAM_HOME ~/Downloads/IAMCli-1.5.0
setenv PATH $AWS_IAM_HOME/bin:$PATH

# See http://docs.aws.amazon.com/IAM/latest/CLIReference/Setup.html#d0e258
# This document contains the AWS access and secret keys for the owner of the GATK AWS account (currently Mark DePristo)
# The only way to get this information is to talk to Mark
setenv AWS_CREDENTIAL_FILE /Users/depristo/aws-account-key

# Should we create these users?
setenv CREATE_USERS false

# WARNING -- running this can invalid the GATK out in the wild.  Use with caution!
setenv UPDATE_USER_KEYS false

# Update the AWS IAM security policies for each user
setenv UPDATE_USER_POLICY true


set GATK_LOG_UPLOADER_USER = "GATK2LogUploader"
set GATK_LOG_DOWNLOADER_USER = "GATK2LogDownloader"

foreach user ($GATK_LOG_UPLOADER_USER $GATK_LOG_DOWNLOADER_USER) 
    echo "#####################################################"
    echo "Sorting out user $user"
    echo "#####################################################"

    # Create the GATK user
    # update the secret key
    if $CREATE_USERS == true then
	echo "Creating AWS user $user"
	iam-usercreate -u $user -k -v > ${user}_cred.txt
    endif

    # the user access and secret keys are stored encrypted in the GATK
    # and must be updated to be the most current ones
    if $UPDATE_USER_KEYS == true then
	echo "Updating user keys"
	set keys = `iam-userlistkeys -u $user | grep -v "Active" | grep -v "IsTruncated"`
	foreach key ( $keys )
	    echo "Deleting key $key for user $user"
	    iam-userdelkey -u $user -k $key
	end
	echo "Adding new AWS key for user $user"
	iam-useraddkey -u $user > ${user}_cred.txt 
	cat ${user}_cred.txt
    endif

    if $UPDATE_USER_POLICY == true then
	echo "Updating user policies"
	set policies = `iam-userlistpolicies -u $user | grep -v "IsTruncated"`
	foreach policy ( $policies )
	    echo "Deleting policy $policy for user $user"
	    iam-userdelpolicy -u $user -p $policy
	end
	iam-useruploadpolicy -u $user -p ${user}Policy -f ${user}Policy.txt
    endif

    # List the user policies for $user
    echo "User policies for $user"
    iam-userlistpolicies -u $user -v 

    echo "#####################################################"
    echo "Done with user $user"
    echo "#####################################################"
end

