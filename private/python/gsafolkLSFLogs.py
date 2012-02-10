import os.path
import sys
from optparse import OptionParser
import gzip
import datetime, calendar
import MySQLdb

host = "calcium.broadinstitute.org"
gsafolk_lsf_table = "gsafolk_lsf"
MANY_SIZE = 100000

# Creates a SQL table in the MySQL server calcium at the Broad that contains only 
# key information about the LSF usage of members of the gsafolk fairshare group
# 
# Does this by first building a list of gsafolk uids, selecting lsf info from matter's
# table, and inserts this information into the gsafolk_lsf queue as part of the 
# GATK schema.  The standard way to run this is with incremental refreshes enabled,
# so that the program only fetches new raw lsf records with timestamps beyond the
# max timestamp present in the GATK LSF table.
#
# The default way to run this program is via cron with 'python private/python/gsafolkLSFLogs.py'
#
def main():
    global OPTIONS
    usage = "usage: %prog [options]"
    parser = OptionParser(usage=usage)
    
    parser.add_option("-v", "--verbose", dest="verbose",
                        action='store_true', default=False,
                        help="If provided, verbose progress will be enabled")         

    parser.add_option("-n", "--dry-run", dest="dryRun",
                        action='store_true', default=False,
                        help="Don't submit to the DB")

    parser.add_option("-f", "--fullRefresh", dest="fullRefresh",
                        action='store_true', default=False,
                        help="If true, we will reload the entire DB")

    parser.add_option("-M", "--maxRecords", dest="maxRecords",
                        type='int', default=None,
                        help="")

    parser.add_option("-Y", "--firstYear", dest="firstYear",
                        type='int', default=2010,
                        help="")
         
    (OPTIONS, args) = parser.parse_args()
    if len(args) != 0:
        parser.error("No arguments can be specified")

    dbGATK = DBConnector("gatk")
    dbMatter = DBConnector("matter")

    if OPTIONS.fullRefresh:
        setupTable(dbGATK)
    gsafolk_uids = select_gsafolk_uids(dbMatter)
    create_gsafolk_lsfdb(dbMatter, dbGATK, gsafolk_uids)
    
    dbGATK.close()
    dbMatter.close()
    
class DBConnector:
    """Encapsilates DB connection, allowing code to execute SQL queries to this connection"""
    def __init__(self, db):
        self.db = db
        if not OPTIONS.dryRun: 
            self.connection = MySQLdb.connect( host=host, db=self.db, user="gsamember", passwd="gsamember" )
            self.dbc = self.connection.cursor() 

    def execute(self, command):
        if OPTIONS.verbose: print "EXECUTING: ", command
        if not OPTIONS.dryRun: 
            self.dbc.execute(command)
        if OPTIONS.verbose: print '  DONE'        

    def executemany(self, command, params):
        if OPTIONS.verbose: print "EXECUTING: ", command
        if not OPTIONS.dryRun: 
            self.dbc.executemany(command, params)
        if OPTIONS.verbose: print '  DONE' 
        
    def get_summary_value(self, command):
        """Execute query and return the first records' first value"""
        self.execute(command)
        val = self.dbc.fetchall()[0][0]
        return val
    
    def close(self):
        if not OPTIONS.dryRun: 
            self.dbc.close()
            self.connection.close()
       
def select_gsafolk_uids(db):
    """Query the lsf_acct DB and return a dictionary mapping uids -> username for all users in gsafolk group"""
    q = "select U.uid, U.username from users U INNER JOIN lsf_fairshare S ON U.username = S.username WHERE S.groupname = 'gsafolk'"
    db.execute(q)
    data = db.dbc.fetchall()
    ids = dict(map(lambda x: (int(x[0]), x[1]), data))
    print 'gsafolk: ', ", ".join(ids.itervalues())
    return ids

def create_gsafolk_lsfdb(dbMatter, dbGATK, uidMap):
    """Workhorse procedure that gets records from dbMatter and inserts them in dbGATK
    
    Gets records from dbMatter that match our time constraints, formats them for dbGATK,
    and inserts them into dbGATK.  Will be incremental, if command line arguments allow this."""
    def insertRecs(records):
        asStrs = map(lambda x: map(str, x), records)
        if OPTIONS.verbose: 
            for rec in asStrs: print rec
        cmd = "INSERT INTO " + gsafolk_lsf_table + " VALUES(%s, %s, %s, %s, %s, %s)"
        dbGATK.executemany(cmd, asStrs)

    def parseRec(rec):
        #print rec
        uid, queue, timestamp, cpu, memMb = rec
        date = datetime.date.fromtimestamp(timestamp)
        if date.year >= OPTIONS.firstYear:
            username = uidMap[uid]
            cpuDays = cpu / (60*60*24.0)
            memGB = memMb/1024.0
            return username, queue, date.isoformat(), timestamp, cpuDays, memGB 
        else:
            return None

    uids = uidMap.keys()

    # count up records for printing
    maxTimeStampYear = calendar.timegm(datetime.date(OPTIONS.firstYear, 1, 1).timetuple())
    print '### timestamp from first year', maxTimeStampYear
    maxTimeStampDb = dbGATK.get_summary_value("select max(timestamp) from " + gsafolk_lsf_table)
    if maxTimeStampDb == None: maxTimeStampDb = datetime.date.min
    print '### timestamp from existing db records', maxTimeStampDb
    maxTimeStamp = max(maxTimeStampYear, maxTimeStampDb) 
    print '### Last timestamp for incremental refresh', maxTimeStamp

    qCondition = "where uid in (" + ",".join(map(str, uids)) + ") and date > " + str(maxTimeStamp)
    if OPTIONS.maxRecords != None:
        qCondition = qCondition + " LIMIT " + str(OPTIONS.maxRecords)
        
    # count up records for printing
    nRecordsToProcess = dbMatter.get_summary_value("select count(uid) from lsf_acct " + qCondition)
    print '### Number of records to be processed', nRecordsToProcess

    # select info for putting into another DB    
    q = "SELECT L.uid, Q.name, L.date, L.cpu, L.mem from lsf_acct L inner join lsf_queues Q on L.qid = Q.id " + qCondition
    dbMatter.execute(q)
    
    nTotal = 0
    nSkipped = 0
    nPassed = 0
    records = dbMatter.dbc.fetchmany(MANY_SIZE)
    while records:
        parsed = map(parseRec, records)
        passing = filter(lambda x: x != None, parsed)
        nTotal = nTotal + len(parsed)
        nPassed = nPassed + len(passing)
        nSkipped = nSkipped + len(records) - len(passing) 
        insertRecs(passing)
        print 'records', nTotal, 'passed', nPassed, 'skipped', nSkipped, '% complete', (100.0*nTotal) / nRecordsToProcess
        records = dbMatter.dbc.fetchmany(MANY_SIZE)
    

def setupTable(db):
    """Configure dbGATK gsafolk_lsf table by dropping existing one and creating it from scratch"""
    try:
        db.execute("DROP TABLE " + gsafolk_lsf_table)
    except Exception, e:
        print 'WARNING: ', e
    
    fields = "(\
        username varchar(32),\
        queue varchar(255),\
        date DATETIME,\
        timestamp BIGINT,\
        cpu_days DOUBLE,\
        mem_gb DOUBLE)"
    db.execute("CREATE TABLE " + gsafolk_lsf_table + " " + fields)

if __name__ == "__main__":
    main()
