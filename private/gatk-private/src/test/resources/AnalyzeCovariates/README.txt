
*.table Files
===========

  The input recalibration report files. The name of the file indicates the role of each, so 'bqsr.table' is the one to be provided using the -BQSR argument and so forth.

after-*.table
===========

  Alternative after.table that present some parameter changes with respect to bqsr.table and before.table files. These after-*.table were generated from 
  a smaller interval (chr20:10000000-11000000) for convenience, but that should not make a difference as far as testing is concern. 

*-IBA.* Files
=============

  These files are not necessary for the test to run but are included so you can check differences manually between these and the output of the corresponding tests.

  IBA stands for I = Input (-BQSR), B = (-before), A = (-after). So basically IBA means that all three files were provided in the command line. 
  There would be *-BA.* *-IA.* *-BA.* files if other on/off combinations where to be tested for the content of their output files. 
  This is not the case in this round of test; perhaps in the future. 

  + expected-IBA.pdf -- The plots file
  + expected-IBA.csv.txt -- The intermediary file (gzipped to save space).


