#! /usr/bin/env python

import os
import sys
import glob
import shutil
import datetime
import re
from subprocess import Popen,call,check_output,STDOUT,PIPE
import dendropy
from dendropy.calculate import treecompare

ver=""

tree_files = ["bestTree", "support", "consensusTreeSTRICT", "consensusTreeMR", 
              "consensusTreeMR80", "consensusTreeMRE"]

def raxml_file(prefix, suffix):
    return ".".join([prefix, "raxml", suffix])

def cmd_outfiles(command):
    files = []
    if command=="check":
        files = ["log"]
    elif command=="parse":
        files = ["log", "rba"]
    elif command=="start":
        files = ["log", "startTree"]
    elif command=="evaluate":
        files = ["log", "startTree", "bestTree", "bestModel"]
    elif command=="sitelh":
        files = ["log", "startTree", "bestTree", "bestModel", "siteLH"]
    elif command=="search":
        files = ["log", "startTree", "bestTree", "bestModel"]
    elif command=="all":
        files = ["log", "startTree", "mlTrees", "bestTree", "bestModel", "bootstraps", "support"]
    elif command=="bootstrap":
        files = ["log", "bootstraps"]
    elif command=="support":
        files = ["log", "support"]
    elif command=="bsconverge":
        files = ["log"]
    elif command=="consense":
        files = ["log", "consensusTreeSTRICT", "consensusTreeMR", "consensusTreeMR80", "consensusTreeMRE"]
    elif command=="terrace":
        files = ["log"]

    return files

def check_files(command, prefix, goldprefix):
    for f in cmd_outfiles(command):
        if not os.path.isfile(raxml_file(prefix, f)):
            return False
    
    return True

def check_loglh(command, prefix, goldprefix):
    return True

def tree_comp(tree1_fname, tree2_fname):
    tree_list = dendropy.TreeList()
    tree_list.read(path=tree1_fname, schema="newick")
    tree_list.read(path=tree2_fname, schema="newick")
    tree1=tree_list[0]
    tree2=tree_list[1]

    dist_rf = treecompare.symmetric_difference(tree1, tree2)
#    dist_eu = treecompare.euclidean_distance(tree1, tree2)
    if tree1.length() > 0 and tree1.length() > 0:
      dist_bs = treecompare.robinson_foulds_distance(tree1, tree2)
    else:
      dist_bs = 0

    print dist_rf, dist_bs

    return dist_rf == 0

def check_tree(command, prefix, goldprefix):
    for suffix in tree_files:
      if suffix in cmd_outfiles(command):
        tree1_fname=raxml_file(prefix, suffix)
        tree2_fname=raxml_file(goldprefix, suffix)
        if os.path.isfile(tree2_fname):
            if not tree_comp(tree1_fname, tree2_fname):
                return False
        else:
           print "WARNING: ground truth not found: ", tree2_fname

    return True
    
def check(test_name, prefix, goldprefix):
    command = test_name.split("_")[0]
    if not check_files(command, prefix, goldprefix):
        return False
    if not check_loglh(command, prefix, goldprefix):
        return False
    if not check_tree(command, prefix, goldprefix):
        return False
    return True

def raxng_ver(raxng):
    rxpipe = Popen([raxng, "-v"], stdout=PIPE)
    pat = 'RAxML-NG v. \K[[:digit:]]+\.[[:digit:]]+\.[[:digit:]]+[\-]*[[:alpha:]]*'
#    pat = 'RAxML-NG v. \K[[:digit:]]+\.[[:digit:]]+[\.]*[[:digit:]]*[\-]*[[:alpha:]]*'
    ver = check_output(["grep", "-oP", pat], stdin=rxpipe.stdout)
#    ver = check_output(raxng + " -v")
#    print ver
#    sys.exit()
#    ver="1.0.3-master"

    return ver.strip()

if __name__ == "__main__":

  if len(sys.argv) == 1:
    print "Usage: runtest.py <raxml-ng-binary> [testname | all] [threads/workers]" 
    sys.exit() 

  raxng=sys.argv[1]
#  ver="0.7.0git"

  if (not os.path.isfile(raxng)):
     print "File not found: ", raxng
     sys.exit()

  ver = raxng_ver(raxng)
  if (not ver):
      print "Failed to get RAxML-NG version!"
      sys.exit()

  threads=1
  workers=1

  if len(sys.argv) > 3:
    s = sys.argv[3].split("/")
    threads=int(s[0])
    if len(s) > 1:
      workers=int(s[1])

  par="T"+str(threads)+"W"+str(workers)

  rootdir=os.path.dirname(os.path.abspath(sys.argv[0]))
  datadir=os.path.join(rootdir, "data")
  outdir=os.path.join(rootdir, "out", ver, par)
  golddir=os.path.join(rootdir, "out", "gold")
  testdir=os.path.join(rootdir, "test")
  log_fname=os.path.join(rootdir, "-".join(["log_raxng",ver]))
  if os.path.isfile(log_fname):
    os.remove(log_fname)

  if not os.path.isdir(outdir):
      os.mkdir(outdir)

  if len(sys.argv) > 2 and sys.argv[2] != "all":
      test_mask=os.path.join(testdir, sys.argv[2] + ".sh")
  else:
      test_mask=os.path.join(testdir, "*.sh")

  script_list = sorted(glob.glob(test_mask))

  print "Testing RAxML-NG v.", ver, ", tests found: ", len(script_list)

#  print test_mask
  errors = 0
  total = 0
  for test_script in script_list:
#      print test_script
      test_name=os.path.basename(test_script).replace(".sh", "")
      testoutdir=os.path.join(outdir, test_name)
      testgolddir=os.path.join(golddir, test_name)
      prefix=os.path.join(testoutdir, "test")
      goldprefix=os.path.join(testgolddir, "test")
   
#      print test_name
#      sys.stdout.write(test_name + "...")

      if os.path.isdir(testoutdir):
        shutil.rmtree(testoutdir)
      os.mkdir(testoutdir)
      my_env = os.environ.copy()
      my_env["RAXNG"] = raxng
      my_env["DATADIR"] = datadir
      my_env["PREFIX"] = prefix
      my_env["THREADS"] = str(threads)
      my_env["WORKERS"] = str(workers)
      if workers > 1:
        my_env["NGARGS"] = "--workers " + str(workers)

#      my_env["NGARGS"] = "--tip-inner on"
      call_str = ["bash", test_script]
#      print call_str
#      call(call_str, env=my_env)
      retval = 0
      with open(log_fname, "a") as fout:
        retval = call(call_str, env=my_env, stdout=fout, stderr=STDOUT)

      total += 1
      sys.stdout.write(test_name + "...")
      if retval == 0 and check(test_name, prefix, goldprefix):
          print "OK"
      else:
          print "ERROR"
          errors += 1

  if errors > 0:
      print "Tests failed: ", errors, " / ", total
  else:
      print "All test completed successfully: ", total

# RAXNG=~/hits/raxml-ng/bin/raxml-ng-static DATADIR=~/hits/ngtest/data PREFIX=~/hits/ngtest/out/0.7.0/search_GTR_default/test


