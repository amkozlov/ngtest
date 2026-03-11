#! /usr/bin/env python3

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

tree_files = ["bestTree", "support", "supportTBE", "supportFBP", "supportSH", "supportRBS", "supportEBG", "supportPS", "supportPBS", 
              "consensusTreeSTRICT", "consensusTreeMR", "consensusTreeMR80", "consensusTreeMRE", "ancestralTree", "mutationMapTree"]

tsv_files = ["mutationMapList", "ancestralStates"]

loglh_commands = ["evaluate", "sitelh", "search", "all", "fast" ]

def raxml_file(prefix, suffix):
    return ".".join([prefix, "raxml", suffix])

def cmd_outfiles(command, opts={}):
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
        files = ["log", "startTree", "mlTrees", "bestTree", "bestModel", "bootstraps"]
    elif command=="bootstrap":
        files = ["log", "bootstraps"]
    elif command=="support":
        files = ["log"]
    elif command=="bsconverge":
        files = ["log"]
    elif command=="consense":
        files = ["log", "consensusTreeSTRICT", "consensusTreeMR", "consensusTreeMR80", "consensusTreeMRE"]
    elif command=="terrace":
        files = ["log"]
    elif command=="ebg":
        files = ["log", "startTree", "bestTree", "bestModel", "supportEBG"]
    elif command=="moose":
        files = ["log", "startTree", "moose.bestModel", "moose.json"]
    elif command=="ancestral":
        files = ["log", "ancestralTree", "ancestralStates", "ancestralProbs"]
    elif command=="mutmap":
        files = ["log", "mutationMapTree", "mutationMapList"]


    if (command=="all" or command=="support") and "bs-metric" in opts:
      metrics = opts["bs-metric"].split("+")
      if len(metrics) < 2 and metrics[0] == "FBP":
        files.append("support")
      else:
        for m in metrics:
          files.append(f"support{m}")

    return files

def test_must_fail(test_name):
  test_toks = test_name.split("_")
  return len(test_toks) > 2 and test_toks[2].startswith("BAD")

def parse_logfile(fname):
    if not os.path.isfile(fname):
      return None
    d = {}
    d["errors"] = set([])
#    max_tokens = 2
    with open(fname, "r") as f:
      for line in f:
        if line.startswith("ERROR:"):
          errmsg = line[7:]
          d["errors"].add(errmsg)
        elif 'Final LogLikelihood:' in line:
          d["likelihood"] = float(line.split()[2])
        elif 'Elapsed time:' in line:
          d["time"] = float(line.split()[2])
#        if len(d) >= max_tokens:
#           break
    return d 

def check_files(command, prefix, goldprefix, opts={}):
    for f in cmd_outfiles(command, opts):
        if not os.path.isfile(raxml_file(prefix, f)):
            return False
    
    return True

def check_loglh(command, d1, d2):
#    print(d1,d2)
    lh_eps = 0.1
    if abs(d1["likelihood"] - d2["likelihood"]) < lh_eps:
      return True
    else:
      print(d1["likelihood"], " <-- ", d2["likelihood"])
      return False

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

    print(dist_rf, dist_bs)

    return dist_rf == 0

def diff_comp(fname1, fname2):
  call_str = ["diff", fname1, fname2]
  retval = call(call_str)
  return True if retval == 0 else False
  

def check_tree(command, prefix, goldprefix):
    for suffix in tree_files:
      if suffix in cmd_outfiles(command):
        tree1_fname=raxml_file(prefix, suffix)
        tree2_fname=raxml_file(goldprefix, suffix)
        if os.path.isfile(tree2_fname):
            if not tree_comp(tree1_fname, tree2_fname):
                return False
        else:
           print("WARNING: ground truth not found: ", tree2_fname)

    return True

def check_tsv(command, prefix, goldprefix):
    for suffix in tsv_files:
      if suffix in cmd_outfiles(command):
        tsv1_fname=raxml_file(prefix, suffix)
        tsv2_fname=raxml_file(goldprefix, suffix)
        if os.path.isfile(tsv2_fname):
            if not diff_comp(tsv1_fname, tsv2_fname):
                return False
        else:
           print("WARNING: ground truth not found: ", tsv2_fname)

    return True
    
def check(test_name, prefix, goldprefix):
    passed = True
    test_toks = test_name.split("_")
    command = test_toks[0]
    opts = {}
    if command == "support" and len(test_toks) > 2:  
      opts["bs-metric"] = test_toks[2].upper()

    suffix = "log"
    log_fname = raxml_file(prefix, suffix)
    gold_log_fname = raxml_file(goldprefix, suffix)
    log_attrs = parse_logfile(log_fname)
    gold_log_attrs = parse_logfile(gold_log_fname)

    if test_must_fail(test_name):
      if log_attrs["errors"] != gold_log_attrs["errors"]:
        print(log_attrs, "\n", gold_log_attrs)
        return False
    elif log_attrs and len(log_attrs["errors"]) > 0:
      return False

    if not check_files(command, prefix, goldprefix, opts):
        return False
    if not check_tree(command, prefix, goldprefix):
        passed = False
    if not check_tsv(command, prefix, goldprefix):
        passed = False
    if command in loglh_commands: 
      if not check_loglh(command, log_attrs, gold_log_attrs):
          passed = False
    return passed

def raxng_ver(raxng):
    rxpipe = Popen([raxng, "-v"], stdout=PIPE)
    pat = r'RAxML-NG v. \K[[:digit:]]+\.[[:digit:]]+\.[[:digit:]]+[\-]*[[:alpha:]]*'
#    pat = 'RAxML-NG v. \K[[:digit:]]+\.[[:digit:]]+[\.]*[[:digit:]]*[\-]*[[:alpha:]]*'
    ver = check_output(["grep", "-oP", pat], stdin=rxpipe.stdout, universal_newlines=True)
#    ver = check_output(raxng + " -v")
#    print ver
#    sys.exit()
#    ver="2.0-beta2"

    return ver.strip()

if __name__ == "__main__":

  if len(sys.argv) == 1:
    print("Usage: runtest.py <raxml-ng-binary> [testname | all] [threads/workers]") 
    sys.exit() 

  raxng=sys.argv[1]
#  ver="0.7.0git"

  if (not os.path.isfile(raxng)):
     print("File not found: ", raxng)
     sys.exit()

  ver = raxng_ver(raxng)
  if (not ver):
      print("Failed to get RAxML-NG version!")
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
  print(rootdir, "out", ver, par)
  outdir=os.path.join(rootdir, "out", ver, par)
  golddir=os.path.join(rootdir, "out", "gold")
  testdir=os.path.join(rootdir, "test")
  log_fname=os.path.join(rootdir, "-".join(["log_raxng",ver]))
  if os.path.isfile(log_fname):
    os.remove(log_fname)

  if not os.path.isdir(outdir):
      os.makedirs(outdir)

  if len(sys.argv) > 2 and sys.argv[2] != "all":
      test_mask=os.path.join(testdir, sys.argv[2] + ".sh")
  else:
      test_mask=os.path.join(testdir, "*.sh")

  script_list = sorted(glob.glob(test_mask))

  print("Testing RAxML-NG v.", ver, ", tests found: ", len(script_list))

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
      expected_retval = 1 if test_must_fail(test_name) else 0
#      print(retval)
      if retval == expected_retval and check(test_name, prefix, goldprefix):
          print("OK")
      else:
          print("ERROR")
          errors += 1

  print("")
  if errors > 0:
      print(" 🔴  Tests failed: ", errors, " / ", total)
  else:
      print(" 🟢  All test completed successfully: ", total)
  print("")

# RAXNG=~/hits/raxml-ng/bin/raxml-ng-static DATADIR=~/hits/ngtest/data PREFIX=~/hits/ngtest/out/0.7.0/search_GTR_default/test


