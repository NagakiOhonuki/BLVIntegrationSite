#!/usr/bin/env python

import sys,subprocess
argvs = sys.argv

bam_f_name = argvs[1]
IS_f_name = argvs[2]
out_dir = argvs[3]

view_cmd = "samtools view %(bam_f_name)s" % {"bam_f_name":bam_f_name}
sam_l = subprocess.check_output(view_cmd, shell=True).strip().split("\n")
sam_d = {}
for line in sam_l:
  line = line.strip().split()
  read_name = line[0]
  if read_name not in sam_d:
    sam_d[read_name] = []
  sam_d[read_name].append(line)

IS_f = open(IS_f_name)
IS_f.next()
read_d = {}
for line in IS_f:
  line = line.strip().split()
  IS_Id = line[0]
  left_read_l = line[5].split(",")
  right_read_l = line[10].split(",")
  read_l = left_read_l + right_read_l
  if '.' not in read_l:
    out_sam_f_name = out_dir + "/" + IS_Id + ".sam"
    out_sam_f = open(out_sam_f_name, "w")
    for read_name in read_l :
      for line in sam_d[read_name]:
        print >> out_sam_f, "\t".join(line)
