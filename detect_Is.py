#!/usr/bin/env python

import sys, subprocess, re
argvs = sys.argv

bam_f_name = argvs[1]
chrom_size_f = argvs[2]
max_read_dist = argvs[3]
min_mapq = argvs[4]

#judging whether mate reads are on LTRs

read_on_5LTR_l = []
view_cmd = "samtools view -f 16 -q %(min_mapq)s %(bam_f_name)s" % {"bam_f_name":bam_f_name, "min_mapq":min_mapq}
sam_l = subprocess.check_output(view_cmd, shell=True).strip().split("\n")
for line in sam_l:
  line = line.strip().split()
  read = line[0]
  seq_name = line[2]
  if seq_name == "BLV_LTR":
    read_on_5LTR_l.append(read)

read_on_3LTR_l = []
view_cmd = "samtools view -F 16 -q %(min_mapq)s %(bam_f_name)s" % {"bam_f_name":bam_f_name, "min_mapq":min_mapq}
sam_l = subprocess.check_output(view_cmd, shell=True).strip().split("\n")
for line in sam_l:
  line = line.strip().split()
  read = line[0]
  seq_name = line[2]
  if seq_name == "BLV_LTR":
    read_on_3LTR_l.append(read)

read_on_5int_l = []
view_cmd = "samtools view -f 16 -q %(min_mapq)s %(bam_f_name)s" % {"bam_f_name":bam_f_name, "min_mapq":min_mapq}
sam_l = subprocess.check_output(view_cmd, shell=True).strip().split("\n")
for line in sam_l:
  line = line.strip().split()
  read = line[0]
  seq_name = line[2]
  if seq_name == "BLV_noLTR":
    read_on_5int_l.append(read)

read_on_3int_l = []
view_cmd = "samtools view -F 16 -q %(min_mapq)s %(bam_f_name)s" % {"bam_f_name":bam_f_name, "min_mapq":min_mapq}
sam_l = subprocess.check_output(view_cmd, shell=True).strip().split("\n")
for line in sam_l:
  line = line.strip().split()
  read = line[0]
  seq_name = line[2]
  if seq_name == "BLV_noLTR":
    read_on_3int_l.append(read)


#generating intermediate files
left_cmd = "samtools view -q %(min_mapq)s -b -F 16 %(bam_f_name)s | bedtools merge -i \'stdin\' -d %(max_read_dist)s -c 1 -o count,collapse | grep -v 'LTR' > %(bam_f_name)s.left.bed" % {"bam_f_name":bam_f_name, "max_read_dist":max_read_dist, "min_mapq":min_mapq}

right_cmd = "samtools view -q %(min_mapq)s -b -f 16 %(bam_f_name)s | bedtools merge -i \'stdin\' -d %(max_read_dist)s -c 1 -o count,collapse | grep -v 'LTR' > %(bam_f_name)s.right.bed" % {"bam_f_name":bam_f_name, "max_read_dist":max_read_dist, "min_mapq":min_mapq}

both_cmd = "bedtools slop -i %(bam_f_name)s.left.bed -g %(chrom_size_f)s -l 0 -r %(max_read_dist)s | bedtools intersect -a  \'stdin\' -b %(bam_f_name)s.right.bed -wa -wb | awk \'BEGIN{OFS=\"\t\"}{print $1,$2,$3-%(max_read_dist)s,$4,$5,$6,$7,$8,$9,$10}\' > %(bam_f_name)s.both.bed" % {"bam_f_name":bam_f_name, "chrom_size_f":chrom_size_f, "max_read_dist":max_read_dist}

subprocess.check_output(left_cmd, shell=True)
subprocess.check_output(right_cmd, shell=True)
subprocess.check_output(both_cmd, shell=True)

left_f_name = "%(bam_f_name)s.left.bed" % {"bam_f_name":bam_f_name}
right_f_name = "%(bam_f_name)s.right.bed" % {"bam_f_name":bam_f_name}
both_f_name = "%(bam_f_name)s.both.bed" % {"bam_f_name":bam_f_name}
both_f = open(both_f_name)
left_f = open(left_f_name)
right_f = open(right_f_name)

left_l = []
right_l = []
res_l = []

for line in both_f:
  line = line.strip().split()
  left = line[:5]
  right = line[5:]
  left_l.append(left)
  right_l.append(right)
  total = int(line[3]) + int(line[8])
  line.append(total)
  res_l.append(line)

for line in left_f:
  line = line.strip().split()
  if line not in left_l:
    for i in range(5):
      line.append('.')
    line.append(line[3])
    res_l.append(line)

for line in right_f:
  line = line.strip().split()
  if line not in right_l:
    for i in range(5):
      line.insert(0,'.')
    line.append(line[8])
    res_l.append(line)

#output
out_f_name = "%(bam_f_name)s.out.txt" % {"bam_f_name":bam_f_name}
out_f = open(out_f_name,"w")
print >> out_f, "IS_Id\tChrom.left\tStart.left\tEnd.left\tReads.left\tReads_list.left\tChrom.right\tStart.right\tEnd.right\tReads.right\tReads_list.right\tReads.total\tStrand\tReads_on_5LTR\tReads_on_3LTR\tReads_on_5int\tReads_on_3int"

i = 0
both_d = {}
for line in res_l:
  IS_Id = "IS_" + str(i)
  left_num_on_5LTR = 0
  left_num_on_3LTR = 0
  left_num_on_5int = 0
  left_num_on_3int = 0

  if line[4] != '.':
    read_l = [re.sub(r'\/[0-9]+','',read) for read in line[4].split(",")]
    line[4] = ",".join(read_l)
    for read in read_l:
      if read in read_on_5LTR_l:
        left_num_on_5LTR += 1
      if read in read_on_3LTR_l:
        left_num_on_3LTR += 1
      if read in read_on_5int_l:
        left_num_on_5int += 1
      if read in read_on_3int_l:
        left_num_on_3int += 1

  right_num_on_5LTR = 0
  right_num_on_3LTR = 0
  right_num_on_5int = 0
  right_num_on_3int = 0
  if line[9] != '.':
    read_l = [re.sub(r'\/[0-9]+','',read) for read in line[9].split(",")]
    line[9] = ",".join(read_l)
    for read in read_l:
      if read in read_on_5LTR_l:
        right_num_on_5LTR += 1
      if read in read_on_3LTR_l:
        right_num_on_3LTR += 1
      if read in read_on_5int_l:
        right_num_on_5int += 1
      if read in read_on_3int_l:
        right_num_on_3int += 1
  num_on_5LTR = left_num_on_5LTR + right_num_on_5LTR
  num_on_3LTR = left_num_on_3LTR + right_num_on_3LTR
  num_on_5int = left_num_on_5int + right_num_on_5int
  num_on_3int = left_num_on_3int + right_num_on_3int

  num_on_plus = left_num_on_5LTR + left_num_on_5int + right_num_on_3LTR + right_num_on_3int
  num_on_minus = left_num_on_3LTR + left_num_on_3int + right_num_on_5LTR + right_num_on_5int
  strand = "nd"
  if num_on_plus > num_on_minus * 5:
    strand = "+"
  elif num_on_minus > num_on_plus * 5:
    strand = "-"
  line = line + [strand,num_on_5LTR,num_on_3LTR,num_on_5int,num_on_3int]
  print >> out_f, IS_Id + "\t" + "\t".join([str(c) for c in line])
  i += 1

#removing intermediate files
rm_cmd = "rm %(bam_f_name)s.*.bed" % {"bam_f_name":bam_f_name}
subprocess.check_output(rm_cmd, shell=True)

