function get_gcf(str) {
  return gensub(/(.*)GCF_([0-9]+)\.([0-9]+)(.*)/, "GCF_\\2.\\3", "g", str)
}
# GCF = correct GCF, results = probabilitys from decoded result
function output(GCF, results, stats) {
  top1_id=0 # get max
  top2_id=0 # get 2nd best max

  correct_hit=0
  ambiguous_hit=0
  for (id in results)
  {
    if (id == GCF && results[id] >= 0.9)
      correct_hit=1
    else if (id != GCF && results[id] >= 0.05)
      ambiguous_hit=1

    if (top1_id == 0 || results[id] > results[top1_id]) {
      top2_id = top1_id
      top1_id = id
    }
    else if (top2_id == 0 || results[id] > results[top2_id]) {
      top2_id = id
    }
  }
  
  if (correct_hit && !ambiguous_hit) {
    #printf "%s:\t%s\n", GCF, "Very Certain"
    stats["VeryCertain"]++
  }
  else if (GCF == top1_id && (top2_id == 0 || results[top1_id] - results[top2_id] >= 0.25)) {
    #printf "%s:\t%s\n", GCF, "Certain"
    stats["Certain"]++
  }
  else if (top1_id == 0 || (top2_id != 0 && results[top1_id] - results[top2_id] < 0.03)) {
    printf "%s:\t%s\n", GCF, "Ambiguous/Unknown"
    stats["Ambiguous"]++
  }
  else if (top1_id == GCF) {
    printf "%s:\t%s\n", GCF, "Uncertain"
    stats["Uncertain"]++
  }
  else if (top1_id != GCF) {
    printf "%s:\t%s\n", GCF, "Wrong"
    stats["Wrong"]++
  }
  stats["Total"]++

  if (!((correct_hit && !ambiguous_hit) || (GCF == top1_id && (top2_id == 0 || results[top1_id] - results[top2_id] > 0.25))))
    for (id in results)
      print id "\t" results[id]
}
BEGIN{
  stats["Total"]=0
  stats["VeryCertain"]=0
  stats["Certain"]=0
  stats["Uncertain"]=0
  stats["Ambiguous"]=0
  stats["Wrong"]=0
}
{
  if (last_filename == FILENAME) {
    if (!($0 ~ /^#/)) {
      $1 = get_gcf($1)
      results[$1] = $2
    }
  }
  else {
    if (NR!=1) {
      output(GCF, results, stats);
      delete results
    }
    GCF=get_gcf(FILENAME)
    last_filename = FILENAME
  }
}
END{
  output(GCF, results, stats);
  print "Total:\t" stats["Total"]
  print "Very Certain:\t" stats["VeryCertain"]
  print "Certain:\t" stats["Certain"]
  print "Uncertain:\t" stats["Uncertain"]
  print "Ambiguous:\t" stats["Ambiguous"]
  print "Wrong:\t" stats["Wrong"]
}

