{   
  str=$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6
  for (i = 7; i <= NF; i+=4)
  {
    if ($i$(i+1) != $(i+2)$(i+3))
      str=str"\t0 0\t0 0";
    else
      str= str"\t"$i" "$(i+1)"\t"$(i+2)" "$(i+3);
  }
  print str
}