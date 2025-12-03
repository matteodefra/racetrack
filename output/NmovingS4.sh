for n in {5..20}; do
  for L in {50..50}; do
     ./fulltest 20 $n $L
  done
done
