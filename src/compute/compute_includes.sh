# Produce compute_includes.hpp file

rm -f compute_includes.hpp
list=`ls *.hpp`
for file in $list; do
  qfile="$file"
  echo "#include \"""./compute/"$file"\"" >> compute_includes.hpp
done

