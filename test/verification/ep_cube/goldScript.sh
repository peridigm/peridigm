#!/bin/bash
# Run this script in the test directory of the build tree
# It creates gold copies of test files and movies them to the trunk directory
v=$(ls)
basename="ep_cube"
gold="gold"
for f in $v; do
  s=${f#"$basename"}
  t="$basename"_$gold$s
  if [ "$s" != "$f" ] ; then
	cp $f /home/jamitch/peridigm_projects/peridigm/test/verification/ep_cube/$t
  fi
done
