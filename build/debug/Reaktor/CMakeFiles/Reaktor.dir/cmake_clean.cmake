FILE(REMOVE_RECURSE
  "../lib/libReaktor-d.pdb"
  "../lib/libReaktor-d.so"
)

# Per-language clean rules from dependency scanning.
FOREACH(lang)
  INCLUDE(CMakeFiles/Reaktor.dir/cmake_clean_${lang}.cmake OPTIONAL)
ENDFOREACH(lang)
