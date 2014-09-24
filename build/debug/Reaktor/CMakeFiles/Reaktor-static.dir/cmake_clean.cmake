FILE(REMOVE_RECURSE
  "../lib/libReaktor-static-d.pdb"
  "../lib/libReaktor-static-d.a"
)

# Per-language clean rules from dependency scanning.
FOREACH(lang)
  INCLUDE(CMakeFiles/Reaktor-static.dir/cmake_clean_${lang}.cmake OPTIONAL)
ENDFOREACH(lang)
