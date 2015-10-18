SET(phreeqc_DATABASE
Amm.dat
frezchem.dat
iso.dat
llnl.dat
minteq.dat
minteq.v4.dat
phreeqc.dat
pitzer.dat
sit.dat
wateq4f.dat
)

IF(WIN32)
  install (FILES ${phreeqc_DATABASE} DESTINATION database)
ELSE()
  install (FILES ${phreeqc_DATABASE} DESTINATION ${CMAKE_INSTALL_DOCDIR}/database)
ENDIF()

