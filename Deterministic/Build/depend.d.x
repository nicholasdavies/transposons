config.o: config.cpp config.h luafunc.h config_def.h
data.o: data.cpp data.h string.hpp
determ.o: determ.cpp config.h luafunc.h config_def.h data.h string.hpp \
  helper.h population_rna.h population_dna.h
helper.o: helper.cpp helper.h config.h luafunc.h config_def.h string.hpp
luafunc.o: luafunc.cpp luafunc.h
population_dna.o: population_dna.cpp population_dna.h helper.h config.h \
  luafunc.h config_def.h data.h string.hpp population_rna.h
population_rna.o: population_rna.cpp population_rna.h helper.h config.h \
  luafunc.h config_def.h data.h string.hpp
