bits.o: Source/bits.cpp Source/bits.h Source/shared.h Source/config.h \
  Source/luafunc.h Source/config_def.h Source/randomizer.h \
  Source/recorder.h Source/data.h Source/string.hpp
config.o: Source/config.cpp Source/config.h Source/luafunc.h \
  Source/config_def.h
data.o: Source/data.cpp Source/data.h Source/string.hpp
experiment.o: Source/experiment.cpp Source/recorder.h Source/data.h \
  Source/string.hpp Source/experiment.h Source/bits.h Source/shared.h \
  Source/config.h Source/luafunc.h Source/config_def.h \
  Source/randomizer.h
luafunc.o: Source/luafunc.cpp Source/luafunc.h
main.o: Source/main.cpp Source/shared.h Source/config.h Source/luafunc.h \
  Source/config_def.h Source/randomizer.h Source/recorder.h \
  Source/data.h Source/string.hpp Source/experiment.h Source/bits.h
randomizer.o: Source/randomizer.cpp Source/randomizer.h
recorder.o: Source/recorder.cpp Source/recorder.h Source/data.h \
  Source/string.hpp Source/bits.h Source/shared.h Source/config.h \
  Source/luafunc.h Source/config_def.h Source/randomizer.h
shared.o: Source/shared.cpp Source/shared.h Source/config.h \
  Source/luafunc.h Source/config_def.h Source/randomizer.h \
  Source/recorder.h Source/data.h Source/string.hpp
