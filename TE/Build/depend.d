./Build/bits.o: Source/bits.cpp Source/bits.h Source/shared.h Source/config.h \
  Source/luafunc.h Source/config_def.h Source/randomizer.h \
  Source/recorder.h Source/data.h Source/string.hpp
./Build/config.o: Source/config.cpp Source/config.h Source/luafunc.h \
  Source/config_def.h
./Build/data.o: Source/data.cpp Source/data.h Source/string.hpp
./Build/experiment.o: Source/experiment.cpp Source/recorder.h Source/data.h \
  Source/string.hpp Source/experiment.h Source/bits.h Source/shared.h \
  Source/config.h Source/luafunc.h Source/config_def.h \
  Source/randomizer.h
./Build/luafunc.o: Source/luafunc.cpp Source/luafunc.h
./Build/main.o: Source/main.cpp Source/shared.h Source/config.h Source/luafunc.h \
  Source/config_def.h Source/randomizer.h Source/recorder.h \
  Source/data.h Source/string.hpp Source/experiment.h Source/bits.h
./Build/randomizer.o: Source/randomizer.cpp Source/randomizer.h
./Build/recorder.o: Source/recorder.cpp Source/recorder.h Source/data.h \
  Source/string.hpp Source/bits.h Source/shared.h Source/config.h \
  Source/luafunc.h Source/config_def.h Source/randomizer.h
./Build/shared.o: Source/shared.cpp Source/shared.h Source/config.h \
  Source/luafunc.h Source/config_def.h Source/randomizer.h \
  Source/recorder.h Source/data.h Source/string.hpp
