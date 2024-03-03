// LuaFunc: simplifies use of runtime-provided lua-scripted functions.
// luafunc.cpp - implementation.

#include "luafunc.h"
#include <cstdlib>

namespace LuaFunc
{
    namespace Internal
    {
        // The Lua state used by LuaFunc.
        lua_State* L = 0;

        void Init()
        {
            if (L == 0)
            {
                L = luaL_newstate();
                luaL_openlibs(L);
                std::atexit(Exit);
            }
        }

        void Exit()
        {
            lua_close(L);
            L = 0;
        }

        RegistryHandle NewHandle()
        {
            static std::uintptr_t ch = 0;
            return (void*)(ch++);
        }
    };
};
