// LuaFunc: simplifies use of runtime-provided lua-scripted functions.
// luafunc.h - API
// (C) 2013 Nicholas G. Davies

// Example usage:
//
// #include <iostream>
// #include "luafunc.h"
// using namespace LuaFunc;
//
// int main()
// {
//     LuaFunction<double> L_sq("x", "x*x");        // function returns double: args [1] and expression [2]
//     LuaFunction<double> L_pow("x,y -> x^y");     // function returns double: compact notation, arrow separating args & expression
//
//     std::cout << L_sq(2) << "\n";     // prints 4
//     std::cout << L_pow(2, 3) << "\n"; // prints 8
// }

#ifndef LUAFUNC_H
#define LUAFUNC_H

#include <memory>
#include <lua/lua.hpp>
#include <string>
#include <iostream>

namespace LuaFunc
{
    namespace Internal
    {
        // The Lua state used by LuaFunc.
        extern lua_State* L;
        void Init();
        void Exit();

        typedef void* RegistryHandle;
        RegistryHandle NewHandle();
        inline void pushhandle(lua_State* L, RegistryHandle h) { lua_pushlightuserdata(L, h); }

        // Push: pushes one or several values onto the Lua stack.
        inline void Push(int v)         { lua_pushinteger(L, v); }
        inline void Push(bool v)        { lua_pushboolean(L, v); }
        inline void Push(double v)      { lua_pushnumber(L, v); }
        inline void Push(const char* v) { lua_pushstring(L, v); }
        inline void Push(std::string v) { lua_pushstring(L, v.c_str()); }

        template <typename T, typename... U>
        void Push(T next, U... remaining)
        {
            Push(next);
            Push(remaining...);
        }

        // To: converts a value on the Lua stack into a C++ type.
        template <typename T>
        inline T To(int index) { (void) index; T().CompilationError("Unregistered Type"); }

        template<> inline int To<int>(int index)                  { return lua_tointeger(L, index); }
        template<> inline bool To<bool>(int index)                { return lua_toboolean(L, index); }
        template<> inline double To<double>(int index)            { return lua_tonumber(L, index); }
        template<> inline const char* To<const char*>(int index)  { return lua_tostring(L, index); }
        template<> inline std::string To<std::string>(int index)  { return lua_tostring(L, index); }

        // RegistryEntry: an entry on the Lua registry, which cleans itself up upon destruction.
        // First, push the desired value onto the Lua stack, then call Register().
        class RegistryEntry
        {
        public:
            // Create a registry entry with the given unique handle h.
            RegistryEntry(RegistryHandle h) : handle(h) { }

            // Delete the registry entry.
            ~RegistryEntry()
            {
                if (Internal::L)
                {
                    pushhandle(Internal::L, handle);                // push registry handle to stack
                    lua_pushnil(Internal::L);                       // push nil to stack
                    lua_settable(Internal::L, LUA_REGISTRYINDEX);   // set registry entry to nil
                }
            }

            // Register the value at the top of the stack. Pops the value from the stack.
            void Register()
            {
                pushhandle(Internal::L, handle);                // push registry handle to stack
                lua_insert(Internal::L, -2);                    // move registry handle below value on stack
                lua_settable(Internal::L, LUA_REGISTRYINDEX);   // set value in registry
            }

            // Put the registered value at the top of the stack.
            void Get()
            {
                pushhandle(Internal::L, handle);                // push handle to stack
                lua_gettable(Internal::L, LUA_REGISTRYINDEX);   // get entry on stack
            }

        private:
            RegistryHandle handle;
        };
    };

    // LuaFunction<ReturnValue>: wraps a Lua function.
    template <typename ReturnValue>
    class LuaFunction
    {
    public:
        // default constructor
        LuaFunction() : n_args_pushed(0)
        { }

        // args of the form "x,y"; expr of the form "x^y+0.5"
        LuaFunction(std::string args, std::string expr)
         : compact_form(args + "->" + expr), n_args_pushed(0)
        {
            Create("return function(" + args + ") return " + expr + " end");
        }

        // compact_notation of the form x,y->x^y+0.5 or of the form 3^7+9 (no arguments, no arrow)
        LuaFunction(std::string compact_notation)
         : compact_form(compact_notation), n_args_pushed(0)
        {
            std::string::size_type arrow = compact_notation.find("->");
            if (arrow == std::string::npos)
            {
                Create("return function() return " + compact_notation + " end");
            }
            else
            {
                Create("return function(" + compact_notation.substr(0, arrow) + ") return " + compact_notation.substr(arrow + 2) + " end");
            }
        }

        // Call function with provided arguments
        template <typename... T>
        ReturnValue operator()(T... args) const
        {
            lua_checkstack(Internal::L, sizeof...(args) + 1);   // ensure there is space enough on the stack
            registry_entry->Get();                              // get Lua function on stack
            Internal::Push(args...);                            // get arguments on stack
            lua_call(Internal::L, sizeof...(args), 1);          // call function, enforce one return value
            ReturnValue r = Internal::To<ReturnValue>(-1);      // get the returned value
            lua_pop(Internal::L, 1);                            // pop the returned value from the stack
            return r;
        }

        // Call function with given vector as arguments
        template <typename ForwardIterator>
        ReturnValue Call(ForwardIterator first, ForwardIterator last) const
        {
            lua_checkstack(Internal::L, last - first + 1);      // ensure there is space enough on the stack
            registry_entry->Get();                              // get Lua function on stack
            for (ForwardIterator i = first; i != last; ++i)     // get arguments on stack
                Internal::Push(*i);
            lua_call(Internal::L, last - first, 1);             // call function, enforce one return value
            ReturnValue r = Internal::To<ReturnValue>(-1);      // get the returned value
            lua_pop(Internal::L, 1);                            // pop the returned value from the stack
            return r;
        }

        // Push argument for later use of Call()
        template <typename Arg>
        LuaFunction<ReturnValue>& operator << (Arg a)
        {
            Internal::Push(a);
            ++n_args_pushed;
            return *this;
        }

        // Call function with previously-pushed arguments
        ReturnValue Call()
        {
            lua_checkstack(Internal::L, n_args_pushed + 1);     // ensure there is space enough on the stack
            registry_entry->Get();                              // get Lua function on stack
            lua_insert(Internal::L, -n_args_pushed - 1);        // move Lua function to correct stack position (before arguments)
            lua_call(Internal::L, n_args_pushed, 1);            // call function, enforce one return value
            n_args_pushed = 0;                                  // reset argument counter
            ReturnValue r = Internal::To<ReturnValue>(-1);      // get the returned value
            lua_pop(Internal::L, 1);                            // pop the returned value from the stack
            return r;
        }

        // Provide compact form of function definition
        std::string CompactForm() const
        {
            return compact_form;
        }

    private:
        std::unique_ptr<Internal::RegistryEntry> registry_entry;
        std::string compact_form;
        int n_args_pushed;

        void Create(std::string command)
        {
            Internal::Init();                                           // initialise Lua, if needed
            registry_entry.reset(new Internal::RegistryEntry(Internal::NewHandle()));   // create registry entry
            luaL_dostring(Internal::L, command.c_str());                // push function to stack
            registry_entry->Register();                                 // store function in registry
        }
    };
};

#endif  // LUAFUNC_H