module dffdd.utils.msgpackrpc;

import std.string;
import std.traits;

//import msgpackrpc;


/+
private
string makeInterfaceMethods(I)(string format)
{
    auto app = appender!string();
    foreach(m; __traits(allMembers, I))
      static if(mixin(`isCallable!(I.` ~ m ~ `)`))
        app.formattedWrite(format, m);

    return app.data;
}



final class MsgpackRPCClient(Interface) : Interface
{
    this(TCPClient client = null)
    {
        _client = client;
    }


    void connect(TCPClient client)
    {
        _client = client;
    }


    mixin(makeInterfaceMethods!Interface(
    q{
        ReturnType!(Interface.%1$s) %1$s(Parameters!(Interface.%1$s) args)
        {
            static immutable name = "%1$s";

          static if(is(typeof(return) == void))
            _client.notify(name, args);
          else
            return _client.call!(typeof(return))(name, args);
        }
    }));


  private:
    TCPClient _client;
}


/*
class MsgpackRPCClientObject(Interface) : Interface
{

}
*/
+/