from mod_python import apache

# Executable modpython server

def handler(req):
    req.content_type = "text/plain"
    req.write("Hello, World!")
    return apache.OK
