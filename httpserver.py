import SimpleHTTPServer
import SocketServer

handl = SimpleHTTPServer.SimpleHTTPRequestHandler
s = SocketServer.TCPServer((192.168.1.141,6000),handl)
s.serve_forever()
