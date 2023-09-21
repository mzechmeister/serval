#! /usr/bin/env python3
from __future__ import division, print_function

# Just a wrapper to display the template in the browser with plotly.

# Javascript only allows to fetch files from the local system by drag&drop.
# Here a local server is setup to sent a file from command line to csv_plotter.
# Furthermore, we also need to allow CORS...
# The server closes finally; thus the browser cannot reload the url/file.

# Example
# SERVALPATH/src/tpl_plot.py <tag>
# srv -tpl <tag>

from http.server import HTTPServer, SimpleHTTPRequestHandler, os, sys
from _thread import start_new_thread
from urllib.parse import quote_plus
import webbrowser
import time

class RequestHandler(SimpleHTTPRequestHandler):
    def do_GET(self):
        if self.path.startswith('/?'):
            self.send_response(200)
            self.send_header('Content-type', 'text/html; charset=utf-8')
            self.end_headers()
            with open(os.path.dirname(os.path.realpath(__file__))+'/tpl_plot.html', 'rb') as file:
                self.wfile.write(file.read()) # Read the file and send the contents
        else:
            SimpleHTTPRequestHandler.do_GET(self)

def plot(tag, args=None):
    qtag = quote_plus(tag)
    filename = qtag+'/'+qtag+'.fits'
    url = 'http://localhost:8000'
    args = args if args else 'title='+qtag+'.fits'
    start_new_thread(HTTPServer(('', 8000), RequestHandler).serve_forever, (1,))   # run in background
    webbrowser.open(f"{url}?file={filename}&{args}")
    time.sleep(3)   # waits a bit to sent the file

if __name__ == '__main__':
    plot(sys.argv[1], sys.argv[2:])

