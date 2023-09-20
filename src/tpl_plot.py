#! /usr/bin/env python3
from __future__ import division, print_function

# Just a wrapper to display the template in the browser with plotly.

# Javascript only allows to fetch files from the local system by drag&drop.
# Here a local server is setup to sent a file from command line to csv_plotter.
# Furthermore, we also need to allow CORS...
# The server closes finally; thus the browser cannot reload the url/file.

# Example
# SERVALPATH/src/tpl_plot.py <filename>

from http.server import HTTPServer, SimpleHTTPRequestHandler, os, sys
from _thread import start_new_thread
from urllib.parse import quote_plus


class CORSRequestHandler(SimpleHTTPRequestHandler):
    def end_headers(self):
        self.send_header('Access-Control-Allow-Origin', '*')
        SimpleHTTPRequestHandler.end_headers(self)

def plot(tag, args=None):
    browser = 'xdg-open'   # or e.g. firefox
    qtag = quote_plus(tag)
    pwd = os.getcwd()
    filename = 'http://localhost:8000' + pwd +'/'+qtag+'/'+qtag+'.fits'
    url = 'http://localhost:8000' + os.path.dirname(os.path.realpath(__file__)) + '/tpl_plot.html'
    args = args if args else 'title='+qtag+'.fits'
    os.chdir("/")   # start the serve in root directory
    start_new_thread(HTTPServer(('', 8000), CORSRequestHandler).serve_forever, (1,))
    os.system(f"{browser} '{url}?file={filename}&{args}' && sleep 6")   # waits a bit to sent the file
    os.chdir(pwd)

if __name__ == '__main__':
    plot(sys.argv[1], sys.argv[2:])

