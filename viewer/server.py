#!/usr/bin/env python3
"""
Albedo Hyperspectral Web Visualizer Server.
Serves static assets and provides lightweight JSON API endpoints for
inspecting synthetic scenes, running live J unmixing pipelines, and viewing diagnostics.
"""

import os
import sys
import http.server
import socketserver
import json
import subprocess

PORT = 8080
VIEWER_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(VIEWER_DIR)

class AlbedoRequestHandler(http.server.SimpleHTTPRequestHandler):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, directory=VIEWER_DIR, **kwargs)

    def do_GET(self):
        if self.path == "/api/status":
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.end_headers()
            status = {
                "engine": "Albedo Hyperspectral Engine",
                "version": "1.0.0",
                "core_language": "J904",
                "bands": 224,
                "range_nm": [400.0, 2500.0],
                "status": "online"
            }
            self.wfile.write(json.dumps(status).encode("utf-8"))
        elif self.path == "/api/run_pipeline":
            # Invoke J pipeline
            jconsole = os.path.join(PROJECT_ROOT, "bin", "jconsole")
            if not os.path.exists(jconsole):
                jconsole = "/Users/aneeshkrishna/.local/bin/jconsole"
            
            cmd = f"""
1!:44 < '{PROJECT_ROOT}'
load 'pipeline/pipeline_unmix.ijs'
'RadCube ReflNoisy A_gt M_true Wvls Plume'=: albedo_sim_generate_scene (8; 8; 30.0)
'BoaCube M_est Abundances dim_k rmse'=: (Wvls; 120; 35.0) albedo_pipeline_unmix_run_unmixing_pipeline RadCube
echo 'PIPELINE_COMPLETE_RMSE:', (": rmse)
exit 0
"""
            try:
                proc = subprocess.Popen([jconsole], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, cwd=PROJECT_ROOT)
                stdout, stderr = proc.communicate(input=cmd, timeout=10)
                self.send_response(200)
                self.send_header("Content-Type", "application/json")
                self.end_headers()
                res = {"success": True, "output": stdout, "error": stderr}
                self.wfile.write(json.dumps(res).encode("utf-8"))
            except Exception as e:
                self.send_response(500)
                self.send_header("Content-Type", "application/json")
                self.end_headers()
                self.wfile.write(json.dumps({"error": str(e)}).encode("utf-8"))
        else:
            super().do_GET()

def main():
    port = int(sys.argv[1]) if len(sys.argv) > 1 else PORT
    with socketserver.TCPServer(("", port), AlbedoRequestHandler) as httpd:
        print(f"Albedo Hyperspectral Visualizer serving at http://127.0.0.1:{port}")
        try:
            httpd.serve_forever()
        except KeyboardInterrupt:
            print("\nShutting down server.")
            httpd.shutdown()

if __name__ == "__main__":
    main()
