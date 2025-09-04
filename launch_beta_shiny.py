import subprocess
import threading
import time
import webview

def run_shiny():
    # Run your Shiny app in the background
    subprocess.run(["shiny", "run", "--port", "8000", "/Users/angusdixon/masters/met591_diss/all_scripts/scripts_wip/deltaplot_wip/python_shiny_gui.py"])

# Start Shiny in a separate thread
thread = threading.Thread(target=run_shiny, daemon=True)
thread.start()

# Give Shiny time to start
time.sleep(2)

# Open app in embedded window
webview.create_window("Delta Plot GUI", "http://127.0.0.1:8000", width=1200, height=800)
webview.start()
