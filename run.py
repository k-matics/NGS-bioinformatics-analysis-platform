"""Development entry point:  python run.py

Production-style start (web + worker + redis): scripts/start_services.sh
"""
import os

from app import create_app

app = create_app()

if __name__ == "__main__":
    # Bind to localhost only; reach it from Windows via an SSH tunnel.
    app.run(host="127.0.0.1", port=int(os.environ.get("NGS_PORT", 5000)), threaded=True, debug=False)
