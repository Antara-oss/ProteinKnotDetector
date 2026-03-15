# Protein Knot Detector Server Deployment Guide

This guide explains:

1. How this project works on a server today
2. How to deploy it on a Linux server
3. What will break or become risky in production
4. What to change when you want to scale

## 1. What This Project Is Right Now

The current repository is a Django web app, not a separated frontend/backend system.

Important parts:

- `knot_web/urls.py` exposes a single route: `/`
- `detector/views.py` handles the form submit and runs the protein analysis
- `detector/templates/detector/index.html` renders the input form, results, and 3D viewer
- `Dockerfile` starts Gunicorn on port `7860`
- `knot_web/settings.py` uses SQLite and currently has development settings enabled

## 2. How It Will Work After Deployment

When deployed in its current form, the request flow is:

1. A user opens your website.
2. Nginx forwards the request to Gunicorn.
3. Gunicorn serves the Django app.
4. Django returns the HTML page from `detector/templates/detector/index.html`.
5. The browser loads `3Dmol.js` from the public CDN.
6. The user pastes a protein sequence and submits the form.
7. Django receives the `POST` request in `detector/views.py`.
8. If the sequence length is `<= 400`, the app sends the sequence to the ESMFold API and writes the returned PDB to a local file.
9. If the sequence length is `> 400`, the app splits it into overlapping windows and repeats the same process for each fragment.
10. `topoly` runs on the generated PDB file and estimates knot type and knot probability.
11. Django renders a new HTML response with the analysis report.
12. For short sequences, the PDB text is embedded into the page and shown in the 3D viewer.

In simple terms:

- The server does the computation inside the web request
- The app depends on outbound internet access to `https://api.esmatlas.com/foldSequence/v1/pdb/`
- The app writes temporary `.pdb` files to local disk
- The browser gets a fully rendered HTML page back, not JSON

## 3. Current Production Risks

This code can run on a server, but it is only suitable for a small, low-traffic deployment in its current form.

### Security and configuration risks

- `knot_web/settings.py` has `DEBUG = True`
- `knot_web/settings.py` has `ALLOWED_HOSTS = ['*']`
- `knot_web/settings.py` contains a hard-coded Django `SECRET_KEY`
- SQLite is being used as the production database

### Runtime risks

- The Docker image runs `gunicorn` with `--workers 1`, so only one request can be processed at a time reliably
- Long computations block the request until the response is finished
- Large proteins can take a long time because ESMFold and `topoly` are both part of the same HTTP request
- Sliding-window analysis writes files like `fragment_201_600.pdb`, which can collide across users
- Fragment files are not cleaned up in the large-sequence path, so disk usage will grow over time
- The template depends on a public CDN for `3Dmol.js`, so internet access from the browser matters too

### What this means operationally

For a class project, demo, or low-traffic private deployment, the current architecture is acceptable.

For public traffic, the main bottlenecks are:

- request timeout risk
- one-analysis-at-a-time behavior
- temporary file collisions
- poor scale-out behavior across multiple containers

## 4. Recommended Deployment Model Right Now

For the current codebase, the most practical deployment is:

- one Linux server
- Docker for packaging
- Gunicorn inside the container
- Nginx as a reverse proxy
- HTTPS with Certbot

This matches the repository as it exists today and avoids extra application rewrites.

## 5. Server Requirements

Minimum for a small demo deployment:

- Ubuntu 22.04 or 24.04
- 2 vCPU
- 4 GB RAM
- 20+ GB disk
- A domain name pointed to the server IP
- Outbound internet access from the server to the ESMFold API

If you expect repeated long analyses, use:

- 4 vCPU
- 8 GB RAM

## 6. Step-by-Step Deployment on Ubuntu

### Step 1: Prepare the server

```bash
sudo apt update
sudo apt install -y git docker.io docker-compose-plugin nginx certbot python3-certbot-nginx
sudo systemctl enable --now docker
sudo systemctl enable --now nginx
```

### Step 2: Clone the project

```bash
git clone <your-repo-url>
cd ProteinKnotDetector
```

### Step 3: Build the Docker image

```bash
docker build -t protein-knot-detector .
```

### Step 4: Run the container

Bind it only to localhost so Nginx is the public entry point:

```bash
docker run -d \
  --name protein-knot-detector \
  --restart unless-stopped \
  -p 127.0.0.1:7860:7860 \
  protein-knot-detector
```

### Step 5: Configure Nginx

Create `/etc/nginx/sites-available/protein-knot-detector`:

```nginx
server {
    listen 80;
    server_name your-domain.com;

    client_max_body_size 10m;

    location / {
        proxy_pass http://127.0.0.1:7860;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
        proxy_read_timeout 180;
        proxy_connect_timeout 30;
        proxy_send_timeout 180;
    }
}
```

Enable the site:

```bash
sudo ln -s /etc/nginx/sites-available/protein-knot-detector /etc/nginx/sites-enabled/protein-knot-detector
sudo nginx -t
sudo systemctl reload nginx
```

### Step 6: Enable HTTPS

```bash
sudo certbot --nginx -d your-domain.com
```

### Step 7: Verify the deployment

Check the container:

```bash
docker ps
docker logs protein-knot-detector
```

Then open:

```text
https://your-domain.com
```

## 7. What Happens on the Server During a Real Request

For a short sequence:

1. Browser sends `POST /`
2. Django writes a temporary PDB file such as `web_result_ab12cd34.pdb`
3. Django calls the external ESMFold API
4. Django stores the returned structure in that file
5. `topoly.alexander(...)` runs on the file
6. Django reads the file contents back into the template for visualization
7. Django tries to delete the file
8. HTML is returned to the browser with the embedded structure data

For a long sequence:

1. Browser sends `POST /`
2. Django splits the sequence into overlapping chunks
3. Each chunk is sent to ESMFold separately
4. Each chunk is written to a fragment file like `fragment_1_400.pdb`
5. `topoly` runs for each fragment
6. The first detected knot is reported back in the rendered HTML
7. Fragment files remain on disk unless you clean them manually

## 8. What You Should Change Before Public Deployment

Before putting this on a real public server, make these code changes:

1. Move `SECRET_KEY`, `DEBUG`, and `ALLOWED_HOSTS` into environment variables
2. Set `DEBUG = False`
3. Restrict `ALLOWED_HOSTS` to your real domain
4. Replace SQLite if you need reliability beyond a simple demo
5. Store temporary PDB files in a per-request temporary directory
6. Delete fragment files after analysis
7. Stop using shared fixed names like `fragment_201_600.pdb` across requests
8. Add request logging and application error logging

## 9. Best Next Architecture When You Want to Scale

If you want multiple users and reliable long-running jobs, move to this architecture:

- Django or FastAPI API service
- Celery worker
- Redis queue
- PostgreSQL database
- object storage for generated PDB files
- a separate frontend if you want a richer UI

The improved flow becomes:

1. User submits a sequence
2. Web app creates a background job
3. API immediately returns a job ID
4. Celery worker calls ESMFold and runs `topoly`
5. Result is stored in database or object storage
6. Frontend polls job status and displays the completed result

This is better because:

- web requests return quickly
- server timeouts stop being the main bottleneck
- multiple analyses can queue safely
- workers can scale separately from the web server

## 10. Recommended Hosting Choices

### Option A: Single VPS

Use this if you want full control and the lowest cost.

Good for:

- demos
- private access
- small research groups

Suggested providers:

- DigitalOcean Droplet
- AWS EC2
- Hetzner Cloud
- Linode/Akamai

### Option B: Docker-based PaaS

Use this if you want simpler ops and are fine with platform limits.

Good for:

- easier deployment
- less server administration

Examples:

- Render
- Railway

Important note: because this app performs long synchronous work, PaaS timeouts can still be a problem unless you add a background worker architecture.

## 11. Bottom Line

You can deploy this project today as a small single-server Dockerized Django app behind Nginx.

That deployment will work like this:

- user opens the Django site
- user submits a protein sequence
- server calls ESMFold
- server runs `topoly`
- server returns an HTML page with the knot result and viewer

But if you want this to serve real public traffic, you should next refactor:

- production settings
- temporary file handling
- background job processing
- concurrency handling

## 12. Suggested Next Work in This Repo

If you want, the next useful code changes are:

1. production-ready Django settings using environment variables
2. `docker-compose.yml` for easier server deployment
3. Nginx config file committed in the repo
4. temp-file cleanup fixes in `detector/views.py`
5. Celery + Redis integration for long-running analyses

