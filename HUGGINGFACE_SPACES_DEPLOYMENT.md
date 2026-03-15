# Deploying This Project on Hugging Face Spaces

This guide is specific to the current `ProteinKnotDetector` repository.

It covers:

1. How the app works on Hugging Face today
2. What errors or deployment blockers exist for Spaces
3. What I changed in this repo to make Spaces deployment safer
4. The exact steps to deploy it as a Docker Space
5. The limits you should keep in mind before putting real users on it

## 1. Which Hugging Face Space Type You Need

Use a **Docker Space**.

Do not use:

- Static Space
- Gradio Space
- Streamlit Space

This repository is a Django application with a `Dockerfile`, so Docker Space is the correct deployment target.

## 2. How It Will Work on Hugging Face

Once deployed to a Docker Space, the runtime flow is:

1. Hugging Face builds the container from your `Dockerfile`.
2. The platform starts the container and exposes the app on the port declared in the Space metadata.
3. Gunicorn starts Django on port `7860`.
4. A user opens the Space URL.
5. Django renders the page at `/`.
6. The user submits a FASTA sequence.
7. Django sanitizes the FASTA input and removes headers/non-amino-acid characters.
8. Django sends the cleaned sequence to the ESMFold API over HTTPS.
9. The returned PDB structure is written to a temporary directory inside the container.
10. `topoly` analyzes that structure and estimates knot probability.
11. Django renders the result page and, for short sequences, embeds the PDB content into the 3D viewer.
12. Temporary files are removed automatically when the request finishes.

In practical terms:

- Hugging Face hosts the Django container
- ESMFold still does the structure prediction remotely
- your Space mainly performs request handling, PDB parsing, and `topoly` analysis
- the app is still synchronous, so the browser waits for the analysis to finish

## 3. Errors and Risks I Found for Hugging Face

These were the main repo-specific problems for Spaces:

### A. Writable file issue in Docker Spaces

The app writes `.pdb` files during analysis.

In the original repo:

- files were written into the application directory
- the Docker image copied files in as root
- Hugging Face Docker Spaces run with user ID `1000`

That combination can cause runtime write failures even when the Space starts correctly.

### B. Missing Space metadata

Docker Spaces need a `README.md` with YAML front matter such as:

- `sdk: docker`
- `app_port: 7860`

This repo did not have that file.

### C. Large Docker build context

The repository contains many generated/sample `.pdb` files and `db.sqlite3`.

Without a `.dockerignore`, those files get copied into the image, which makes Spaces builds slower and larger than necessary.

### D. Production Django settings were not proxy-aware

The original settings had:

- `DEBUG = True`
- `ALLOWED_HOSTS = ['*']`
- a hard-coded `SECRET_KEY`

That is not appropriate for a public Space.

### E. FASTA input handling was too loose

The page says FASTA input is supported, but the Django view previously just stripped whitespace. A header like `>protein_1` could leak invalid characters into the request sent to ESMFold.

### F. Long-request behavior is still a real limitation

This is not a bug in the repo, but it matters on Spaces:

- short sequences should work
- long sequences trigger sliding-window processing
- sliding-window processing makes multiple ESMFold calls plus multiple `topoly` runs
- the request remains synchronous the whole time

That means users will wait on one open browser request until the full job finishes.

## 4. Changes Already Made in This Repo

I updated the repository for safer Docker Space deployment:

- `Dockerfile`
  - creates a user with UID `1000`
  - uses a writable home directory
  - copies project files with the correct ownership
  - exposes port `7860`
- `detector/views.py`
  - sanitizes FASTA input
  - writes temporary PDB files into a request-scoped temporary directory
  - cleans up fragment files automatically
  - avoids file collisions across users
  - adds a request timeout for ESMFold calls
- `knot_web/settings.py`
  - supports environment variables for `SECRET_KEY`, `DEBUG`, and hosts
  - adds reverse-proxy-friendly settings for Spaces
- `.dockerignore`
  - excludes `.pdb`, SQLite, cache, and git metadata from the image build
- `requirements.txt`
  - pins Django to the Python 3.11-compatible release line used by the container
- `README.md`
  - adds the Docker Space metadata Hugging Face expects

## 5. Step-by-Step Deployment

### Step 1: Make sure the repo contains these files

Required for this project:

- `Dockerfile`
- `README.md` with Docker Space YAML front matter
- `requirements.txt`
- Django project files

Already added:

- `README.md`
- `.dockerignore`

### Step 2: Create a new Space

In Hugging Face:

1. Go to `https://huggingface.co/new-space`
2. Choose your account or organization
3. Give the Space a name
4. Select **Docker** as the SDK
5. Choose visibility: public or private
6. Create the Space

### Step 3: Add your code to the Space repo

You can do this either through Git or through the web UI.

Using Git:

```bash
git remote add hf https://huggingface.co/spaces/<your-username>/<your-space-name>
git push hf main
```

If your branch is not `main`, push the branch you are actually using.

### Step 4: Set Space variables and secrets

Open the Space settings and add:

Variables:

- `DJANGO_DEBUG=False`
- `DJANGO_ALLOWED_HOSTS=.hf.space,huggingface.co`
- `DJANGO_CSRF_TRUSTED_ORIGINS=https://*.hf.space,https://huggingface.co`

Secrets:

- `DJANGO_SECRET_KEY=<a long random secret>`

If you later add custom domains, databases, or storage, add those here too.

### Step 5: Let Hugging Face build the image

After the push:

1. Open the Space
2. Go to the build logs
3. Wait for the Docker image to build
4. Confirm Gunicorn starts successfully

You want to see the container come up without import errors and without permission errors.

### Step 6: Test the root page

Check that:

1. `/` loads successfully
2. the HTML form appears
3. the page renders without Django host or CSRF errors

### Step 7: Test a short sequence first

Start with a small sequence under 400 residues.

Expected behavior:

1. form submits
2. Space calls the ESMFold API
3. result page returns
4. knot report appears
5. 3D viewer renders the returned structure

### Step 8: Test a long sequence second

Now test a sequence over 400 residues.

Expected behavior:

1. the app switches to sliding-window mode
2. multiple chunk analyses run
3. the page returns a knot result if one is detected

This will be much slower than the short-sequence path.

## 6. Recommended Space Hardware

Start with CPU hardware.

Why:

- the app itself is a Django server plus CPU-bound `topoly`
- structure prediction is offloaded to the external ESMFold API
- you do not need a GPU to launch the current version

Start small, then upgrade only if request times or memory usage justify it.

## 7. What to Watch Closely on Spaces

### Request duration

This app still performs the full analysis inside the web request.

Keep in mind:

- long sequences will feel slow
- users may think the page is stuck
- long-running requests are the main risk for a public Space

### Outbound network dependency

The Space must be able to reach:

- `https://api.esmatlas.com/foldSequence/v1/pdb/`

If that API is slow, unavailable, or rate-limited, your Space will fail even if your own code is correct.

### Non-persistent local files

Do not rely on local disk for long-term storage in Spaces.

For this app that means:

- temporary PDB files are fine
- SQLite is not a good permanent production database
- saved results should eventually move to external storage or an external database

### Cold starts and sleeping

If you use free or lower-cost hardware, the Space may sleep when idle and take longer to start on the next request.

That is normal platform behavior and it matters because this app already has non-trivial request time.

### Public internet exposure

Even on Hugging Face, this is still a public web deployment if the Space is public.

Treat it like production:

- keep `DEBUG` off
- keep secrets in Space secrets, not in git
- tighten hosts/origins

## 8. Troubleshooting

### Build fails during `pip install`

Check:

- whether `topoly` resolved correctly for the selected Python version
- whether a dependency version changed upstream

If this becomes unstable, pin versions in `requirements.txt`.

### Space starts, but submit fails

Most likely causes:

- ESMFold API request failed
- CSRF or host configuration is wrong
- the sequence contains invalid FASTA content

Check the runtime logs first.

### Permission denied while writing files

That usually means the container filesystem path is not writable by the runtime user.

This repo was updated to avoid that by:

- running with UID `1000`
- writing analysis files to a temporary directory

### Very slow responses for long proteins

That is expected with the current synchronous design.

The real fix is architectural:

- queue the job
- process it in a worker
- poll for result status

## 9. Best Next Improvement After the First Hugging Face Deploy

If the first version works and you want it to be more reliable, the next change should be background processing.

A better architecture would be:

1. user submits sequence
2. Space creates a job ID
3. analysis runs asynchronously
4. frontend polls status
5. result appears when complete

That change matters more than UI polish because your bottleneck is long-running synchronous compute.

## 10. Bottom Line

Yes, this project can run on Hugging Face Spaces today as a **Docker Space**.

The current best path is:

1. push this repo to a Docker Space
2. set the Django environment variables in Space settings
3. test short sequences first
4. expect long sequences to be slow
5. move to async jobs if you want reliability at larger scale

## 11. Useful Official References

- Spaces overview: `https://huggingface.co/docs/hub/spaces-overview`
- Docker Spaces: `https://huggingface.co/docs/hub/spaces-sdks-docker`
- Space configuration reference: `https://huggingface.co/docs/hub/spaces-config-reference`
