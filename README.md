---
title: Protein Knot Detector
sdk: docker
app_port: 7860
short_description: Django app for protein knot detection with ESMFold and Topoly.
---

# Protein Knot Detector

This repository contains a Docker-based Django Space for protein knot analysis.

The web app:

- accepts FASTA protein sequences
- sends sequences to the ESMFold API
- runs `topoly` on the returned PDB structure
- reports knot probability and knot type

For deployment details, see:

- `HUGGINGFACE_SPACES_DEPLOYMENT.md`
- `SERVER_DEPLOYMENT_GUIDE.md`
