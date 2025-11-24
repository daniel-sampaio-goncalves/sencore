# sencore
*A user‑friendly web application to query and analyze sequencing data from cellular senescence models.*


![Python](https://img.shields.io/badge/python-3.13.5-blue)
![Streamlit](https://img.shields.io/badge/streamlit-1.51.0-blue)
![docker](https://img.shields.io/badge/docker-29.0.1-blue)
![License](https://img.shields.io/badge/license-Apache%202.0-blue)
![Version](https://img.shields.io/badge/version-1.0.0-blue)

---

## Short Description
*sencore* is a lightweight Streamlit web app that lets biologists explore a curated database of 44 bulk
RNA‑sequencing experiments comparing senescent cell models to non‑senescent controls.
With an intuitive interface and Docker support, you can start querying and visualising data in minutes.
*sencore* aggregates these datasets into a single SQLite database and exposes them through an interactive web UI.
Users can filter experiments, view differential expression results, and download summary tables—all without
writing code.

---

## Table of Contents
- [Key Features](#key-features)
- [Prerequisites / Dependencies](#prerequisites--dependencies)
- [Installation / Setup](#installation--setup)
- [Docker](#installation-using-docker)
- [Demo](#demo)
- [Contributing](#contributing)
- [License](#license)
- [FAQ](#faq)

---

## Key Features
- **Streamlit‑based UI**: No command line, just a browser.
- **SQLite backend**: Lightweight, zero‑configuration database.
- **Docker support**: Run the app in an isolated container.
---

## Prerequisites / Dependencies
| Component | tested version |
|-----------|-----------------|
| Python | 3.13.5 |
| Streamlit | 1.51.0
| pip | 25.1 |
| Docker (optional) | 29.0.1+ |

Python dependencies are listed in `requirements.txt`.

---

## Installation / Setup

### 1. Clone the Repository
```bash
git clone https://github.com/daniel-sampaio-goncalves/sencore.git
cd sencore
```

### 2. Install Python Dependencies
```bash
pip install -r requirements.txt
```

### 3. Download the Database
> **Note:** The database is stored at [Open Science Framework (OSF)](https://osf.io/) in a gzipped format; download it and ensure you decompress it before launching the app.
```bash
curl -L -o SENCORE_08_11_2025.db.gz "https://osf.io/download/p84gq/?view_only=9539cd4703fe42c09e8fb1cf0567cc6b"
gzip -d SENCORE_08_11_2025.db.gz
```



### 4. Run the App (Local)
```bash
streamlit run app.py \
  --server.port=8501 \
  --server.address=0.0.0.0 \
  --theme.base=dark \
  --browser.gatherUsageStats=false \
  --global.showWarningOnDirectExecution=true \
  --server.headless=true \
  --server.fileWatcherType='none' \
  --server.showEmailPrompt=false \
  --client.toolbarMode=minimal
```

Open your browser at `http://localhost:8501`.

## Installation using Docker
### 1. Run the App (Docker)
> **Note:** Docker is needed for this, see https://docs.docker.com/engine/install/, alternatively you can install the [docker desktop application](https://docs.docker.com/desktop/)

```bash
#after docker installation
git clone https://github.com/daniel-sampaio-goncalves/sencore.git
cd sencore/docker
docker compose up -d
```

Your app will be available at `http://localhost:8501`.

## DEMO

The app is also freely available at https://sencore.koyeb.app/
> **Note:** The app is hosted on a lightweight server, which can run into memory limitations when processing large data queries.

---

## Contributing
*sencore* is an open‑source project, but at the moment I’m not working on it further until maybe a later time. Feel free to fork the repository, modify the app, or use the data for your own analyses. If you discover bugs or have comments, feel free to open an issue on GitHub.

---

## License
Apache 2.0 – see the [LICENSE](LICENSE) file for details.

---

## FAQ

| Question | Answer |
|----------|--------|
| **Where can I find the raw sequencing data?** | The raw data are not bundled with sencore; the app only uses the pre‑processed differential expression tables. Most of the data is referenced in our two preprints [Senescent cells exhibit features of developmental signaling centres](https://doi.org/10.1101/2025.10.30.685565)  and [Developmental senescence profiling reveals Eda2r as a common mediator in the core senescence program](https://doi.org/10.1101/2025.10.30.685553)|
| **Is there any data privacy concern?** | The database contains only processed gene‑level statistics, no raw sequencing reads or personally identifiable information and is freely open to explore|

---
