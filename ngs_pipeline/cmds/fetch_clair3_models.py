# SPDX-FileCopyrightText: 2026 Hidayat Trimarsanto <trimarsanto@gmail.com>
# SPDX-License-Identifier: MIT

from __future__ import annotations

__copyright__ = "(C) 2026 Hidayat Trimarsanto <trimarsanto@gmail.com>"
__author__ = "trimarsanto@gmail.com"
__license__ = "MIT"

# to improve the responsiveness during bash autocomplete, do not import heavy
# modules (such as numpy, pandas, etc) here, but instead import them within the
# functions that requires respective heavy modules

import sys
import os
from ngs_pipeline import cerr, cexit, arg_parser, check_VVG_BASEDIR

# this command will download Clair3 models from the official release page,
# and place them in the expected location under $VVG_BASEDIR/opt/clair3-models/
# for use by the Clair3 variant caller

from html.parser import HTMLParser
import urllib.parse

source_urls = [
    "https://www.bio8.cs.hku.hk/clair3/clair3_models_rerio_pytorch/",
    "https://www.bio8.cs.hku.hk/clair3/clair3_models_pytorch/",
]


def init_argparser():
    p = arg_parser(desc="download Clair3 models (PyTorch-based)")
    p.add_argument(
        "--model",
        required=True,
        help="model name to download (e.g. 'clair3-ont', 'clair3-pacbio', etc)",
    )
    return p


class IndexLinkParser(HTMLParser):
    """A clean, built-in parser to harvest links from HTML anchor tags."""

    def __init__(self):
        super().__init__()
        self.links = []

    def handle_starttag(self, tag, attrs):
        if tag == "a":
            for attr, value in attrs:
                if attr == "href":
                    self.links.append(value)


def web_path_exists(url: str, verify_ssl: bool = True) -> bool:
    """Checks if a remote file or directory URL exists using an HTTP HEAD request."""

    import requests

    headers = {"User-Agent": "Mozilla/5.0"}
    try:
        response = requests.head(url, headers=headers, verify=verify_ssl, timeout=10)
        print(f"Checking {url} - Status Code: {response.status_code}")
        # Treat any 2xx/3xx response as existing
        return 200 <= response.status_code < 400
    except requests.RequestException as e:
        print(f"Connection check failed for {url}: {e}")
        return False


def download_web_directory(
    url: str, target_dir: str, verify_ssl: bool = True, base_path: str | None = None
) -> None:
    """Recursively downloads a web directory structure using requests and HTMLParser."""

    import requests

    # local imports to avoid slowing down CLI auto-complete

    if not url.endswith("/"):
        url += "/"

    os.makedirs(target_dir, exist_ok=True)
    headers = {"User-Agent": "Mozilla/5.0"}

    try:
        response = requests.get(url, headers=headers, verify=verify_ssl)
        response.raise_for_status()
        html_text = response.text
    except requests.RequestException as e:
        print(f"Skipping directory {url}: {e}")
        return

    # Parse HTML using standard library HTMLParser
    parser = IndexLinkParser()
    parser.feed(html_text)

    parsed_url = urllib.parse.urlparse(url)
    if base_path is None:
        base_path = parsed_url.path

    for href in parser.links:
        # Skip parent navigation, server parameters, and root links
        if not href or href in ["../", "/"] or href.startswith("?"):
            continue

        # Skip fully-qualified links to other hosts
        parsed_href = urllib.parse.urlparse(href)
        if parsed_href.netloc and parsed_href.netloc != parsed_url.netloc:
            continue

        # if pointed to parent directory, skip to avoid infinite recursion
        if parsed_url.path.startswith(href):
            continue

        full_url = urllib.parse.urljoin(url, href)
        parsed_full = urllib.parse.urlparse(full_url)
        decoded_name = urllib.parse.unquote(href)

        # Compute local path relative to the original base_path so we don't
        # recreate the full upstream URL path; only create paths under the
        # model root (the last directory in the initial URL).
        rel_path = os.path.relpath(parsed_full.path, start=base_path)
        if rel_path == ".":
            rel_path = ""
        # If the resolved path escapes the base, fall back to using only the
        # basename to avoid creating parent URL directories locally.
        if rel_path.startswith(".."):
            rel_path = os.path.basename(parsed_full.path)

        local_path = os.path.normpath(os.path.join(target_dir, rel_path.lstrip("/")))

        if href.endswith("/"):
            print(f"Directory -> {decoded_name}")
            os.makedirs(local_path, exist_ok=True)
            download_web_directory(full_url, target_dir, verify_ssl, base_path)
        else:
            print(f"Downloading -> {decoded_name}")
            try:
                os.makedirs(os.path.dirname(local_path), exist_ok=True)

                # Stream the download data chunks
                with requests.get(
                    full_url,
                    headers=headers,
                    verify=verify_ssl,
                    stream=True,
                ) as r:
                    r.raise_for_status()
                    with open(local_path, "wb") as f:
                        for chunk in r.iter_content(chunk_size=8192):
                            if chunk:
                                f.write(chunk)
            except requests.RequestException as e:
                print(f"Failed to download {decoded_name}: {e}")


def list_available_clair3_models(verify_ssl: bool = True) -> list[str]:
    """Return a sorted list of available Clair3 model names from all source URLs."""

    import requests

    headers = {"User-Agent": "Mozilla/5.0"}
    models: set[str] = set()

    for source_url in source_urls:
        try:
            response = requests.get(
                source_url, headers=headers, verify=verify_ssl, timeout=20
            )
            response.raise_for_status()
        except requests.RequestException as e:
            print(
                f"Warning: failed to read model index {source_url}: {e}",
                file=sys.stderr,
            )
            continue

        parser = IndexLinkParser()
        parser.feed(response.text)

        for href in parser.links:
            if not href or href in ["../", "/"] or href.startswith("?"):
                continue

            parsed_href = urllib.parse.urlparse(href)

            # skip external hosts
            if parsed_href.netloc:
                continue

            name = urllib.parse.unquote(href).strip("/")

            # keep only top-level directory-like entries
            if not name or "/" in name:
                continue

            models.add(name)
    return sorted(models)


def fetch_clair3_models(args):

    import pathlib

    vvg_basedir = check_VVG_BASEDIR()

    base_dir = pathlib.Path(vvg_basedir) / "opt" / "clair3-models"
    model = args.model.strip().replace(".", "").replace("@", "_")
    model_dir = base_dir / model

    for source_url in source_urls:
        model_url = source_url + model
        if web_path_exists(model_url):
            print(f"Fetching Clair3 models from {model_url}...")
            download_web_directory(model_url, model_dir.as_posix())
            print(f"Finished downloading Clair3 models to {model_dir}")
            break
    else:
        print(f"Error: Model '{model}' not found in any source URLs.")
        sys.exit(1)


def check_and_fetch_clair3_model(args):

    # lock the directory to prevent concurrent downloads using flufl-lock
    from flufl.lock import Lock
    import pathlib

    vvg_basedir = check_VVG_BASEDIR()
    base_dir = pathlib.Path(vvg_basedir) / "opt" / "clair3-models"
    model = args.model.strip().replace(".", "").replace("@", "_")
    model_dir = base_dir / model

    lock_file = base_dir / ".fetch_clair3_model.lock"
    with Lock(lock_file, lifetime=60 * 60, acquire_timeout=10):
        if model_dir.exists() and any(model_dir.iterdir()):
            print(f"Model '{model}' already exists at {model_dir}. Skipping download.")
        else:
            fetch_clair3_models(args)


def main(args):
    fetch_clair3_models(args)


# TODO:
# - use flufl-lock to check_and_fetch to prevent multiple processes downloading the same model
#   simultaneously, which can lead to corrupted downloads.
#   This is especially important in share environments where multiple jobs may attempt to fetch
#   the same model concurrently.

# EOF
