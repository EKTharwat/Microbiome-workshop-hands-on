#!/usr/bin/env python3
"""Offline-capable mothur installer and MiSeq SOP runner (Python 3.10+)."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import urllib.request
import zipfile

ROOT = Path(__file__).resolve().parents[1]

def asset_name(system=None, machine=None):
    system = system or platform.system()
    machine = (machine or platform.machine()).lower()
    arch = {"aarch64": "arm64", "arm64": "arm64", "amd64": "x86_64", "x86_64": "x86_64"}.get(machine)
    if not arch:
        raise ValueError(f"Unsupported CPU: {machine}")
    if system == "Windows":
        if arch != "x86_64":
            raise ValueError("Windows ARM is not supported by this bundle")
        return "Mothur.win.zip"
    if system == "Darwin":
        return f"Mothur.OSX_{arch}.zip"
    if system == "Linux":
        return f"Mothur.linux_{arch}.noReadline.zip"
    raise ValueError(f"Unsupported OS: {system}")

def extract(archive, destination):
    destination.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(archive) as z:
        for member in z.infolist():
            target = (destination / member.filename).resolve()
            if not target.is_relative_to(destination.resolve()):
                raise ValueError("Unsafe archive path")
        z.extractall(destination)

def resource(name):
    manifest = json.loads((ROOT / "resources.json").read_text())
    item = next(x for x in manifest["resources"] if x["name"] == name)
    target = ROOT / "vendor" / name
    if not target.exists():
        target.parent.mkdir(exist_ok=True)
        partial = target.with_suffix(".partial")
        print(f"Downloading {name}", flush=True)
        urllib.request.urlretrieve(item["url"], partial)
        partial.replace(target)
    with target.open("rb") as handle:
        digest = hashlib.file_digest(handle, "sha256").hexdigest() if hasattr(hashlib, "file_digest") else hashlib.sha256(handle.read()).hexdigest()
    if digest != item["sha256"]:
        raise ValueError(f"Checksum mismatch: {target}; remove it and retry")
    return target

def install():
    name = asset_name()
    destination = ROOT / ".tools" / name.removesuffix(".zip")
    extract(resource(name), destination)
    executable = next(f for f in destination.rglob("mothur.exe" if os.name == "nt" else "mothur") if f.is_file())
    for f in destination.rglob("*"):
        if f.is_file() and (f.name in ("mothur", "vsearch", "uchime") or f.suffix == ".so"):
            f.chmod(f.stat().st_mode | 0o111)
    (ROOT / ".tools" / "executable.txt").write_text(str(executable.resolve()))
    print(f"Installed {executable}")
    return executable

def main():
    parser = argparse.ArgumentParser(description="Install the bundled mothur version for your OS.")
    parser.add_argument("command", choices=["install"])
    parser.parse_args()
    try:
        install()
    except (ValueError, OSError, RuntimeError, StopIteration) as exc:
        parser.exit(1, f"Error: {exc}\n")

if __name__ == "__main__":
    main()
